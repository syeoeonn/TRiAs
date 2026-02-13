# src/trias/allele/extract.py
from __future__ import annotations
import sys, os, shutil, re, gzip, subprocess, tempfile, atexit, logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional

import pandas as pd

log = logging.getLogger("trias")
if not log.handlers:
    h = logging.StreamHandler(sys.stderr)
    h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    log.addHandler(h)
log.setLevel(logging.INFO)

# ------------ tiny helpers ------------

def run_command(cmd_list: List[str]) -> Optional[str]:
    """Execute a command and return stdout, or None if failed."""
    env = dict(os.environ, OMP_NUM_THREADS="1")
    r = subprocess.run(cmd_list, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        log.error("Command failed: %s\n%s", " ".join(cmd_list), r.stderr)
        return None
    return r.stdout

def reverse_complement(seq: str) -> str:
    """Return reverse complement of a DNA sequence."""
    comp = {'A':'T','T':'A','G':'C','C':'G','a':'t','t':'a','g':'c','c':'g','N':'N','n':'n'}
    return "".join(comp.get(b, 'N') for b in seq[::-1])

def parse_fasta_output(fasta_string: str) -> tuple[str, str]:
    """Parse FASTA string → (header, sequence)."""
    lines = fasta_string.strip().splitlines()
    if not lines:
        return "", ""
    header = lines[0].lstrip(">")
    seq = "".join(lines[1:])
    return header, seq

# ------------ grouping state ------------

SAMPLE_GROUPS: Dict[str, tuple[str, int]] = {}  # base sample -> (group_name, ploidy)
GROUP_ORDER: List[str] = []                      # group display order
REF_GROUP: Optional[str] = None                  # first group -> treated as REF

def load_sample_groups(group_files: List[str]) -> None:
    """
    Load sample group information from TSV files.

    Each file is a group. The first file is the *reference group*.
    File format: base_sample_name<TAB>ploidy
    """
    global SAMPLE_GROUPS, GROUP_ORDER, REF_GROUP
    SAMPLE_GROUPS = {}
    GROUP_ORDER = []
    REF_GROUP = None

    for i, path in enumerate(group_files):
        if not os.path.exists(path):
            log.error("Group file does not exist: %s", path)
            sys.exit(1)

        name = os.path.basename(path)
        if name.endswith(".tsv"):
            name = name[:-4]
        GROUP_ORDER.append(name)
        if i == 0:
            REF_GROUP = name  # first file is the reference group (must match the VCF REF)

        try:
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    s = line.strip()
                    if not s:
                        continue
                    parts = s.split("\t")
                    if len(parts) >= 2:
                        base = parts[0].strip()
                        ploidy = int(parts[1].strip())
                        SAMPLE_GROUPS[base] = (name, ploidy)
        except Exception as e:
            log.error("Failed to read %s: %s", path, e)
            sys.exit(1)

    log.info("Loaded groups: %s", GROUP_ORDER)
    log.info("Reference group: %s", REF_GROUP)
    log.info("Total samples: %d", len(SAMPLE_GROUPS))

def classify_sample(sample_name: str) -> str:
    """Return the group name for a sample (Unknown if not matched)."""
    base = sample_name.rsplit(".", 1)[0] if sample_name.endswith((".1",".2")) else sample_name
    if base in SAMPLE_GROUPS:
        return SAMPLE_GROUPS[base][0]
    log.warning("Sample %s (base: %s) not found in any group.", sample_name, base)
    return "Unknown"

def get_haploid_samples() -> List[str]:
    return [base for base, (_, ploidy) in SAMPLE_GROUPS.items() if ploidy == 1]

# ------------ VCF helpers ------------

def split_diploid_vcf(in_vcf: str, out_vcf: str, skip_samples: Optional[List[str]] = None) -> None:
    """
    Split a diploid VCF into a haploid-like VCF where each diploid sample becomes sample.1 and sample.2.
    Samples in skip_samples are kept as single column (e.g., haploid).
    """
    if skip_samples is None:
        skip_samples = get_haploid_samples()
    tmp_uncompressed = out_vcf.replace(".gz", "")

    opener_in = gzip.open if in_vcf.endswith(".gz") else open
    with opener_in(in_vcf, "rt") as fin, open(tmp_uncompressed, "w", encoding="utf-8") as fout:
        sample_names: List[str] = []
        for line in fin:
            line = line.rstrip("\n")
            if line.startswith("##"):
                fout.write(line + "\n")
                continue
            if line.startswith("#CHROM"):
                parts = line.split("\t")
                fixed = parts[:9]
                original = parts[9:]
                new_cols: List[str] = []
                for s in original:
                    if s in skip_samples:
                        new_cols.append(s)
                    else:
                        new_cols.append(s + ".1")
                        new_cols.append(s + ".2")
                sample_names = original
                fout.write("\t".join(fixed + new_cols) + "\n")
                continue

            parts = line.split("\t")
            fixed = parts[:9]
            smp_cols = parts[9:]
            new_vals: List[str] = []
            for i, sname in enumerate(sample_names):
                val = smp_cols[i]
                if sname in skip_samples:
                    new_vals.append(val)
                else:
                    sub = val.split(":")
                    gt_full = sub[0] if sub else "."
                    sep = "|" if "|" in gt_full else ("/" if "/" in gt_full else None)
                    if sep:
                        alleles = gt_full.split(sep)
                        a1, a2 = (alleles if len(alleles) == 2 else ("0", "0"))
                    else:
                        a1, a2 = gt_full, "0"
                    new_vals.append(a1)
                    new_vals.append(a2)

            fout.write("\t".join(fixed + new_vals) + "\n")

    subprocess.run(["bgzip", "-f", tmp_uncompressed], check=True)
    subprocess.run(["tabix", "-p", "vcf", out_vcf], check=True)

# ------------ core region processing ------------

def process_one_region(row_dict: dict, big_vcf_path: str, base_tmp_dir: str, idx: int) -> List[dict]:
    """Process one TR interval and return a list of result dicts (rows)."""
    chrom = row_dict.get("CHROM", "")
    original_start = int(row_dict.get("START", 0))
    original_end   = int(row_dict.get("END", 0))
    ref_allele = row_dict.get("TR_allele", "")
    motif = row_dict.get("MOTIF", "")

    is_neg_strand = (original_start > original_end)
    start, end = sorted([original_start, original_end])

    region_start_1based = start + 1
    region_end_1based = end
    region_str = f"{chrom}:{region_start_1based}-{region_end_1based}"

    with tempfile.TemporaryDirectory(dir=base_tmp_dir) as tmpdir:
        region_vcf_gz = os.path.join(tmpdir, "region.vcf.gz")
        if run_command(["bcftools", "view", "--threads", "1", "-r", region_str, "-Oz", "-o", region_vcf_gz, big_vcf_path]) is None:
            return []
        if run_command(["tabix", "-p", "vcf", region_vcf_gz]) is None:
            return []

        variant_lines = run_command(["bcftools", "view", "--threads", "1", "-H", region_vcf_gz])
        if not variant_lines or not variant_lines.strip():
            log.info("Region %s: no variants, skipping.", region_str)
            return []

        splitted_vcf = os.path.join(tmpdir, "region.split.vcf.gz")
        split_diploid_vcf(region_vcf_gz, splitted_vcf)

        orig_samples_str = run_command(["bcftools", "query", "-l", big_vcf_path]) or ""
        original_samples = [s.strip() for s in orig_samples_str.strip().splitlines() if s.strip()]
        haploid_samples = get_haploid_samples()

        expected_samples: List[str] = []
        for s in original_samples:
            if s in haploid_samples:
                expected_samples.append(s)
            else:
                expected_samples.append(s + ".1")
                expected_samples.append(s + ".2")

        actual_samples_str = run_command(["bcftools", "query", "-l", splitted_vcf]) or ""
        actual_samples = [s.strip() for s in actual_samples_str.strip().splitlines() if s.strip()]

        bcftools_ref_allele = reverse_complement(ref_allele) if is_neg_strand else ref_allele
        temp_fa = os.path.join(tmpdir, "TEMP.fa")
        with open(temp_fa, "w", encoding="utf-8") as fw:
            fw.write(f">{chrom}:{region_start_1based}-{region_end_1based}\n{bcftools_ref_allele}\n")

        sample_signatures = {s: [] for s in expected_samples}
        missing_samples = set()

        gt_lines = run_command(["bcftools", "query", "-f", "[%GT\\t]\\n", splitted_vcf]) or ""
        if gt_lines:
            for line in gt_lines.strip().splitlines():
                cols = line.split("\t")
                # collect GT only for actual samples
                for i, s in enumerate(actual_samples):
                    if i < len(cols):
                        gt = cols[i]
                        if gt == ".":
                            gt = "0"   # treat missing as REF
                            missing_samples.add(s)
                        sample_signatures[s].append(gt)
                # expected-but-missing samples → default 0
                for s in expected_samples:
                    if s not in actual_samples:
                        sample_signatures[s].append("0")

        normalized_signatures: Dict[str, str] = {}
        for s in expected_samples:
            if not sample_signatures[s]:
                sample_signatures[s] = ["0"]
            normalized_signatures[s] = "|".join(sample_signatures[s])

        ref_variant_count = max(1, len(next(iter(sample_signatures.values()))))
        ref_signature = "|".join(["0"] * ref_variant_count)
        if REF_GROUP is None:
            log.error("REF_GROUP is not set. First group file must be the reference group.")
            return []
        normalized_signatures[REF_GROUP] = ref_signature
        sample_signatures[REF_GROUP] = ["0"] * ref_variant_count

        # group by signature
        sig_groups: Dict[str, List[str]] = defaultdict(list)
        for s, sig in normalized_signatures.items():
            sig_groups[sig].append(s)

        allele_counts: Dict[str, Dict[str, int]] = {}
        allele_samples: Dict[str, Dict[str, List[str]]] = {}
        grp_idx = 0

        for sig, group_samps in sig_groups.items():
            grp_idx += 1
            # representative sample: REF group if present, else the first sample in the signature group
            rep_sample = REF_GROUP if REF_GROUP in group_samps else group_samps[0]
            is_all_ref = all(gt == "0" for gt in sig.split("|"))

            if is_all_ref or rep_sample == REF_GROUP:
                final_seq = reverse_complement(bcftools_ref_allele) if is_neg_strand else bcftools_ref_allele
            else:
                try:
                    group_name = f"G{grp_idx}"
                    raw_vcf = os.path.join(tmpdir, f"group_{grp_idx}_raw.vcf")
                    env = dict(os.environ, OMP_NUM_THREADS="1")
                    proc = subprocess.run(
                        ["bcftools", "view", "--threads", "1", "-s", rep_sample, "-o", "-", splitted_vcf],
                        capture_output=True, text=True, env=env
                    )
                    if proc.returncode != 0:
                        raise Exception(f"bcftools view failed: {proc.stderr}")
                    with open(raw_vcf, "w", encoding="utf-8") as fw:
                        fw.write(proc.stdout)

                    reheader_txt = os.path.join(tmpdir, f"group_{grp_idx}_reheader.txt")
                    with open(reheader_txt, "w", encoding="utf-8") as fh:
                        fh.write(f"{rep_sample}\t{group_name}\n")

                    reheader_vcf = os.path.join(tmpdir, f"group_{grp_idx}_reheader.vcf")
                    if run_command(["bcftools", "reheader", "--threads", "1", "-s", reheader_txt, "-o", reheader_vcf, raw_vcf]) is None:
                        raise Exception("bcftools reheader failed")
                    if run_command(["bgzip", "-f", reheader_vcf]) is None:
                        raise Exception("bgzip failed")
                    if run_command(["tabix", "-p", "vcf", reheader_vcf + ".gz"]) is None:
                        raise Exception("tabix failed")

                    group_vcf = reheader_vcf + ".gz"
                    cmd_cons = f"cat {temp_fa} | bcftools consensus -s {group_name} {group_vcf}"
                    r = subprocess.run(cmd_cons, shell=True, text=True, capture_output=True, env=env)
                    if r.returncode != 0:
                        raise Exception(f"bcftools consensus failed: {r.stderr}")
                    if not r.stdout:
                        raise Exception("bcftools consensus produced no output")

                    _, plus_strand_seq = parse_fasta_output(r.stdout)
                    final_seq = reverse_complement(plus_strand_seq) if is_neg_strand else plus_strand_seq

                except Exception as e:
                    log.error("Group %s processing failed: %s", grp_idx, e)
                    return []

            if final_seq not in allele_counts:
                allele_counts[final_seq] = {g: 0 for g in GROUP_ORDER}
                allele_samples[final_seq] = {g: [] for g in GROUP_ORDER}

            for s_ in group_samps:
                # mark missing genotype samples with +...+ for visibility
                display_name = f"+{s_}+" if s_ in missing_samples else s_
                gname = classify_sample(s_)
                if gname in allele_counts[final_seq]:
                    allele_counts[final_seq][gname] += 1
                    allele_samples[final_seq][gname].append(display_name)

        total_expected = len(expected_samples) + 1  # +1 for REF
        total_found = sum(sum(cnts.values()) for cnts in allele_counts.values())
        if total_found != total_expected:
            log.error("Region %s: sample count mismatch (expected=%d, found=%d)", region_str, total_expected, total_found)
            return []

        out_rows: List[dict] = []
        for hap_seq, group_dict in allele_counts.items():
            # flatten sample names following GROUP_ORDER
            all_samples: List[str] = []
            for g in GROUP_ORDER:
                all_samples.extend(sorted(allele_samples[hap_seq][g]))
            sample_str = ",".join(all_samples)

            row = {
                "__idx": idx,
                "CHROM": chrom,
                "START": original_start,
                "END":   original_end,
                "MOTIF": motif,
                "TR_allele": hap_seq,
                "Samples": sample_str
            }
            for g in GROUP_ORDER:
                row[g] = group_dict[g]
            out_rows.append(row)

        df_sub = pd.DataFrame(out_rows)
        if not df_sub.empty:
            df_sub["__len"] = df_sub["TR_allele"].apply(len)
            df_sub.sort_values(by="__len", inplace=True)
            df_sub.drop(columns="__len", inplace=True)

        return df_sub.to_dict(orient="records")

# ------------ public API (called by CLI) ------------

def vcf_alleles(
    vcf: str,
    input_tsv: str,
    output_tsv: str,
    group_files: List[str],
    processes: int = 32,
    chunk_size: int = 10_000,
    tmpdir: str = "./tmp"
) -> None:
    """
    Extract population tandem-repeat alleles from a VCF over input TR intervals.

    IMPORTANT: The **first group file** is treated as the *reference group* and must correspond
    to the VCF REF allele. The representative sequence for that group will be the REF sequence.

    Parameters
    ----------
    vcf : str
        Indexed VCF(.gz) with bcftools-ready index (.tbi/.csi).
    input_tsv : str
        TSV with at least columns: CHROM, START, END, MOTIF, TR_allele.
        START/END are interpreted in BED-like 0-based (END-exclusive) coordinates.
    output_tsv : str
        Output TSV path.
    group_files : List[str]
        List of TSV files, each defining a group. The **first file is REF**.
        Format per line: base_sample_name<TAB>ploidy
    processes : int
        Number of worker processes.
    chunk_size : int
        Number of intervals per chunk for parallel processing.
    tmpdir : str
        Base temporary directory (will be cleaned up at exit).
    """
    Path(tmpdir).mkdir(parents=True, exist_ok=True)
    atexit.register(lambda: shutil.rmtree(tmpdir, ignore_errors=True))

    load_sample_groups(group_files)

    df = pd.read_csv(input_tsv, sep="\t")
    df["START"] = df["START"].astype(int)
    df["END"]   = df["END"].astype(int)

    all_rows: List[dict] = []

    n = len(df)
    n_chunks = (n + chunk_size - 1) // chunk_size
    log.info("Total intervals: %d (chunk_size=%d, chunks=%d, processes=%d)", n, chunk_size, n_chunks, processes)

    for ci in range(n_chunks):
        now_str = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log.info("[%s] Processing chunk %d/%d ...", now_str, ci+1, n_chunks)
        df_chunk = df.iloc[ci*chunk_size : min((ci+1)*chunk_size, n)]

        with ProcessPoolExecutor(max_workers=processes) as ex:
            futs = []
            for idx_row, row in df_chunk.iterrows():
                futs.append(ex.submit(process_one_region, row.to_dict(), vcf, tmpdir, int(idx_row)))
            for fut in as_completed(futs):
                sub = fut.result() or []
                all_rows.extend(sub)

    df_out = pd.DataFrame(all_rows)
    if df_out.empty:
        log.warning("No results produced.")
        # still write an empty file with header for reproducibility?
        df_out = pd.DataFrame(columns=(GROUP_ORDER + ["CHROM","START","END","MOTIF","TR_allele","Samples"]))

    if "__idx" in df_out.columns:
        df_out.sort_values(by="__idx", kind="stable", inplace=True)
        df_out.drop(columns="__idx", inplace=True)

    col_order = GROUP_ORDER + ["CHROM","START","END","MOTIF","TR_allele","Samples"]
    col_order = [c for c in col_order if c in df_out.columns]
    df_out = df_out[col_order]
    outp = Path(output_tsv).expanduser()
    if str(outp.parent) not in ("", "."):
        outp.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(str(outp), sep="\t", index=False)
    log.info("Saved: %s", output_tsv)

# ------------ CLI adapter ------------

def vcf_alleles_cmd(args) -> None:
    """
    CLI entry for `trias vcf-alleles`.
    """
    if not args.groups or len(args.groups) == 0:
        log.error("At least one group file must be provided. The first one is the REF group.")
        sys.exit(1)
    vcf_alleles(
        vcf=args.vcf,
        input_tsv=args.input,
        output_tsv=args.output,
        group_files=args.groups,
        processes=args.processes,
        chunk_size=args.chunk_size,
        tmpdir=args.tmpdir
    )
