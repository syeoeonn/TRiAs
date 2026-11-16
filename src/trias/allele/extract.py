# src/trias/allele/extract.py
from __future__ import annotations
import sys, os, shutil, re, gzip, subprocess, tempfile, atexit, logging
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional

import pandas as pd

from ..merge.manifest import (
    STATUS_BCFTOOLS_FAILED, STATUS_CONSENSUS_FAILED, STATUS_COUNT_MISMATCH,
    STATUS_NO_REF_GROUP, STATUS_NO_VARIANTS, STATUS_OK,
    build_group_payload, file_stat, write_loci_status, write_manifest,
)

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

def check_vcf_readable(vcf: str) -> None:
    """Exit if the VCF or its index is missing, or its header cannot be read.

    Without this, a bad path makes every locus fail one by one as
    bcftools_failed, and the run still ends with exit code 0.
    """
    if not os.path.isfile(vcf):
        log.error("VCF does not exist: %s", vcf)
        sys.exit(1)
    if not (os.path.isfile(vcf + ".tbi") or os.path.isfile(vcf + ".csi")):
        log.error("VCF index (.tbi/.csi) does not exist: %s", vcf)
        sys.exit(1)
    if run_command(["bcftools", "view", "-h", vcf]) is None:
        log.error("Cannot read VCF header: %s", vcf)
        sys.exit(1)

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

# ------------ haplotype building ------------

# Insertions this close to a locus edge, made mostly of the locus motif, are
# counted as part of the locus (see apply_haplotype).
EDGE_PAD = 20


class HaplotypeError(Exception):
    """A haplotype could not be built from the records it carries."""


class UnresolvedHaplotype(HaplotypeError):
    """The locus cannot be cut out of a record that runs past its edge."""


def _locate_flank(seq: str, flank: str, from_right: bool) -> Optional[int]:
    """Find a reference flank in an ALT sequence; return the cut point or None.

    Tries the full flank, then shorter pieces of it (closest to the locus first),
    and accepts only a unique hit. The cut point is where the locus starts
    (left flank) or ends (right flank) inside `seq`.
    """
    seq_u = seq.upper()
    for n in (len(flank), 16, 12, 8):
        if n > len(flank) or n < 8:
            continue
        piece = (flank[-n:] if not from_right else flank[:n]).upper()
        hits = [m.start() for m in re.finditer(f"(?={re.escape(piece)})", seq_u)]
        if len(hits) == 1:
            return hits[0] + n if not from_right else hits[0]
    return None


def motif_fraction(seq: str, motifs: List[str]) -> float:
    """Fraction of `seq` covered by copies of the motifs (any rotation, either strand)."""
    units = set()
    for m in motifs:
        m = m.strip().upper()
        if not m:
            continue
        for x in (m, reverse_complement(m)):
            units.update(x[i:] + x[:i] for i in range(len(x)))
    if not seq or not units:
        return 0.0
    seq = seq.upper()
    lens = sorted({len(u) for u in units}, reverse=True)
    covered, i = 0, 0
    while i < len(seq):
        for L in lens:
            if seq[i:i + L] in units:
                covered += L
                i += L
                break
        else:
            i += 1
    return covered / len(seq)


def apply_haplotype(locus_seq: str, locus_start: int, records: List[Tuple[int, str, str]],
                    pad: int = 0, motifs: Optional[List[str]] = None,
                    min_motif_fraction: float = 0.5) -> str:
    """Apply one haplotype's variants to a locus and return the locus sequence.

    locus_seq   : reference sequence of the locus (+ strand), first base at locus_start
    locus_start : 1-based position of locus_seq[0]
    records     : (POS, REF, ALT) this haplotype carries, ALT already chosen

    Replaces `bcftools consensus`, which crashes ("double free or corruption")
    or drops variants on the nested, overlapping records graph callers emit.

    Rules
      - Records are applied in POS order; one that overlaps an already applied
        record is skipped, as bcftools consensus does.
      - Each record is reduced to the REF bases it changes (shared prefix,
        then shared suffix stripped). Only the changed bases inside the locus
        are replaced; a change that straddles an edge keeps its replacement.
      - An insertion on either edge of the locus belongs to the locus. Repeat
        insertions are left-aligned, so they are often anchored on the base just
        before the locus, which fetching from the first locus base missed.
      - An insertion up to `pad` bp outside the locus is moved onto the nearest
        edge when at least `min_motif_fraction` of it is made of the locus
        motif. Graph callers sometimes place a repeat expansion a few bases
        away from the BED boundary and delete the reference repeat instead.
    """
    locus_end = locus_start + len(locus_seq) - 1
    edits: List[Tuple[int, int, str]] = []   # (start idx, end idx exclusive, replacement) in locus coords
    last_end = -1
    for pos, ref, alt in sorted(records, key=lambda r: r[0]):
        if alt in ("*", ".") or alt == ref:
            continue
        if alt.startswith("<") or not re.fullmatch(r"[ACGTNacgtn]+", alt) or not re.fullmatch(r"[ACGTNacgtn]+", ref):
            raise HaplotypeError(f"unsupported allele at {pos}: {ref}>{alt}")
        rec_end = pos + len(ref) - 1
        if pos <= last_end:
            continue                                   # overlaps an applied record
        # Reduce the record to the bases that actually change. Strip the shared
        # prefix first (the VCF anchor base and anything after it that agrees),
        # then the shared suffix of what is left.
        R_, A_ = ref.upper(), alt.upper()
        overlaps_locus = pos <= locus_end and rec_end >= locus_start
        # Inside a repeat the change could sit anywhere; for a record that
        # overlaps the locus, stop stripping before the change leaves it.
        p_max = min(len(R_), len(A_))
        if overlaps_locus:
            p_max = min(p_max, max(0, locus_end - pos))
        p = 0
        while p < p_max and R_[p] == A_[p]:
            p += 1
        q_max = min(len(R_), len(A_)) - p
        if overlaps_locus:
            q_max = min(q_max, max(0, rec_end - locus_start))
        q = 0
        while q < q_max and R_[len(R_) - 1 - q] == A_[len(A_) - 1 - q]:
            q += 1
        core = alt[p:len(alt) - q]
        chg_start = pos + p                            # first changed REF base
        chg_end = rec_end - q                          # last changed REF base (< chg_start: pure insertion)
        last_end = rec_end
        if chg_end < chg_start:
            # pure insertion between chg_start-1 and chg_start. Inside the locus,
            # or on either edge of it, the inserted bases belong to the locus.
            if locus_start <= chg_start <= locus_end + 1:
                i = chg_start - locus_start
                edits.append((i, i, core))
            elif (pad and motifs and locus_start - pad <= chg_start <= locus_end + 1 + pad
                  and motif_fraction(core, motifs) >= min_motif_fraction):
                i = 0 if chg_start < locus_start else len(locus_seq)
                edits.append((i, i, core))
            continue
        if chg_end < locus_start or chg_start > locus_end:
            # outside the locus; a motif-made gain within the padding still counts
            if (pad and motifs and len(core) > chg_end - chg_start + 1
                    and locus_start - pad <= chg_start and chg_end <= locus_end + pad
                    and motif_fraction(core, motifs) >= min_motif_fraction):
                i = 0 if chg_end < locus_start else len(locus_seq)
                edits.append((i, i, core))
            continue
        s_idx = max(chg_start, locus_start) - locus_start
        e_idx = min(chg_end, locus_end) - locus_start + 1
        if core and (chg_start < locus_start or chg_end > locus_end):
            # A replacement that straddles an edge: find where the locus begins
            # and ends inside ALT by matching the reference flank next to it.
            i0 = max(locus_start, pos) - pos             # locus span inside REF
            i1 = min(locus_end, rec_end) - pos + 1
            # Each edge, in order of trust: ALT keeps the outside REF bases as
            # they are; the flank is found once in ALT; else cut by position,
            # which is what bcftools consensus did.
            n_out_right = len(ref) - i1
            a0, a1 = 0, len(alt)
            if pos < locus_start:
                if alt[:i0].upper() == ref[:i0].upper():
                    a0 = i0
                else:
                    a0 = _locate_flank(alt, ref[max(0, i0 - 20):i0], from_right=False)
                    if a0 is None:
                        # by position; an ALT shorter than the outside part means
                        # the deletion runs into the locus, so it gives no bases
                        a0 = min(i0, len(alt))
            if rec_end > locus_end:
                tail = ref[i1:]
                if len(alt) - a0 >= n_out_right and alt[len(alt) - n_out_right:].upper() == tail.upper():
                    a1 = len(alt) - n_out_right
                else:
                    cut = _locate_flank(alt[a0:], ref[i1:i1 + 20], from_right=True)
                    a1 = a0 + cut if cut is not None else max(a0, len(alt) - n_out_right)
            if a1 < a0:
                raise UnresolvedHaplotype(f"record at {pos} spans the locus edge and its ALT cannot be placed")
            edits.append((i0 + pos - locus_start, i1 + pos - locus_start, alt[a0:a1]))
            continue
        edits.append((s_idx, e_idx, core))

    out, cur = [], 0
    for s_idx, e_idx, rep_seq in sorted(edits, key=lambda e: (e[0], e[1])):
        if s_idx < cur:
            raise HaplotypeError("edits overlap after trimming")
        out.append(locus_seq[cur:s_idx])
        out.append(rep_seq)
        cur = e_idx
    out.append(locus_seq[cur:])
    return "".join(out)


# ------------ core region processing ------------

def process_one_region(row_dict: dict, big_vcf_path: str, base_tmp_dir: str, idx: int) -> List[dict]:
    """Process one TR interval and return a list of result dicts (rows)."""
    rows, _status, _key = _process_one_region_impl(row_dict, big_vcf_path, base_tmp_dir, idx)
    return rows


def _process_one_region_impl(
    row_dict: dict, big_vcf_path: str, base_tmp_dir: str, idx: int
) -> Tuple[List[dict], str, tuple]:
    """As process_one_region, but also say WHY a region produced no rows.

    An empty result means either "the VCF has no variant here, so every sample
    is reference" or "this locus was attempted and lost". Only the first is
    safe to treat as reference later on, so the reason is reported to the
    caller and recorded in the run manifest.
    """
    chrom = row_dict.get("CHROM", "")
    original_start = int(row_dict.get("START", 0))
    original_end   = int(row_dict.get("END", 0))
    ref_allele = row_dict.get("TR_allele", "")
    motif = row_dict.get("MOTIF", "")
    key = (chrom, str(original_start), str(original_end), motif)

    is_neg_strand = (original_start > original_end)
    start, end = sorted([original_start, original_end])

    region_start_1based = start + 1
    region_end_1based = end
    region_str = f"{chrom}:{region_start_1based}-{region_end_1based}"

    # Not TemporaryDirectory: on NFS, deleting a file that is still open leaves a
    # .nfsXXXX placeholder, and the cleanup then raises "Directory not empty" and
    # kills the whole run. Ignore cleanup failures instead.
    tmpdir = tempfile.mkdtemp(dir=base_tmp_dir)
    try:
        region_vcf_gz = os.path.join(tmpdir, "region.vcf.gz")
        # Start one base early so an insertion anchored on the base before the
        # locus (how left-aligned repeat expansions are written) is fetched too.
        fetch_str = f"{chrom}:{max(1, region_start_1based - 1 - EDGE_PAD)}-{region_end_1based + EDGE_PAD}"
        if run_command(["bcftools", "view", "--threads", "1", "-r", fetch_str, "-Oz", "-o", region_vcf_gz, big_vcf_path]) is None:
            return [], STATUS_BCFTOOLS_FAILED, key
        if run_command(["tabix", "-p", "vcf", region_vcf_gz]) is None:
            return [], STATUS_BCFTOOLS_FAILED, key

        variant_lines = run_command(["bcftools", "view", "--threads", "1", "-H", region_vcf_gz])
        if not variant_lines or not variant_lines.strip():
            log.info("Region %s: no variants, skipping.", region_str)
            return [], STATUS_NO_VARIANTS, key

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

        sample_signatures = {s: [] for s in expected_samples}
        missing_samples = set()

        rec_lines = run_command(["bcftools", "query", "-f", "%POS\\t%REF\\t%ALT[\\t%GT]\\n", splitted_vcf]) or ""
        region_records: List[Tuple[int, str, List[str]]] = []   # (POS, REF, ALTs) in query order
        if rec_lines:
            for line in rec_lines.strip().splitlines():
                cols = line.split("\t")
                r_pos, r_ref, r_alts = int(cols[0]), cols[1], cols[2].split(",")
                # Records in the padding only matter when they insert bases;
                # keeping the rest would just split samples into more groups.
                touches = r_pos + len(r_ref) - 1 >= region_start_1based - 1 and r_pos <= region_end_1based
                if not touches and not any(len(a) > len(r_ref) and motif_fraction(a, motif.split(",")) >= 0.5
                                           for a in r_alts):
                    continue
                region_records.append((r_pos, r_ref, r_alts))
                cols = cols[3:]
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

        if not region_records:
            # only padding records that insert nothing: same as no variant here
            log.info("Region %s: no variants, skipping.", region_str)
            return [], STATUS_NO_VARIANTS, key

        normalized_signatures: Dict[str, str] = {}
        for s in expected_samples:
            if not sample_signatures[s]:
                sample_signatures[s] = ["0"]
            normalized_signatures[s] = "|".join(sample_signatures[s])

        ref_variant_count = max(1, len(next(iter(sample_signatures.values()))))
        ref_signature = "|".join(["0"] * ref_variant_count)
        if REF_GROUP is None:
            log.error("REF_GROUP is not set. First group file must be the reference group.")
            return [], STATUS_NO_REF_GROUP, key
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
                    carried = []
                    for (pos, ref, alts), gt in zip(region_records, sig.split("|")):
                        if gt in ("0", "."):
                            continue
                        carried.append((pos, ref, alts[int(gt) - 1]))
                    try:
                        plus_strand_seq = apply_haplotype(bcftools_ref_allele, region_start_1based, carried,
                                                          pad=EDGE_PAD, motifs=motif.split(","))
                    except UnresolvedHaplotype as e:
                        # Same treatment as a missing GT: assume REF and flag the
                        # samples as +name+ so the call is visibly an assumption.
                        log.warning("Region %s: group %s assumed REF: %s", region_str, grp_idx, e)
                        plus_strand_seq = bcftools_ref_allele
                        missing_samples.update(group_samps)
                    final_seq = reverse_complement(plus_strand_seq) if is_neg_strand else plus_strand_seq
                except (HaplotypeError, ValueError, IndexError) as e:
                    log.error("Region %s: group %s processing failed: %s", region_str, grp_idx, e)
                    return [], STATUS_CONSENSUS_FAILED, key

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
            return [], STATUS_COUNT_MISMATCH, key

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

        return df_sub.to_dict(orient="records"), STATUS_OK, key
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

# ------------ public API (called by CLI) ------------

def vcf_alleles(
    vcf: str,
    input_tsv: str,
    output_tsv: str,
    group_files: List[str],
    processes: int = 32,
    chunk_size: int = 10_000,
    tmpdir: str = "./tmp",
    write_run_manifest: bool = True
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
    check_vcf_readable(vcf)
    Path(tmpdir).mkdir(parents=True, exist_ok=True)
    atexit.register(lambda: shutil.rmtree(tmpdir, ignore_errors=True))

    load_sample_groups(group_files)

    df = pd.read_csv(input_tsv, sep="\t")
    df["START"] = df["START"].astype(int)
    df["END"]   = df["END"].astype(int)

    all_rows: List[dict] = []
    status_counts: Dict[str, int] = defaultdict(int)
    absent_loci: List[tuple] = []

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
                futs.append(ex.submit(_process_one_region_impl, row.to_dict(), vcf, tmpdir, int(idx_row)))
            for fut in as_completed(futs):
                sub, status, key = fut.result()
                all_rows.extend(sub or [])
                status_counts[status] += 1
                if status != STATUS_OK:
                    absent_loci.append((*key, status))

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

    log.info(
        "Loci: %d emitted, %d with no variants, %d lost to errors (of %d)",
        status_counts.get(STATUS_OK, 0), status_counts.get(STATUS_NO_VARIANTS, 0),
        sum(v for k, v in status_counts.items() if k not in (STATUS_OK, STATUS_NO_VARIANTS)),
        n,
    )

    if write_run_manifest:
        # Records what merge-alleles would otherwise have to guess: the roster,
        # the reference group, and which loci are absent because the VCF was
        # quiet there rather than because the locus failed.
        absent_loci.sort(key=lambda r: int(r[1]) if str(r[1]).lstrip("-").isdigit() else 0)
        status_file = write_loci_status(output_tsv, absent_loci)
        write_manifest(output_tsv, {
            "stage": "vcf-alleles",
            "level": "alleles",
            "chrom": str(df["CHROM"].iloc[0]) if not df.empty else "",
            "output": {"path": os.path.abspath(str(outp)), "n_rows": int(len(df_out))},
            "input_bed": file_stat(input_tsv),
            "vcf": file_stat(vcf),
            "group_order": list(GROUP_ORDER),
            "ref_group": REF_GROUP,
            "ref_pseudo_sample": REF_GROUP,
            "groups": build_group_payload(group_files, SAMPLE_GROUPS, GROUP_ORDER),
            "expected_tokens_per_locus": sum(p for _g, p in SAMPLE_GROUPS.values()),
            "loci_status": dict(sorted(status_counts.items())),
            "loci_status_file": os.path.basename(status_file) if status_file else None,
        })

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
        tmpdir=args.tmpdir,
        write_run_manifest=not getattr(args, "no_manifest", False)
    )
