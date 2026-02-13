# src/trias/bed/prepare.py
from __future__ import annotations
import os, sys, gzip, logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import subprocess
import pysam
from typing import Dict, List, Tuple, Iterable, Optional

log = logging.getLogger("trias")
if not log.handlers:
    h = logging.StreamHandler(sys.stderr)
    h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    log.addHandler(h)
log.setLevel(logging.INFO)

# -------- utilities --------

def open_maybe_gzip(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r", encoding="utf-8")

def ensure_dir(p: str | Path) -> Path:
    p = Path(p); p.mkdir(parents=True, exist_ok=True); return p

def check_dependencies():
    """Verify external dependencies (samtools) are available."""
    try:
        subprocess.run(["samtools", "--version"], capture_output=True, check=True)
    except Exception:
        log.error("samtools is required and must be available on PATH.")
        sys.exit(1)

def index_fasta(fasta: str | Path):
    """Create FASTA index if missing."""
    fai = str(fasta) + ".fai"
    if not os.path.exists(fai):
        log.info(f"Creating FASTA index: {fasta}")
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)

def load_mapping(mapping_path: Optional[str]) -> Dict[str, str]:
    """
    Load chromosome name mapping (BED_name<TAB>FASTA_name).
    Lines starting with '#' are ignored.
    """
    m: Dict[str, str] = {}
    if not mapping_path:
        return m
    if not os.path.exists(mapping_path):
        log.error(f"Mapping file not found: {mapping_path}")
        sys.exit(1)
    with open(mapping_path, "r", encoding="utf-8") as f:
        for ln, line in enumerate(f, 1):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) >= 2:
                m[parts[0]] = parts[1]
            else:
                log.warning(f"Malformed mapping (line {ln}): {s}")
    if m:
        log.info(f"Loaded {len(m)} chromosome mappings.")
    return m

# -------- core FASTA extraction --------

def extract_chunk(
    entries: List[Tuple[int, str, int, int, str]],
    fasta_path: str,
    chrom_map: Dict[str, str]
) -> List[Tuple[int, str, int, int, str, str]]:
    """
    entries: (idx, chrom, start, end, motif)
    returns: (idx, chrom, start, end, motif, sequence)
    """
    out: List[Tuple[int, str, int, int, str, str]] = []
    fa = pysam.FastaFile(fasta_path)
    for idx, chrom, start, end, motif in entries:
        try:
            fasta_chrom = chrom_map.get(chrom, chrom)
            seq = fa.fetch(fasta_chrom, start, end).upper()  # BED is 0-based, end-exclusive
        except Exception as e:
            log.warning(f"Sequence fetch failed {chrom}:{start}-{end}: {e}")
            seq = "N/A"
        out.append((idx, chrom, start, end, motif, seq))
    fa.close()
    return out

def write_tr_bed(path: str | Path, rows: Iterable[Tuple[str, int, int, str, str]]):
    """Write a TR BED-like table with header."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("CHROM\tSTART\tEND\tMOTIF\tTR_allele\n")
        for chrom, start, end, motif, seq in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{motif}\t{seq}\n")

# -------- BED parsing & per-chromosome processing --------

def parse_bed_rows(bed_path: str, motif_col_1based: int, has_header: bool) -> List[Tuple[str, int, int, str]]:
    """
    Parse BED (optionally gzipped) and return rows as (chrom, start, end, motif).
    - motif_col_1based: 1-based column index for motif; empty string if missing.
    - has_header: skip the first line if true.
    """
    rows: List[Tuple[str, int, int, str]] = []
    col_idx = motif_col_1based - 1
    with open_maybe_gzip(bed_path) as f:
        for ln, line in enumerate(f, 1):
            if ln == 1 and has_header:
                continue
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            parts = s.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end   = int(parts[2])
            except ValueError:
                # Likely a header or malformed row
                continue
            motif = parts[col_idx] if 0 <= col_idx < len(parts) else ""
            rows.append((chrom, start, end, motif))
    return rows

def unique_chroms(rows: List[Tuple[str, int, int, str]]) -> List[str]:
    """Collect unique chromosome names and return a sorted list."""
    seen = {}
    for chrom, _, _, _ in rows:
        seen[chrom] = 1
    return sorted(seen.keys(), key=lambda c: (c.lstrip("chr"), c))

def process_one_chrom(
    chrom: str,
    rows: List[Tuple[str, int, int, str]],
    fasta: str,
    chrom_map: Dict[str, str],
    out_path: str,
    threads: int,
    chunk: int
):
    # Filter rows for the target chromosome
    sub = [(i, c, s, e, m) for i, (c, s, e, m) in enumerate(rows) if c == chrom]
    if not sub:
        return
    # Split into chunks
    chunks = [sub[i:i+chunk] for i in range(0, len(sub), chunk)]
    results: List[Tuple[int, str, int, int, str, str]] = []
    with ThreadPoolExecutor(max_workers=threads) as ex:
        futs = {ex.submit(extract_chunk, ck, fasta, chrom_map): k for k, ck in enumerate(chunks)}
        for fut in as_completed(futs):
            results.extend(fut.result())
    # Restore original order
    results.sort(key=lambda x: x[0])
    # Write file
    write_tr_bed(out_path, ((chrom, s, e, m, seq) for _, _, s, e, m, seq in results))
    log.info(f"Saved: {out_path}  (n={len(results)})")

def bed_prepare(
    bed: str,
    fasta: str,
    outdir: str,
    threads: int = 4,
    chunk: int = 1000,
    mapping: Optional[str] = None,
    motif_col_1based: int = 10,
    has_header: bool = True
):
    """
    Generate per-chromosome TR BED files from BED(+gz) and a FASTA.
    Output files: <outdir>/<CHROM>_TR.bed with header:
      CHROM  START  END  MOTIF  TR_allele
    """
    check_dependencies()
    index_fasta(fasta)
    chrom_map = load_mapping(mapping)
    ensure_dir(outdir)

    log.info(f"BED: {bed}")
    log.info(f"FASTA: {fasta}")
    log.info(f"OUTDIR: {outdir}")
    if mapping:
        log.info(f"MAPPING: {mapping}")
    log.info(f"threads={threads} chunk={chunk} motif_col={motif_col_1based} header={has_header}")

    rows = parse_bed_rows(bed, motif_col_1based, has_header)
    if not rows:
        log.error("No valid rows found in the BED. Aborting.")
        sys.exit(1)

    chroms = unique_chroms(rows)
    preview = ", ".join(chroms[:6]) + (" ..." if len(chroms) > 6 else "")
    log.info(f"Chromosomes to process: {len(chroms)} ({preview})")

    for chrom in chroms:
        out_path = str(Path(outdir) / f"{chrom}_TR.bed")
        process_one_chrom(chrom, rows, fasta, chrom_map, out_path, threads, chunk)

# ----- CLI wrapper -----

def bed_prepare_cmd(args):
    bed_prepare(
        bed=args.bed,
        fasta=args.fasta,
        outdir=args.outdir,
        threads=args.threads,
        chunk=args.chunk,
        mapping=args.mapping,
        motif_col_1based=args.motif_col,
        has_header=args.has_header
    )
