from __future__ import annotations
import math
import logging
from pathlib import Path
from typing import List, Tuple, Dict, Optional

import pandas as pd
import numpy as np
import multiprocessing as mp

log = logging.getLogger("trias")
if not log.handlers:
    import sys as _sys
    h = logging.StreamHandler(_sys.stderr)
    h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    log.addHandler(h)
log.setLevel(logging.INFO)

# Core columns (excluded from auto-detected count columns)
CORE_COLS = {
    "CHROM","START","END","MOTIF","TR_allele",
    "total_repeat_dna","max_consec_dna","Samples"
}

# ---------- file/group helpers ----------
def read_group_list(path: Optional[str]) -> List[str]:
    """Read group names (one per line). Lines starting with '#' or empty lines are ignored."""
    if not path:
        return []
    p = Path(path).expanduser()
    if not p.exists():
        log.warning("Group list file not found: %s (ignored)", path)
        return []
    names = []
    with p.open() as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            names.append(s)
    return names

def label_from_path(path: str) -> str:
    """Derive group label from file basename without extension."""
    return Path(path).stem

def detect_group_count_columns(df: pd.DataFrame) -> List[str]:
    """Auto-detect numeric-like count columns excluding CORE_COLS."""
    cols = []
    for c in df.columns:
        if c in CORE_COLS:
            continue
        ser = pd.to_numeric(df[c], errors="coerce")
        if ser.notna().any():
            cols.append(c)
    return cols

# ---------- reverse complement & counts parsing ----------
def reverse_complement(seq: str) -> str:
    if pd.isna(seq) or seq == "." or seq == "":
        return seq
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]

def parse_motif_counts(count_str):
    if pd.isna(count_str) or count_str == "." or count_str == "":
        return []
    try:
        s = str(count_str)
        if "_" in s:
            return [int(x) for x in s.split("_")]
        else:
            return [int(s)]
    except:
        return []

def format_min_max_range(values_list):
    if not values_list:
        return "."
    if len(values_list) == 1:
        return str(values_list[0])
    mn, mx = min(values_list), max(values_list)
    return str(mn) if mn == mx else f"{mn}-{mx}"

def reverse_counts_string(count_str: str, motif_count: int) -> str:
    vals = parse_motif_counts(count_str)
    if not vals:
        return "."
    if motif_count <= 1:
        return "_".join(str(x) for x in vals)
    if len(vals) == motif_count:
        vals = vals[::-1]
    else:
        chunks = [vals[i:i+motif_count] for i in range(0, len(vals), motif_count)]
        vals = [x for ch in (ch[::-1] for ch in chunks) for x in ch]
    return "_".join(str(x) for x in vals)

def normalize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """Fix START>END and sync MOTIF + per-motif count strings accordingly."""
    df = df.copy()
    df["START"] = pd.to_numeric(df["START"], errors="coerce")
    df["END"]   = pd.to_numeric(df["END"], errors="coerce")

    mask = df["START"] > df["END"]
    if not mask.any():
        return df

    flipped = df.loc[mask].copy()

    # swap coords
    flipped[["START","END"]] = flipped[["END","START"]].values

    # reverse-complement motif and order
    def rc_motif_list(motif_str: str) -> str:
        if pd.isna(motif_str) or motif_str == "." or motif_str == "":
            return motif_str
        toks = str(motif_str).split(",")
        rc = [reverse_complement(t) for t in toks][::-1]
        return ",".join(rc)

    old = flipped["MOTIF"].astype(str).fillna(".")
    motif_counts = old.apply(lambda s: 0 if s in [".",""] else len(s.split(",")))
    flipped["MOTIF"] = old.apply(rc_motif_list)

    # reverse per-motif count strings
    for col in ("total_repeat_dna","max_consec_dna"):
        if col in flipped.columns:
            flipped[col] = [
                reverse_counts_string(val, mc)
                for val, mc in zip(flipped[col].astype(str).fillna("."), motif_counts)
            ]

    df.loc[mask] = flipped
    return df

# ---------- math helpers ----------
def format_multi_motif_range(values: List[int], motif_count: int) -> str:
    if not values:
        return "."
    if motif_count == 1:
        return format_min_max_range(values)
    parts = []
    for i in range(motif_count):
        sub = [values[j] for j in range(i, len(values), motif_count)]
        parts.append(format_min_max_range(sub))
    return "_".join(parts)

def weighted_mean(values_with_weights: List[tuple[int,int]], motif_count: int) -> str:
    if not values_with_weights:
        return "."
    if motif_count == 1:
        num = sum(v*w for v,w in values_with_weights)
        den = sum(w for v,w in values_with_weights)
        return "." if den == 0 else f"{num/den:.2f}"
    parts = []
    for i in range(motif_count):
        sub = [(v,w) for idx,(v,w) in enumerate(values_with_weights) if idx % motif_count == i]
        if not sub:
            parts.append(".")
            continue
        num = sum(v*w for v,w in sub); den = sum(w for v,w in sub)
        parts.append("." if den == 0 else f"{num/den:.2f}")
    return "_".join(parts)

# ---------- per-TR group summarizer ----------
def summarize_one_group(group_df: pd.DataFrame,
                        pop_cols: List[str],
                        custom_groups: Dict[str, List[str]]) -> Dict[str, object]:
    first = group_df.iloc[0]
    chrom = first["CHROM"]; start = int(first["START"]); end = int(first["END"]); motif = first["MOTIF"]

    # prepare numeric safe frame for weights
    safe = group_df.copy()
    for c in set(pop_cols + [c for cols in custom_groups.values() for c in cols]):
        if c in safe.columns:
            safe[c] = pd.to_numeric(safe[c], errors="coerce").fillna(0)

    # motif multiplicity
    motif_count = 1
    if isinstance(motif, str) and motif not in (".","") and "," in motif:
        motif_count = len(motif.split(","))

    def collect_for_columns(cols: List[str], label: str) -> Dict[str, object]:
        if not cols:
            return {
                f"{label}_allele_count": 0,
                f"{label}_total_repeat_dna_range": ".",
                f"{label}_max_consec_dna_range": ".",
                f"{label}_total_repeat_dna_mean": ".",
                f"{label}_max_consec_dna_mean": ".",
            }
        # mask rows that have any alleles in these cols
        mask = safe[cols].sum(axis=1) > 0
        rows = group_df.loc[mask]
        if rows.empty:
            return {
                f"{label}_allele_count": 0,
                f"{label}_total_repeat_dna_range": ".",
                f"{label}_max_consec_dna_range": ".",
                f"{label}_total_repeat_dna_mean": ".",
                f"{label}_max_consec_dna_mean": ".",
            }

        weights = safe.loc[mask, cols].sum(axis=1).astype(int)

        totals_flat: List[int] = []
        maxs_flat: List[int]   = []
        totals_w: List[tuple[int,int]] = []
        maxs_w: List[tuple[int,int]]   = []

        for (idx, row), w in zip(rows.iterrows(), weights):
            if w <= 0:
                continue
            tlist = parse_motif_counts(row.get("total_repeat_dna", "."))
            mlist = parse_motif_counts(row.get("max_consec_dna", "."))
            if tlist:
                totals_flat.extend(tlist)
                totals_w.extend((v, w) for v in tlist)
            if mlist:
                maxs_flat.extend(mlist)
                maxs_w.extend((v, w) for v in mlist)

        allele_count = int(safe.loc[mask, cols].sum(axis=1).sum())

        return {
            f"{label}_allele_count": allele_count,
            f"{label}_total_repeat_dna_range": format_multi_motif_range(totals_flat, motif_count),
            f"{label}_max_consec_dna_range":   format_multi_motif_range(maxs_flat,   motif_count),
            f"{label}_total_repeat_dna_mean":  weighted_mean(totals_w, motif_count),
            f"{label}_max_consec_dna_mean":    weighted_mean(maxs_w,   motif_count),
        }

    # Population = all numeric count columns
    res = {
        "CHROM": chrom, "START": start, "END": end, "MOTIF": motif,
    }
    res.update(collect_for_columns(pop_cols, "Population"))
    # Custom groups
    for label, cols in custom_groups.items():
        res.update(collect_for_columns(cols, label))
    return res

# ---------- top-level worker for multiprocessing ----------
def _process_chunk_worker(args):
    """
    Top-level worker for multiprocessing (must be picklable).
    args = (chunk, pop_cols, custom_groups)
    """
    chunk, pop_cols, custom_groups = args
    out = []
    for _, gdf in chunk:
        try:
            out.append(summarize_one_group(gdf, pop_cols, custom_groups))
        except Exception as e:
            log.error("Group failed: %s", e)
    return out

# ---------- main API ----------
def motif_summary(
    input_path: str,
    output_path: str,
    group_files: Optional[List[str]] = None,
    n_processes: Optional[int] = None,
    chunk_size: Optional[int] = None,
):
    """
    Summarize per-TR motif counts.

    Population: sum over ALL auto-detected numeric count columns (excluding core columns).
    Custom groups: provided via --groups <file>, label = basename(file), columns = names listed in the file.
    Required input columns: CHROM, START, END, MOTIF, total_repeat_dna, max_consec_dna (+ arbitrary count columns)
    """
    in_path = Path(input_path).expanduser()
    out_path = Path(output_path).expanduser()
    if str(out_path.parent) not in ("","."):
        out_path.parent.mkdir(parents=True, exist_ok=True)

    log.info("Reading: %s", in_path)
    df = pd.read_csv(in_path, sep="\t")
    log.info("Loaded %d rows", len(df))

    # normalize
    log.info("Normalizing coordinates/MOTIF")
    df = normalize_dataframe(df)

    # auto-detect population columns
    pop_cols = detect_group_count_columns(df)
    log.info("Population columns detected: %s", ", ".join(pop_cols) if pop_cols else "(none)")

    # build custom groups
    custom_groups: Dict[str, List[str]] = {}
    group_files = group_files or []
    for gf in group_files:
        label = label_from_path(gf)
        cols = [c for c in read_group_list(gf) if c in df.columns]
        if not cols:
            log.warning("Group '%s' has no matching columns in input; it will be zeros.", label)
        custom_groups[label] = cols
        log.info("Group '%s' columns: %s", label, ", ".join(cols) if cols else "(none)")

    # group by TR key
    key_cols = ["CHROM","START","END","MOTIF"]
    groups = list(df.groupby(key_cols, sort=False))
    total = len(groups)
    log.info("Processing %d TR groups", total)

    if n_processes is None:
        n_processes = mp.cpu_count()
    if chunk_size is None:
        chunk_size = max(1, total // (n_processes * 4) if total else 1)

    # chunking
    chunks: List[List[Tuple[tuple, pd.DataFrame]]] = []
    for i in range(0, total, chunk_size):
        chunks.append(groups[i:i+chunk_size])
    log.info("Created %d chunks (avg ~%d groups/chunk)", len(chunks), (total//len(chunks) if chunks else 0))

    # run
    results: List[Dict[str, object]] = []
    job_args = [(ch, pop_cols, custom_groups) for ch in chunks]

    if n_processes == 1 or len(chunks) <= 1:
        log.info("Using single-process mode")
        for args in job_args:
            results.extend(_process_chunk_worker(args))
    else:
        try:
            with mp.Pool(processes=n_processes) as pool:
                for part in pool.map(_process_chunk_worker, job_args):
                    results.extend(part)
        except Exception as e:
            log.error("Multiprocessing error: %s; falling back to single-process.", e)
            results.clear()
            for args in job_args:
                results.extend(_process_chunk_worker(args))

    # assemble output
    out_df = pd.DataFrame(results)

    # Column order: Population first, then each custom group, then metadata/ranges/means
    columns = [
        "Population_allele_count", "CHROM","START","END","MOTIF",
        "Population_total_repeat_dna_range","Population_max_consec_dna_range",
        "Population_total_repeat_dna_mean","Population_max_consec_dna_mean",
    ]
    # add each custom group’s fields
    for label in custom_groups.keys():
        columns.insert(1, f"{label}_allele_count")  # keep counts up front
        columns += [
            f"{label}_total_repeat_dna_range", f"{label}_max_consec_dna_range",
            f"{label}_total_repeat_dna_mean",  f"{label}_max_consec_dna_mean",
        ]
    # keep only present
    columns = [c for c in columns if c in out_df.columns]
    out_df = out_df[columns]

    log.info("Writing: %s", out_path)
    out_df.to_csv(str(out_path), sep="\t", index=False)
    log.info("Done.")
    return str(out_path)

# ---------- CLI adapter ----------
def motif_summary_cmd(args):
    return motif_summary(
        input_path=args.input,
        output_path=args.output,
        group_files=args.groups or [],
        n_processes=args.processes,
        chunk_size=args.chunk_size,
    )
