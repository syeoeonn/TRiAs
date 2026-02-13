# src/trias/motif/count.py
from __future__ import annotations
import re
import math
import logging
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, List, Tuple, Iterable, Optional
from pathlib import Path

import pandas as pd

log = logging.getLogger("trias")
if not log.handlers:
    import sys as _sys
    h = logging.StreamHandler(_sys.stderr)
    h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    log.addHandler(h)
log.setLevel(logging.INFO)

# ----------------- regex cache -----------------
_regex_cache: Dict[str, re.Pattern] = {}

def get_regex(motif: str) -> re.Pattern:
    """Get or compile and cache a regex pattern for a motif."""
    pat = _regex_cache.get(motif)
    if pat is None:
        pat = re.compile(f"(?:{re.escape(motif)})+")
        _regex_cache[motif] = pat
    return pat

# ----------------- single-motif compress -----------------
def compress_text_single_motif(text: str, motif: str):
    """
    Compress a sequence using a single motif and return:
      (compressed_text, total_repeats, max_consec)
    """
    if not motif:
        return text, 0, 0
    regex = get_regex(motif)
    total_repeats = 0
    max_consec = 0

    def repl(m: re.Match):
        nonlocal total_repeats, max_consec
        block = m.group(0)
        count = len(block) // len(motif)
        total_repeats += count
        max_consec = max(max_consec, count)
        clean = motif.strip("()")
        return f"({clean})" if count == 1 else f"({clean}){count}"

    compressed_text = regex.sub(repl, text)
    return compressed_text, total_repeats, max_consec

# ----------------- interval utils -----------------
def _merge_intervals(ivals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Merge overlapping intervals (end exclusive)."""
    if not ivals: return []
    ivals = sorted(ivals)
    out = [ivals[0]]
    for s, e in ivals[1:]:
        if s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out

def _invert_intervals(n: int, ivals: List[Tuple[int,int]]) -> List[Tuple[int,int]]:
    """Return gaps (complement) within [0, n) given non-overlapping intervals."""
    gaps: List[Tuple[int,int]] = []
    cur = 0
    for s, e in ivals:
        if cur < s:
            gaps.append((cur, s))
        cur = max(cur, e)
    if cur < n:
        gaps.append((cur, n))
    return gaps

# ----------------- WIS (weighted interval scheduling) -----------------
def _select_non_overlapping_best(matches: List[Tuple[int,int,str,int]],
                                 weight_fn=lambda m: (m[1]-m[0], len(m[2]), m[3], 1)):
    """
    matches: list of (start, end, motif, count), end exclusive.
    Returns chosen list with max lexicographic weight:
      (covered_bp, motif_len, count, blocks)
    """
    if not matches: return []
    matches_sorted = sorted(matches, key=lambda m: m[1])
    ends = [m[1] for m in matches_sorted]

    import bisect
    p = []
    for i, (s, e, _, _) in enumerate(matches_sorted):
        j = bisect.bisect_right(ends, s) - 1
        p.append(j)

    def t_add(a, b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3])
    def t_gt(a, b):
        if a[0] != b[0]: return a[0] > b[0]
        if a[1] != b[1]: return a[1] > b[1]
        if a[2] != b[2]: return a[2] > b[2]
        return a[3] > b[3]

    n = len(matches_sorted)
    W = [weight_fn(m) for m in matches_sorted]
    dp = [(0,0,0,0)] * (n+1)
    choose = [False] * n

    for i in range(1, n+1):
        opt1 = dp[i-1]
        prev = p[i-1]
        opt2 = t_add(W[i-1], dp[prev+1] if prev >= 0 else (0,0,0,0))
        if t_gt(opt2, opt1):
            dp[i] = opt2; choose[i-1] = True
        else:
            dp[i] = opt1

    chosen = []
    i = n
    while i > 0:
        if choose[i-1]:
            chosen.append(matches_sorted[i-1])
            i = p[i-1] + 1
        else:
            i -= 1
    chosen.reverse()
    return chosen

# ----------------- gap greedy -----------------
def _compress_gap_greedy(text_gap: str, motifs: List[str], regex_cache: Dict[str,re.Pattern]):
    """
    Compress a gap greedily with all motifs (longer first), non-overlapping maximal runs.
    Returns: (compressed_gap_text, delta_stats_per_motif)
    """
    if not text_gap:
        return "", {}
    motifs_sorted = sorted(motifs, key=len, reverse=True)
    matches: List[Tuple[int,int,str,int]] = []
    for motif in motifs_sorted:
        pat = regex_cache.get(motif)
        if pat is None:
            pat = re.compile(f"(?:{re.escape(motif)})+")
            regex_cache[motif] = pat
        for m in pat.finditer(text_gap):
            s, e = m.span()
            cnt = (e - s) // len(motif)
            if cnt >= 1:
                matches.append((s, e, motif, cnt))

    matches.sort(key=lambda x: (x[0], -len(x[2]), -x[3]))  # start asc, longer motif first, larger count first
    chosen = []
    last_end = -1
    for s, e, motif, cnt in matches:
        if s >= last_end:
            chosen.append((s, e, motif, cnt))
            last_end = e

    chosen.sort(key=lambda x: x[0])
    res = []
    cur = 0
    delta = {m: {'total': 0, 'max_consec': 0} for m in motifs}
    for s, e, motif, cnt in chosen:
        if s > cur: res.append(text_gap[cur:s])
        clean = motif.strip("()")
        res.append(f"({clean})" if cnt == 1 else f"({clean}){cnt}")
        delta[motif]['total'] += cnt
        delta[motif]['max_consec'] = max(delta[motif]['max_consec'], cnt)
        cur = e
    if cur < len(text_gap): res.append(text_gap[cur:])
    return "".join(res), delta

# ----------------- multi-motif compressor -----------------
def compress_with_multi_motif_separated(text: str, motif_str: str):
    """
    Multi-motif compression + per-motif separated statistics.
    Pass-1: global optimal non-overlapping tiling via WIS (lexicographic weight).
    Pass-2: greedily compress gaps with all motifs (longer first).
    Returns:
      (compressed_text, total_repeats, max_consec, motif_stats)
      where motif_stats: {motif: {'total': int, 'max_consec': int}}
    """
    if "," not in motif_str:
        comp, total, mx = compress_text_single_motif(text, motif_str.strip())
        return comp, total, mx, {motif_str.strip(): {'total': total, 'max_consec': mx}}

    motifs = [m.strip() for m in motif_str.split(",") if m.strip()]
    if not motifs:
        return text, 0, 0, {}

    # pre-compile
    patterns: List[Tuple[str, re.Pattern]] = []
    for motif in motifs:
        if motif not in _regex_cache:
            _regex_cache[motif] = re.compile(f"(?:{re.escape(motif)})+")
        patterns.append((motif, _regex_cache[motif]))

    # collect matches
    matches: List[Tuple[int,int,str,int]] = []
    for motif, pat in patterns:
        for m in pat.finditer(text):
            s, e = m.span()
            cnt = (e - s) // len(motif)
            if cnt >= 1:
                matches.append((s, e, motif, cnt))

    def _weight(m):
        s,e,motif,cnt = m
        return (e-s, len(motif), cnt, 1)  # covered bp, motif length, repeat count, blocks

    chosen = _select_non_overlapping_best(matches, weight_fn=_weight)
    chosen.sort(key=lambda x: x[0])

    locked = _merge_intervals([(s,e) for s,e,_,_ in chosen])
    gaps = _invert_intervals(len(text), locked)

    gap_snippets: Dict[Tuple[int,int], str] = {}
    gap_stats = {m: {'total': 0, 'max_consec': 0} for m in motifs}
    for g0, g1 in gaps:
        comp_gap, delta = _compress_gap_greedy(text[g0:g1], motifs, _regex_cache)
        gap_snippets[(g0,g1)] = comp_gap
        for m in motifs:
            gap_stats[m]['total'] += delta.get(m, {}).get('total', 0)
            gap_stats[m]['max_consec'] = max(gap_stats[m]['max_consec'], delta.get(m, {}).get('max_consec', 0))

    result = []
    last_end = 0
    total_repeats = 0
    max_consec = 0
    motif_stats = {m: {'total': 0, 'max_consec': 0} for m in motifs}

    def _flush_gap(upto: int):
        nonlocal last_end
        if last_end < upto:
            result.append(gap_snippets.get((last_end,upto), text[last_end:upto]))
            last_end = upto

    for s, e, motif, cnt in chosen:
        _flush_gap(s)
        clean = motif.strip("()")
        result.append(f"({clean})" if cnt == 1 else f"({clean}){cnt}")
        motif_stats[motif]['total'] += cnt
        motif_stats[motif]['max_consec'] = max(motif_stats[motif]['max_consec'], cnt)
        total_repeats += cnt
        max_consec = max(max_consec, cnt)
        last_end = e

    _flush_gap(len(text))

    # merge gap stats
    for m in motifs:
        motif_stats[m]['total'] += gap_stats[m]['total']
        motif_stats[m]['max_consec'] = max(motif_stats[m]['max_consec'], gap_stats[m]['max_consec'])
        total_repeats += gap_stats[m]['total']
        max_consec = max(max_consec, gap_stats[m]['max_consec'])

    return "".join(result), total_repeats, max_consec, motif_stats

def format_motif_separated_stats(motif_str: str, motif_stats: Dict[str, Dict[str,int]]):
    """
    Convert per-motif stats to '_' separated strings following motif order in motif_str.
    Returns: (total_str, max_consec_str)
    """
    if "," not in motif_str:
        m = motif_str.strip()
        s = motif_stats.get(m, {'total': 0, 'max_consec': 0})
        return str(s['total']), str(s['max_consec'])

    motifs = [m.strip() for m in motif_str.split(",") if m.strip()]
    total_parts, max_parts = [], []
    for m in motifs:
        s = motif_stats.get(m, {'total': 0, 'max_consec': 0})
        total_parts.append(str(s['total']))
        max_parts.append(str(s['max_consec']))
    return "_".join(total_parts), "_".join(max_parts)

# ----------------- per-group worker -----------------
def process_group(args):
    """Process a single DNA-only group of tandem repeat data."""
    group_key, group_df = args
    group_df = group_df.copy()

    rep_motif_dna = str(group_df["MOTIF"].iloc[0])

    for idx, row in group_df.iterrows():
        if row["TR_allele"] != "." and rep_motif_dna not in (".", "nan"):
            comp_dna, total_dna, maxconsec_dna, stats_dna = compress_with_multi_motif_separated(
                str(row["TR_allele"]), rep_motif_dna
            )
            total_dna_sep, max_dna_sep = format_motif_separated_stats(rep_motif_dna, stats_dna)
        else:
            comp_dna = row["TR_allele"]
            total_dna_sep = max_dna_sep = "."

        group_df.at[idx, "TR_allele"] = comp_dna
        group_df.at[idx, "total_repeat_dna"] = total_dna_sep
        group_df.at[idx, "max_consec_dna"] = max_dna_sep

    return group_df

def chunkify(lst, chunk_size):
    """Yield chunks of a list."""
    for i in range(0, len(lst), chunk_size):
        yield lst[i:i + chunk_size]

# ----------------- public API -----------------
def motif_count(
    input_path: str,
    output_path: str,
    max_workers: int = 100,
    group_chunk_size: int = 8000
):
    """
    Count per-motif repeats from DNA-level TR alleles using multi-motif compression.

    Required columns:
      CHROM, START, END, MOTIF, TR_allele
    Grouping is performed over: CHROM, START, END, MOTIF

    All other columns (e.g., user-defined population/group columns) are preserved as-is.
    """
    df = pd.read_csv(input_path, sep="\t")

    # Preserve original column order (keeps user-defined group columns)
    original_cols = list(df.columns)

    # Ensure result columns exist (if absent)
    for col in ["total_repeat_dna", "max_consec_dna"]:
        if col not in df.columns:
            df[col] = None

    # DNA-only grouping columns
    group_cols = ["CHROM", "START", "END", "MOTIF"]
    groups = list(df.groupby(group_cols, sort=False))
    total_groups = len(groups)
    total_chunks = math.ceil(total_groups / group_chunk_size)
    log.info("Total number of groups: %d, total chunks: %d", total_groups, total_chunks)

    first_chunk = True
    outp = Path(output_path).expanduser()
    if str(outp.parent) not in ("", "."):
        outp.parent.mkdir(parents=True, exist_ok=True)

    for idx, group_chunk in enumerate(chunkify(groups, group_chunk_size), start=1):
        log.info("[%s] Processing chunk %d/%d (groups in this chunk: %d)",
                 datetime.now().strftime("%Y-%m-%d %H:%M:%S"), idx, total_chunks, len(group_chunk))
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(process_group, group_chunk))
        processed_chunk = pd.concat(results, ignore_index=True)

        # Build output column order dynamically:
        # original columns (unchanged) + our new result columns placed after TR_allele if possible.
        cols = list(original_cols)
        for extra in ["total_repeat_dna", "max_consec_dna"]:
            if extra not in cols:
                if "TR_allele" in cols:
                    insert_at = cols.index("TR_allele") + 1
                    cols.insert(insert_at, extra)
                else:
                    cols.append(extra)

        # Drop any columns that may not exist due to upstream input variability
        cols = [c for c in cols if c in processed_chunk.columns]
        processed_chunk = processed_chunk[cols]

        if first_chunk:
            processed_chunk.to_csv(str(outp), sep="\t", index=False, mode="w")
            first_chunk = False
        else:
            processed_chunk.to_csv(str(outp), sep="\t", index=False, mode="a", header=False)

        log.info("[%s] Finished chunk %d/%d",
                 datetime.now().strftime("%Y-%m-%d %H:%M:%S"), idx, total_chunks)

# ----------------- CLI adapter -----------------
def motif_count_cmd(args):
    """
    CLI entry for `trias motif-count`.
    """
    motif_count(
        input_path=args.input,
        output_path=args.output,
        max_workers=args.max_workers,
        group_chunk_size=args.group_chunk_size
    )
