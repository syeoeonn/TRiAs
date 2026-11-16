# src/trias/common/loci.py
"""Locus and chromosome helpers shared by `dntre` and `merge-alleles`.

These used to live in dntre/dntre.py; they are here so more than one
subcommand can use them without importing the dnTRE caller.
"""
from __future__ import annotations

import os
import re
from typing import List, Tuple

# The 4-tuple that identifies a TR locus everywhere in the pipeline. Compared
# as strings, never as numbers, so that "10000" and "10000.0" cannot diverge.
LOCUS_KEY_COLS: Tuple[str, ...] = ("CHROM", "START", "END", "MOTIF")


def chr_key(chr_value: str) -> int:
    """Sort key that puts chr2 before chr10 and the sex/MT contigs last."""
    v = str(chr_value).replace("chr", "").upper()
    if v == "X":
        return 23
    if v == "Y":
        return 24
    if v in ("M", "MT"):
        return 25
    try:
        return int(v)
    except Exception:
        return 1000


CHR_RE = re.compile(
    r"(?:^|[^A-Za-z0-9])chr(?P<chr>[0-9]+|X|Y|M|MT)(?:[^A-Za-z0-9]|$)",
    re.IGNORECASE,
)


def extract_chr_from_path(path: str) -> str:
    """Derive the chromosome from a file NAME (not from its contents).

    `dntre` relies on this, which is why merged files must keep the basename
    of the file they were merged from.
    """
    m = CHR_RE.search(os.path.basename(path))
    if not m:
        return ""
    raw = m.group("chr").upper()
    if raw == "MT":
        raw = "M"
    return "chr" + raw


def locus_sort_key(chrom: str, start: str, end: str, motif: str):
    """Order loci the way bed-prepare emits them.

    START > END marks a negative-strand locus, so the span is normalised with
    min/max; the stored coordinate strings themselves are never rewritten.
    """
    try:
        a, b = int(start), int(end)
    except ValueError:
        a = b = 0
    lo, hi = (a, b) if a <= b else (b, a)
    return (chr_key(chrom), lo, hi, motif)


def strip_plus(token: str) -> Tuple[str, bool]:
    """Split a Samples token into (name, was_marked_missing).

    `+HG002.1+` means at least one variant row in the interval had GT '.' for
    that haplotype, so the call is low-confidence. Unlike dntre.samples_to_set,
    this keeps that fact instead of discarding it -- merging needs it to tell a
    real call from an assumed one.
    """
    t = token.strip()
    if len(t) >= 3 and t.startswith("+") and t.endswith("+"):
        return t[1:-1], True
    return t, False


def split_samples(samples_cell: str) -> List[Tuple[str, bool]]:
    """Parse a Samples cell into [(token, was_marked_missing), ...]."""
    s = (samples_cell or "").strip()
    if not s or s == ".":
        return []
    out = []
    for raw in s.split(","):
        name, weak = strip_plus(raw)
        if name:
            out.append((name, weak))
    return out


def format_token(name: str, weak: bool) -> str:
    """Inverse of strip_plus."""
    return f"+{name}+" if weak else name


_MOTIF_RUN_RE = re.compile(r"\(([A-Za-z]+)\)(\d*)")


def raw_allele_length(allele: str) -> int:
    """Length of the ORIGINAL sequence behind a (possibly compressed) allele.

    `motif-count` rewrites TR_allele in place as e.g. `TCCC(CTCCCA)2TACT`, so
    len() of a stored mc allele is not the sequence length. vcf-alleles sorts
    the rows of a locus by sequence length, and a merged table has to reproduce
    that order at both levels, so the compressed form is expanded here.
    Uncompressed alleles contain no parentheses and fall through unchanged.
    """
    s = str(allele)
    if "(" not in s:
        return len(s)
    total = 0
    pos = 0
    for m in _MOTIF_RUN_RE.finditer(s):
        total += m.start() - pos
        n = m.group(2)
        total += len(m.group(1)) * (int(n) if n else 1)
        pos = m.end()
    total += len(s) - pos
    return total


def max_motif_count(value) -> int:
    """Largest per-motif number in a motif-count cell such as '2_12'; -1 if none."""
    best = -1
    for part in str(value).split("_"):
        part = part.strip()
        if part.isdigit():
            best = max(best, int(part))
    return best


def allele_order_key(compressed_allele: str, max_consec: str, total_repeat: str) -> Tuple[int, int, int, int]:
    """Sort key for the rows of one locus; every component ascends.

    1. length of the original sequence
    2. largest max_consec_dna across motifs
    3. largest total_repeat_dna across motifs
    4. length of the compressed allele string

    Rows equal on all four keep the order they arrived in, so callers must
    sort stably. `motif-count` and `merge-alleles` share this key, which is
    what keeps a merged table in the same row order as a freshly counted one.
    """
    allele = str(compressed_allele)
    return (
        raw_allele_length(allele),
        max_motif_count(max_consec),
        max_motif_count(total_repeat),
        len(allele),
    )
