# src/trias/merge/tsvio.py
"""Streaming reader/writer for TR allele tables.

Deliberately plain stdlib rather than pandas: merging is a line-oriented set
operation with no numeric work, and a 384 MB chr1 table costs ~2 GB as pandas
objects. Files are read in binary so that block offsets stay meaningful.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

log = logging.getLogger("trias.merge")

LOCUS_COLS = ("CHROM", "START", "END", "MOTIF")
ALLELE_COL = "TR_allele"
SAMPLES_COL = "Samples"
MC_COLS = ("max_consec_dna", "total_repeat_dna")

LocusKey = Tuple[str, str, str, str]
Row = List[str]


class TableFormatError(ValueError):
    """Raised when an input table does not look like a TRiAs allele table."""


@dataclass(frozen=True)
class Schema:
    """Column layout of one allele table.

    `vcf-alleles` writes GROUP_ORDER first, then the locus columns, then
    Samples last (extract.py:428), so the columns before CHROM *are* the group
    order and the first of them is the reference group (extract.py:70-72).
    """

    path: str
    columns: Tuple[str, ...]
    group_cols: Tuple[str, ...]
    extra_cols: Tuple[str, ...]     # mc columns, in the order this file has them

    @property
    def ref_group(self) -> str:
        return self.group_cols[0]

    @property
    def level(self) -> str:
        return "mc" if self.extra_cols else "alleles"

    def index(self, name: str) -> int:
        return self.columns.index(name)


def parse_header(path: str, header_line: str) -> Schema:
    cols = header_line.rstrip("\r\n").split("\t")
    for required in (*LOCUS_COLS, ALLELE_COL, SAMPLES_COL):
        if required not in cols:
            raise TableFormatError(f"{path}: missing required column {required!r}")
    i_chrom = cols.index("CHROM")
    if i_chrom == 0:
        raise TableFormatError(f"{path}: no group count columns before CHROM")
    group_cols = tuple(cols[:i_chrom])
    if cols[i_chrom : i_chrom + 4] != list(LOCUS_COLS):
        raise TableFormatError(
            f"{path}: expected {'/'.join(LOCUS_COLS)} to be adjacent, got {cols[i_chrom:i_chrom+4]}"
        )
    if cols[-1] != SAMPLES_COL:
        raise TableFormatError(f"{path}: {SAMPLES_COL} must be the last column")
    tail = cols[cols.index(ALLELE_COL) + 1 : -1]
    for c in tail:
        if c not in MC_COLS:
            raise TableFormatError(f"{path}: unexpected column {c!r} between TR_allele and Samples")
    return Schema(path=path, columns=tuple(cols), group_cols=group_cols, extra_cols=tuple(tail))


def read_schema(path: str) -> Schema:
    with open(path, "rb") as fh:
        first = fh.readline()
    if not first:
        raise TableFormatError(f"{path}: file is empty")
    return parse_header(path, first.decode("utf-8"))


def _locus_key(schema: Schema, row: Row) -> LocusKey:
    i = schema.index("CHROM")
    return (row[i], row[i + 1], row[i + 2], row[i + 3])


def check_row(schema: Schema, row: Row, path: str, lineno: int) -> None:
    """Guard the invariant every valid table satisfies.

    Within a row the Samples tokens are exactly the concatenation of each
    group's tokens, so their count must equal the sum of the group columns
    (extract.py:326-344). Verified to hold on every row of the real chrM and
    chr21 outputs. A violation means a truncated or hand-edited file, and
    silently merging one would corrupt the counts.
    """
    if len(row) != len(schema.columns):
        raise TableFormatError(
            f"{path}:{lineno}: expected {len(schema.columns)} columns, found {len(row)}"
        )
    total = 0
    for gi in range(len(schema.group_cols)):
        try:
            total += int(row[gi])
        except ValueError:
            raise TableFormatError(
                f"{path}:{lineno}: group column {schema.group_cols[gi]!r} is not an integer: {row[gi]!r}"
            ) from None
    cell = row[-1].strip()
    n_tokens = 0 if not cell or cell == "." else cell.count(",") + 1
    if n_tokens != total:
        raise TableFormatError(
            f"{path}:{lineno}: {n_tokens} sample tokens but group counts sum to {total}"
        )


def iter_blocks(path: str, validate: bool = True) -> Iterator[Tuple[LocusKey, List[Row]]]:
    """Yield (locus key, rows) for each run of consecutive rows sharing a locus.

    Both `vcf-alleles` (sorted by BED row index) and `motif-count` (grouped by
    locus) keep a locus' rows contiguous, so no buffering beyond one locus is
    needed.
    """
    schema = read_schema(path)
    key: Optional[LocusKey] = None
    block: List[Row] = []
    with open(path, "rb") as fh:
        fh.readline()  # header
        lineno = 1
        for raw in fh:
            lineno += 1
            line = raw.decode("utf-8").rstrip("\r\n")
            if not line:
                continue
            row = line.split("\t")
            if validate:
                check_row(schema, row, path, lineno)
            k = _locus_key(schema, row)
            if key is None:
                key, block = k, [row]
            elif k == key:
                block.append(row)
            else:
                yield key, block
                key, block = k, [row]
    if key is not None:
        yield key, block


def build_block_index(path: str) -> Tuple[Schema, Dict[LocusKey, Tuple[int, int]]]:
    """Map each locus to the (offset, length) of its block in the file.

    Only keys and offsets are held, never the rows, so indexing a 384 MB table
    costs tens of MB instead of gigabytes. Blocks are re-read on demand.
    """
    schema = read_schema(path)
    index: Dict[LocusKey, Tuple[int, int]] = {}
    key: Optional[LocusKey] = None
    start = 0
    with open(path, "rb") as fh:
        fh.readline()
        pos = fh.tell()
        lineno = 1
        for raw in iter(fh.readline, b""):
            lineno += 1
            nxt = pos + len(raw)
            line = raw.decode("utf-8").rstrip("\r\n")
            if not line:
                pos = nxt
                continue
            row = line.split("\t")
            check_row(schema, row, path, lineno)
            k = _locus_key(schema, row)
            if k != key:
                if key is not None:
                    index[key] = (start, pos - start)
                if k in index:
                    raise TableFormatError(
                        f"{path}:{lineno}: locus {k} appears in more than one place; "
                        "rows of a locus must be contiguous"
                    )
                key, start = k, pos
            pos = nxt
        if key is not None:
            index[key] = (start, pos - start)
    return schema, index


def read_block_at(fh, offset: int, length: int) -> List[Row]:
    fh.seek(offset)
    chunk = fh.read(length).decode("utf-8")
    return [ln.split("\t") for ln in chunk.split("\n") if ln.strip()]


def partition_by_counts(
    schema: Schema, row: Row, path: str = "", lineno: int = 0
) -> List[Tuple[str, List[str]]]:
    """Split a row's Samples cell into per-group token lists.

    The cell lists each group's tokens in header order (extract.py:326-331),
    and each group column holds that group's token count, so slicing the cell
    by the counts recovers which sample belongs to which group -- no group file
    required.
    """
    cell = row[-1].strip()
    tokens = [] if not cell or cell == "." else cell.split(",")
    out: List[Tuple[str, List[str]]] = []
    pos = 0
    for gi, gname in enumerate(schema.group_cols):
        n = int(row[gi])
        out.append((gname, tokens[pos : pos + n]))
        pos += n
    if pos != len(tokens):
        raise TableFormatError(
            f"{path or schema.path}:{lineno}: group counts consume {pos} of {len(tokens)} tokens"
        )
    return out


def first_block_roster(path: str) -> Tuple[Schema, Dict[str, str]]:
    """Read one locus and return its complete token -> group mapping.

    Every locus in a valid run covers every sample: a run whose counts do not
    add up drops the locus entirely (extract.py:319-323), so one locus is the
    whole roster. Callers still union over all of a run's files, which catches
    a truncated file.
    """
    from ..common.loci import strip_plus

    schema = read_schema(path)
    roster: Dict[str, str] = {}
    for _key, rows in iter_blocks(path):
        for row in rows:
            for gname, tokens in partition_by_counts(schema, row):
                for tok in tokens:
                    name, _weak = strip_plus(tok)
                    if name:
                        roster[name] = gname
        break
    return schema, roster


def build_header(group_cols: Sequence[str], extra_cols: Sequence[str]) -> List[str]:
    return [*group_cols, *LOCUS_COLS, ALLELE_COL, *extra_cols, SAMPLES_COL]


class TableWriter:
    """Write a merged table to a temporary path and publish it atomically."""

    def __init__(self, final_path: str, header: Sequence[str], tmp_suffix: str = ".partial"):
        self.final_path = final_path
        self.tmp_path = final_path + tmp_suffix
        parent = os.path.dirname(final_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        self._fh = open(self.tmp_path, "w", encoding="utf-8", newline="\n")
        self._fh.write("\t".join(header) + "\n")
        self.n_rows = 0

    def write(self, row: Sequence[str]) -> None:
        self._fh.write("\t".join(row) + "\n")
        self.n_rows += 1

    def commit(self) -> None:
        self._fh.close()
        os.replace(self.tmp_path, self.final_path)

    def abort(self) -> None:
        try:
            self._fh.close()
        finally:
            if os.path.exists(self.tmp_path):
                os.remove(self.tmp_path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self.commit()
        else:
            self.abort()
        return False
