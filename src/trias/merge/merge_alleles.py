# src/trias/merge/merge_alleles.py
"""Merge TR allele tables from separate `vcf-alleles` runs into one cohort.

The unit of work is a *token*: one haplotype of one sample, as it appears in
the Samples column (`HG00438.1`, or a bare name for a ploidy-1 sample).

Two properties drive the design:

* Group counts are recomputed from the union of tokens rather than added up.
  The reference pseudo-sample appears in every run, so adding would count it
  twice; a union keeps it at one and makes re-merging a run a no-op.
* A locus missing from a run means the VCF had no variant there, i.e. every
  sample of that run is REF, so those tokens are back-filled onto the locus'
  reference row. That preserves the invariant that a locus covers every
  sample, which `motif-summary` relies on for its allele counts.
"""
from __future__ import annotations

import glob
import logging
import os
import shutil
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Sequence, Tuple

from ..common.loci import allele_order_key, extract_chr_from_path, format_token, locus_sort_key, strip_plus
from . import tsvio
from .tsvio import LocusKey, Row, Schema, TableFormatError

log = logging.getLogger("trias.merge")


# --------------------------------------------------------------------------
# row order within a locus
# --------------------------------------------------------------------------

def _compress_like_motif_count(allele: str, motif: str) -> Tuple[str, str, str]:
    """(compressed allele, max_consec_dna, total_repeat_dna) exactly as motif-count writes them."""
    from ..motif.count import compress_with_multi_motif_separated, format_motif_separated_stats

    if allele != "." and motif not in (".", "nan"):
        comp, _total, _max, stats = compress_with_multi_motif_separated(allele, motif)
        total_sep, max_sep = format_motif_separated_stats(motif, stats)
        return comp, max_sep, total_sep
    return allele, ".", "."


def _order_keys(
    alleles: Sequence[str], values: Dict[str, Dict[str, str]], motif: str, is_mc: bool
) -> Dict[str, Tuple[int, int, int, int]]:
    """Map each allele of a locus to the shared row-order key.

    At the motif-count level the counts are already in the table. An
    uncompressed table has none, but only alleles that tie on length need them,
    so just those are compressed -- the same computation motif-count would do,
    which keeps both levels in the same order.
    """
    if is_mc:
        return {
            a: allele_order_key(
                a, values[a].get("max_consec_dna", "."), values[a].get("total_repeat_dna", ".")
            )
            for a in alleles
        }
    lengths = Counter(len(a) for a in alleles)
    keys: Dict[str, Tuple[int, int, int, int]] = {}
    for a in alleles:
        if lengths[len(a)] < 2:
            keys[a] = (len(a), -1, -1, len(a))
        else:
            keys[a] = allele_order_key(*_compress_like_motif_count(a, motif))
    return keys

CONFLICT_COLUMNS = [
    "CHROM", "START", "END", "MOTIF", "Sample", "Group", "Kind",
    "Base_allele", "Other_allele", "Other_source", "Chosen",
]

REPORT_COLUMNS = [
    "CHROM", "base_loci", "add_loci", "merged_loci", "backfilled_loci",
    "backfilled_tokens", "resolved_missing_gt", "new_alleles", "conflicts",
]


class MergeError(RuntimeError):
    """Inputs cannot be merged, or a conflict aborted the run."""


# --------------------------------------------------------------------------
# run discovery
# --------------------------------------------------------------------------

class RunView:
    """One `vcf-alleles`/`motif-count` run: its files, schema and roster."""

    def __init__(self, role: str, pattern: str):
        self.role = role
        self.pattern = pattern
        self.paths: List[str] = sorted(glob.glob(os.path.expanduser(pattern)))
        if not self.paths:
            raise MergeError(f"{role}: no files matched {pattern!r}")
        self.by_chrom: Dict[str, str] = {}
        for p in self.paths:
            chrom = extract_chr_from_path(p) or os.path.basename(p)
            if chrom in self.by_chrom:
                raise MergeError(
                    f"{role}: {chrom} matched by two files, {self.by_chrom[chrom]} and {p}"
                )
            self.by_chrom[chrom] = p
        self.schema: Schema = tsvio.read_schema(self.paths[0])
        self.roster: Dict[str, str] = self._build_roster()

    def _build_roster(self) -> Dict[str, str]:
        """Union the first locus of every file; each already holds the full roster."""
        roster: Dict[str, str] = {}
        for p in self.paths:
            schema, part = tsvio.first_block_roster(p)
            if schema.group_cols != self.schema.group_cols:
                raise MergeError(
                    f"{self.role}: {p} has group columns {schema.group_cols} but "
                    f"{self.paths[0]} has {self.schema.group_cols}"
                )
            if schema.extra_cols != self.schema.extra_cols:
                raise MergeError(
                    f"{self.role}: {p} is a {schema.level} table but "
                    f"{self.paths[0]} is a {self.schema.level} table"
                )
            for tok, grp in part.items():
                prev = roster.get(tok)
                if prev is not None and prev != grp:
                    raise MergeError(
                        f"{self.role}: sample {tok} is in group {prev} in one file and {grp} in {p}"
                    )
                roster[tok] = grp
            if roster and part and set(part) != set(roster):
                missing = sorted(set(roster) - set(part))[:5]
                extra = sorted(set(part) - set(roster))[:5]
                log.warning(
                    "%s: %s covers a different sample set (missing %s, extra %s); "
                    "files may come from different runs",
                    self.role, p, missing, extra,
                )
        if not roster:
            raise MergeError(f"{self.role}: no samples found; are the tables empty?")
        return roster


def _plan_runs(base_pattern: str, add_patterns: Sequence[str], require_disjoint: bool):
    runs = [RunView("base", base_pattern)]
    for i, pat in enumerate(add_patterns, start=1):
        runs.append(RunView(f"add{i}", pat))

    base = runs[0]
    for r in runs[1:]:
        if r.schema.ref_group != base.schema.ref_group:
            raise MergeError(
                f"{r.role}: reference group is {r.schema.ref_group!r} but base uses "
                f"{base.schema.ref_group!r}; the runs are not on the same reference"
            )
        if r.schema.level != base.schema.level:
            raise MergeError(
                f"{r.role}: is a {r.schema.level} table but base is a {base.schema.level} "
                "table; merge tables from the same pipeline stage"
            )

    roster: Dict[str, str] = {}
    for r in runs:
        for tok, grp in r.roster.items():
            prev = roster.get(tok)
            if prev is not None and prev != grp:
                raise MergeError(
                    f"sample {tok} is in group {prev} in one run and {grp} in {r.role}"
                )
            roster[tok] = grp

    if require_disjoint:
        ref_tok = base.schema.ref_group
        seen: Dict[str, str] = {}
        for r in runs:
            for tok in r.roster:
                if tok == ref_tok:
                    continue
                if tok in seen:
                    raise MergeError(
                        f"--require-disjoint-samples: {tok} appears in both {seen[tok]} and {r.role}"
                    )
                seen[tok] = r.role

    groups: List[str] = []
    for r in runs:
        for g in r.schema.group_cols:
            if g not in groups:
                groups.append(g)

    extra: List[str] = []
    for r in runs:
        for c in r.schema.extra_cols:
            if c not in extra:
                extra.append(c)

    return runs, roster, groups, extra


# --------------------------------------------------------------------------
# per-locus merge
# --------------------------------------------------------------------------

class _Conflict(Exception):
    def __init__(self, row: List[str]):
        super().__init__("; ".join(row))
        self.row = row


def _merge_locus(
    key: LocusKey,
    blocks: List[Tuple[str, Schema, List[Row]]],
    roster: Dict[str, str],
    groups: Sequence[str],
    extra_cols: Sequence[str],
    opts: dict,
    conflicts: List[List[str]],
    stats: Dict[str, int],
) -> List[List[str]]:
    chrom, start, end, motif = key
    ref_group = blocks[0][1].ref_group

    alleles: Dict[str, Dict[str, str]] = {}          # allele -> {mc col: value}
    assign: Dict[str, Tuple[str, bool, str]] = {}    # token -> (allele, weak, source)
    ref_allele: Optional[str] = None

    for source, schema, rows in blocks:
        i_allele = schema.index(tsvio.ALLELE_COL)
        for row in rows:
            allele = row[i_allele]
            values = {c: row[schema.index(c)] for c in schema.extra_cols}
            known = alleles.get(allele)
            if known is None:
                alleles[allele] = values
                if source != blocks[0][0]:
                    stats["new_alleles"] += 1
            else:
                for c, v in values.items():
                    if c in known and known[c] != v:
                        raise MergeError(
                            f"{chrom}:{start}-{end} {motif}: allele has {c}={known[c]!r} in "
                            f"{blocks[0][0]} but {v!r} in {source}; the tables come from "
                            "different motif-count runs"
                        )
                known.update(values)

            for gname, tokens in tsvio.partition_by_counts(schema, row):
                for disp in tokens:
                    name, weak = strip_plus(disp)
                    if gname == ref_group and name == ref_group:
                        if ref_allele is None:
                            ref_allele = allele
                        elif ref_allele != allele:
                            raise MergeError(
                                f"{chrom}:{start}-{end} {motif}: the reference allele differs "
                                f"between runs ({blocks[0][0]} vs {source}); the runs use a "
                                "different reference genome or a different BED"
                            )
                    prev = assign.get(name)
                    if prev is None:
                        assign[name] = (allele, weak, source)
                        continue
                    p_allele, p_weak, p_source = prev
                    if p_allele == allele:
                        assign[name] = (allele, p_weak and weak, p_source)
                        continue
                    if p_weak != weak:
                        # One side is a real call, the other only assumed REF
                        # because a GT was '.'; the real call wins.
                        stats["resolved_missing_gt"] += 1
                        conflicts.append([
                            chrom, start, end, motif, name, roster.get(name, "?"),
                            "resolved_missing_gt", p_allele, allele, source,
                            "add" if p_weak else "base",
                        ])
                        if p_weak:
                            assign[name] = (allele, weak, source)
                        continue
                    policy = opts["on_conflict"]
                    chosen = "abort" if policy == "error" else policy
                    crow = [
                        chrom, start, end, motif, name, roster.get(name, "?"),
                        "allele_disagreement", p_allele, allele, source, chosen,
                    ]
                    stats["conflicts"] += 1
                    if policy == "error":
                        raise _Conflict(crow)
                    conflicts.append(crow)
                    if policy == "add":
                        assign[name] = (allele, weak, source)

    missing = [t for t in roster if t not in assign]
    if missing:
        if opts["fill_missing"] == "none":
            return []
        if ref_allele is None:
            raise MergeError(
                f"{chrom}:{start}-{end} {motif}: no reference row found (no {ref_group!r} "
                "token), so the samples missing from this locus cannot be filled in"
            )
        weak = opts["backfill_mark"] == "plus"
        for name in missing:
            assign[name] = (ref_allele, weak, "backfill")
        stats["backfilled_loci"] += 1
        stats["backfilled_tokens"] += len(missing)

    by_allele: Dict[str, Dict[str, List[str]]] = defaultdict(lambda: defaultdict(list))
    for name, (allele, weak, _src) in assign.items():
        by_allele[allele][roster[name]].append(format_token(name, weak))

    ordered = [a for a in alleles if a in by_allele]
    keys = _order_keys(ordered, alleles, motif, bool(extra_cols))
    ordered.sort(key=keys.__getitem__)  # stable: full ties keep first-seen order

    out: List[List[str]] = []
    for allele in ordered:
        counts: List[str] = []
        tokens: List[str] = []
        for g in groups:
            part = sorted(by_allele[allele].get(g, []))
            counts.append(str(len(part)))
            tokens.extend(part)
        values = alleles[allele]
        out.append([
            *counts, chrom, start, end, motif, allele,
            *(values.get(c, ".") for c in extra_cols),
            ",".join(tokens),
        ])

    total = sum(int(c) for row in out for c in row[: len(groups)])
    if total != len(roster):
        raise MergeError(
            f"{chrom}:{start}-{end} {motif}: merged locus covers {total} haplotypes but the "
            f"cohort has {len(roster)}"
        )
    return out


# --------------------------------------------------------------------------
# per-chromosome worker
# --------------------------------------------------------------------------

def _merge_chromosome(task: dict) -> dict:
    chrom = task["chrom"]
    sources: List[Tuple[str, str]] = task["sources"]  # [(role, path), ...]
    roster: Dict[str, str] = task["roster"]
    groups: List[str] = task["groups"]
    extra_cols: List[str] = task["extra_cols"]
    opts: dict = task["opts"]

    schemas: Dict[str, Schema] = {}
    indexes: Dict[str, Dict[LocusKey, Tuple[int, int]]] = {}
    for role, path in sources:
        schemas[role], indexes[role] = tsvio.build_block_index(path)

    all_keys = set()
    for idx in indexes.values():
        all_keys.update(idx)
    ordered_keys = sorted(all_keys, key=lambda k: locus_sort_key(*k))

    stats = {k: 0 for k in (
        "merged_loci", "backfilled_loci", "backfilled_tokens",
        "resolved_missing_gt", "new_alleles", "conflicts",
    )}
    stats["base_loci"] = len(indexes[sources[0][0]])
    stats["add_loci"] = sum(len(indexes[role]) for role, _ in sources[1:])

    conflicts: List[List[str]] = []
    header = tsvio.build_header(groups, extra_cols)
    handles = {role: open(path, "rb") for role, path in sources}
    writer = tsvio.TableWriter(task["out_path"], header)
    try:
        for key in ordered_keys:
            blocks: List[Tuple[str, Schema, List[Row]]] = []
            for role, _path in sources:
                loc = indexes[role].get(key)
                if loc is not None:
                    blocks.append((role, schemas[role], tsvio.read_block_at(handles[role], *loc)))
            rows = _merge_locus(
                key, blocks, roster, groups, extra_cols, opts, conflicts, stats
            )
            if rows:
                stats["merged_loci"] += 1
                for row in rows:
                    writer.write(row)
    except _Conflict as exc:
        writer.abort()
        return {"chrom": chrom, "stats": stats, "conflicts": [exc.row], "aborted": str(exc)}
    except BaseException:
        writer.abort()
        raise
    else:
        writer.commit()
    finally:
        for fh in handles.values():
            fh.close()

    return {"chrom": chrom, "stats": stats, "conflicts": conflicts, "aborted": None}


# --------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------

def merge_alleles(
    base: str,
    add: Sequence[str],
    outdir: str,
    fill_missing: str = "ref",
    backfill_mark: str = "plain",
    on_conflict: str = "error",
    require_disjoint_samples: bool = False,
    jobs: Optional[int] = None,
    report: Optional[str] = None,
) -> None:
    runs, roster, groups, extra_cols = _plan_runs(base, add, require_disjoint_samples)
    log.info("Merging %d runs: %s", len(runs), ", ".join(f"{r.role}={r.pattern}" for r in runs))
    log.info("Reference group: %s", runs[0].schema.ref_group)
    log.info("Level: %s", runs[0].schema.level)
    log.info("Merged groups: %s", ", ".join(groups))
    log.info("Cohort size: %d haplotypes across %d groups", len(roster), len(groups))

    chroms = sorted(
        {c for r in runs for c in r.by_chrom},
        key=lambda c: locus_sort_key(c, "0", "0", ""),
    )
    for r in runs[1:]:
        absent = [c for c in runs[0].by_chrom if c not in r.by_chrom]
        if absent:
            log.warning(
                "%s has no file for %s; every locus there will be filled in as reference",
                r.role, ", ".join(absent),
            )

    outdir = os.path.expanduser(outdir)
    tmpdir = os.path.join(outdir, ".tmp")
    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)
    os.makedirs(tmpdir, exist_ok=True)

    opts = {
        "fill_missing": fill_missing,
        "backfill_mark": backfill_mark,
        "on_conflict": on_conflict,
    }

    tasks = []
    for chrom in chroms:
        sources = [(r.role, r.by_chrom[chrom]) for r in runs if chrom in r.by_chrom]
        basename = os.path.basename(sources[0][1])
        tasks.append({
            "chrom": chrom,
            "sources": sources,
            "out_path": os.path.join(tmpdir, basename),
            "roster": roster,
            "groups": groups,
            "extra_cols": extra_cols,
            "opts": opts,
        })
    # Largest first, so the slowest chromosome is never the last one started.
    tasks.sort(key=lambda t: os.path.getsize(t["sources"][0][1]), reverse=True)

    if jobs is None:
        jobs = min(8, os.cpu_count() or 1, len(tasks))
    jobs = max(1, min(jobs, len(tasks)))

    results: List[dict] = []
    try:
        if jobs == 1:
            for t in tasks:
                results.append(_merge_chromosome(t))
        else:
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                futs = {ex.submit(_merge_chromosome, t): t["chrom"] for t in tasks}
                for fut in as_completed(futs):
                    results.append(fut.result())
        aborted = [r for r in results if r["aborted"]]
        if aborted:
            _write_conflicts(os.path.join(outdir, "merge_conflicts.tsv"),
                             [row for r in results for row in r["conflicts"]])
            raise MergeError(
                f"aborted on a sample conflict ({aborted[0]['aborted']}). "
                f"See {os.path.join(outdir, 'merge_conflicts.tsv')}. "
                "Use --on-conflict base or add to keep going."
            )
    except BaseException:
        shutil.rmtree(tmpdir, ignore_errors=True)
        raise

    for name in sorted(os.listdir(tmpdir)):
        os.replace(os.path.join(tmpdir, name), os.path.join(outdir, name))
    shutil.rmtree(tmpdir, ignore_errors=True)

    results.sort(key=lambda r: locus_sort_key(r["chrom"], "0", "0", ""))
    _write_conflicts(os.path.join(outdir, "merge_conflicts.tsv"),
                     [row for r in results for row in r["conflicts"]])
    _write_report(report or os.path.join(outdir, "merge_report.tsv"), results)

    total = {k: sum(r["stats"].get(k, 0) for r in results) for k in REPORT_COLUMNS[1:]}
    log.info(
        "Merged %d loci (%d filled in as reference, %d new alleles, %d conflicts) into %s",
        total["merged_loci"], total["backfilled_loci"], total["new_alleles"],
        total["conflicts"], outdir,
    )


def _write_conflicts(path: str, rows: List[List[str]]) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(CONFLICT_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def _write_report(path: str, results: List[dict]) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(REPORT_COLUMNS) + "\n")
        for r in results:
            fh.write("\t".join([r["chrom"], *(str(r["stats"].get(c, 0)) for c in REPORT_COLUMNS[1:])]) + "\n")


def merge_alleles_cmd(args) -> None:
    """CLI entry for `trias merge-alleles`."""
    try:
        merge_alleles(
            base=args.base,
            add=args.add,
            outdir=args.outdir,
            fill_missing=args.fill_missing,
            backfill_mark=args.backfill_mark,
            on_conflict=args.on_conflict,
            require_disjoint_samples=args.require_disjoint_samples,
            jobs=args.jobs,
            report=args.report,
        )
    except (MergeError, TableFormatError) as exc:
        log.error("%s", exc)
        sys.exit(1)
