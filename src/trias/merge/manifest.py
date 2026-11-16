# src/trias/merge/manifest.py
"""Sidecar files describing how an allele table was produced.

Written next to the output but named from its *stem*, so that a sidecar can
never be picked up by a `TR_population_chr*.tsv` glob and fed to the next
stage by accident:

    TR_population_chr1.tsv
    TR_population_chr1.manifest.json
    TR_population_chr1.loci_status.tsv.gz

Everything here is optional. `merge-alleles` reads a manifest when one is
present and derives the same facts from the table itself when it is not, so
tables produced before 0.2.0 merge unchanged.
"""
from __future__ import annotations

import gzip
import json
import logging
import os
import sys
import time
from typing import Dict, Iterable, List, Optional, Sequence

from .. import __version__

log = logging.getLogger("trias.merge")

SCHEMA_VERSION = 1

# Why a locus is absent from the output. "no_variants" means the VCF had no
# record in the interval, i.e. every sample is reference. The others mean the
# locus was attempted and lost, which back-filling cannot tell apart on its own.
STATUS_NO_VARIANTS = "no_variants"
STATUS_CONSENSUS_FAILED = "consensus_failed"
STATUS_BCFTOOLS_FAILED = "bcftools_failed"
STATUS_COUNT_MISMATCH = "count_mismatch"
STATUS_NO_REF_GROUP = "no_ref_group"
STATUS_OK = "ok"


def stem_path(output: str, suffix: str) -> str:
    root, _ext = os.path.splitext(os.path.expanduser(output))
    return root + suffix


def manifest_path(output: str) -> str:
    return stem_path(output, ".manifest.json")


def loci_status_path(output: str) -> str:
    return stem_path(output, ".loci_status.tsv.gz")


def write_manifest(output: str, payload: Dict) -> str:
    path = manifest_path(output)
    payload = dict(payload)
    payload.setdefault("schema_version", SCHEMA_VERSION)
    payload.setdefault("tool", "trias")
    payload.setdefault("trias_version", __version__)
    payload.setdefault("created_utc", time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    payload.setdefault("command", list(sys.argv))
    tmp = path + ".partial"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=False)
        fh.write("\n")
    os.replace(tmp, path)
    return path


def read_manifest(output: str) -> Optional[Dict]:
    path = manifest_path(output)
    if not os.path.exists(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError) as exc:
        log.warning("Ignoring unreadable manifest %s: %s", path, exc)
        return None


def write_loci_status(output: str, rows: Iterable[Sequence[str]]) -> Optional[str]:
    """Record the loci that produced no output, and why."""
    rows = list(rows)
    if not rows:
        return None
    path = loci_status_path(output)
    tmp = path + ".partial"
    with gzip.open(tmp, "wt", encoding="utf-8", newline="\n") as fh:
        fh.write("CHROM\tSTART\tEND\tMOTIF\tstatus\n")
        for row in rows:
            fh.write("\t".join(str(c) for c in row) + "\n")
    os.replace(tmp, path)
    return path


def read_loci_status(output: str) -> Dict[tuple, str]:
    path = loci_status_path(output)
    if not os.path.exists(path):
        return {}
    out: Dict[tuple, str] = {}
    try:
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            fh.readline()
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) == 5:
                    out[tuple(parts[:4])] = parts[4]
    except OSError as exc:
        log.warning("Ignoring unreadable locus status file %s: %s", path, exc)
        return {}
    return out


def build_group_payload(group_files: Sequence[str], sample_groups: Dict[str, tuple], group_order: Sequence[str]) -> List[Dict]:
    """Describe each --groups file and the samples it contributed."""
    by_group: Dict[str, List[Dict]] = {g: [] for g in group_order}
    for base, (group, ploidy) in sample_groups.items():
        if group in by_group:
            by_group[group].append({"name": base, "ploidy": int(ploidy)})
    out = []
    for i, (path, group) in enumerate(zip(group_files, group_order)):
        out.append({
            "path": os.path.abspath(os.path.expanduser(path)),
            "group": group,
            "is_ref": i == 0,
            "samples": sorted(by_group.get(group, []), key=lambda s: s["name"]),
        })
    return out


def file_stat(path: Optional[str]) -> Optional[Dict]:
    if not path:
        return None
    p = os.path.expanduser(path)
    info: Dict = {"path": os.path.abspath(p)}
    try:
        st = os.stat(p)
        info["bytes"] = st.st_size
        info["mtime"] = int(st.st_mtime)
    except OSError:
        pass
    return info
