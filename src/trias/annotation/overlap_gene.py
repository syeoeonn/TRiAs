#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import csv
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

# Type hints for clarity
Row = List[str]
AnnEntry = Tuple[int, int, str, str, str]  # (start, end, region_type, strand, gene)
AnnDict = Dict[str, List[AnnEntry]]
RegionEntry = Tuple[str, int, int, Row]    # (chrom, start, end, original_row)


def parse_annotation_row(row: Row) -> Tuple[str, int, int, str, str, str]:
    """
    Parse one line from the annotation file.
    Expected columns:
      0: chrom
      1: start (1-based)
      2: end   (1-based, inclusive)
      3: region_type
      4: strand
      5: gene (e.g., PLEKHN1(ENSG...))
      6: optional source column
    """
    chrom = row[0]
    start = int(row[1])
    end = int(row[2])
    rtype = row[3]
    strand = row[4]
    gene = row[5]
    return chrom, start, end, rtype, strand, gene


def load_annotation(path: str) -> AnnDict:
    """
    Load annotation file and build a dictionary: chrom -> list of (start, end, type, strand, gene).
    Tabs and mixed whitespace are supported.
    """
    ann: AnnDict = defaultdict(list)
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row:
                continue
            # Handle lines with only spaces instead of tabs
            if len(row) == 1 and row[0].strip():
                row = row[0].split()
            if len(row) < 6:
                continue
            # Skip headers or malformed rows where start/end are not integers
            try:
                _ = int(row[1])
                _ = int(row[2])
            except Exception:
                continue

            chrom, start, end, rtype, strand, gene = parse_annotation_row(row)
            ann[chrom].append((start, end, rtype, strand, gene))

    # Sort by start position for each chromosome
    for c in ann:
        ann[c].sort(key=lambda x: x[0])
    return ann


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """
    Compute the overlap length between two genomic intervals (1-based, inclusive).
    Returns 0 if no overlap.
    """
    s = max(a_start, b_start)
    e = min(a_end, b_end)
    return max(0, e - s + 1)


def load_regions(path: str) -> Tuple[List[RegionEntry], Optional[Row]]:
    """
    Load regions file.
    - If a header exists, it must contain 'Chrom/Chr', 'Start', and 'End' columns (case-insensitive).
    - If no header exists, the first 3 columns are assumed to be chrom/start/end.

    Returns:
        regions: List of tuples (chrom, start, end, original_row)
        header:  Header row if exists, otherwise None
    """
    regions: List[RegionEntry] = []
    with open(path, newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        colmap = {}
        has_header = False

        if header and any(h.strip().lower() in ("chrom", "chr") for h in header):
            has_header = True
            for i, h in enumerate(header):
                hl = h.strip().lower()
                if hl in ("chrom", "chr"):
                    colmap['chrom'] = i
                elif hl == "start":
                    colmap['start'] = i
                elif hl == "end":
                    colmap['end'] = i
            if not {'chrom', 'start', 'end'} <= set(colmap):
                raise ValueError("The regions file must contain Chrom/Start/End columns in the header.")
        else:
            # Treat the first line as data when no header is detected
            if header:
                row = header
                regions.append((row[0], int(row[1]), int(row[2]), row))
            colmap = {'chrom': 0, 'start': 1, 'end': 2}

        # Parse remaining lines
        for row in reader:
            if not row or len(row) < 3:
                continue
            chrom = row[colmap['chrom']]
            start = int(row[colmap['start']])
            end = int(row[colmap['end']])
            regions.append((chrom, start, end, row))

    return regions, (header if has_header else None)


def annotate_aggregate(regions: List[RegionEntry], ann: AnnDict) -> List[Row]:
    """
    For each region, collect all overlapping annotation entries.
    Append four columns:
      - Ann_Genes
      - Ann_RegionTypes
      - Ann_Strands
      - Ann_Overlap_bps
    Values are joined by ';'. If no overlap exists, all are set to '.'.
    """
    out: List[Row] = []
    for chrom, rs, re, raw in regions:
        if chrom not in ann:
            out.append(raw + [".", ".", ".", "."])
            continue

        overlaps = []
        # Annotation list is pre-sorted by start position
        for (as_, ae, rtype, strand, gene) in ann[chrom]:
            if as_ > re:
                break
            if ae < rs:
                continue
            ov = overlap_len(rs, re, as_, ae)
            if ov > 0:
                overlaps.append((ov, gene, rtype, strand))

        if not overlaps:
            out.append(raw + [".", ".", ".", "."])
        else:
            # Sort by overlap length (descending)
            overlaps.sort(key=lambda x: -x[0])
            genes = ";".join([g for _, g, _, _ in overlaps])
            rtypes = ";".join([t for _, _, t, _ in overlaps])
            strands = ";".join([s for _, _, _, s in overlaps])
            ovbps = ";".join([str(ov) for ov, _, _, _ in overlaps])
            out.append(raw + [genes, rtypes, strands, ovbps])
    return out


def write_output(output_path: str, rows: List[Row], regions_header: Optional[Row]) -> None:
    """
    Write results to output file.
    Output columns = original region columns + 4 new annotation columns.
    """
    with open(output_path, "w", newline="") as f:
        w = csv.writer(f, delimiter='\t')
        if regions_header:
            header_out = list(regions_header) + ["Ann_Genes", "Ann_RegionTypes", "Ann_Strands", "Ann_Overlap_bps"]
            w.writerow(header_out)
        w.writerows(rows)


def annotate_paths(annotation_path: str, regions_path: str, output_path: str) -> None:
    """
    Library entry point: run full annotation pipeline using file paths.
    """
    ann = load_annotation(annotation_path)
    regions, regions_header = load_regions(regions_path)
    out_rows = annotate_aggregate(regions, ann)
    write_output(output_path, out_rows, regions_header)


def annotate_cmd(args) -> None:
    """
    CLI entry point (compatible with other TRiAs subcommands).
    """
    annotate_paths(args.annotation, args.regions, args.output)
