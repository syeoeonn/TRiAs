#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations
import csv
from typing import Dict, List, Tuple, Optional, Iterable, Set

Row = List[str]

# ----------------------------
# Utilities
# ----------------------------

def _normalize_gene_token(token: str) -> str:
    """
    Normalize a gene token possibly containing an Ensembl ID in parentheses,
    e.g., "PLEKHN1(ENSG000001...)" -> "PLEKHN1".
    """
    if not token:
        return ""
    t = token.strip()
    # Split off anything after the first '('
    if "(" in t:
        t = t.split("(", 1)[0]
    return t.strip()


def _split_genes(cell: str) -> List[str]:
    """
    Split a semicolon-joined gene list like "GENE1;GENE2(ENSG...);GENE3"
    and normalize into plain symbols.
    """
    if not cell or cell == ".":
        return []
    parts = [p.strip() for p in cell.split(";") if p.strip()]
    return [_normalize_gene_token(p) for p in parts if p]


def _find_col_index(header: List[str], candidates: Iterable[str]) -> Optional[int]:
    """
    Find the index of a header column among a set of candidate names (case-insensitive).
    Returns None if not found.
    """
    lower = [h.strip().lower() for h in header]
    cand = [c.strip().lower() for c in candidates]
    for i, h in enumerate(lower):
        if h in cand:
            return i
    return None


def _dedup_order_preserving(items: Iterable[str]) -> List[str]:
    """
    Deduplicate while preserving order.
    """
    seen: Set[str] = set()
    out: List[str] = []
    for x in items:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


# ----------------------------
# Mapping loader
# ----------------------------

def load_moi_mapping(path: str) -> Dict[str, Dict[str, List[str]]]:
    """
    Load a gene-to-MOI/disease mapping TSV.

    The function flexibly detects columns by header:
      - Gene symbol column (candidates): gene, gene_symbol, symbol
      - MOI column (candidates): moi, moi_title, mode_of_inheritance
      - Disease column (candidates): disease, diseases, phenotype, disorder, condition, omim_disease

    Returns:
        mapping: dict keyed by UPPERCASE gene symbol ->
                 {"moi": [moi strings...], "disease": [disease strings...]}
    """
    mapping: Dict[str, Dict[str, List[str]]] = {}
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError("MOI mapping file appears to be empty.")

        gene_idx = _find_col_index(header, ("gene", "gene_symbol", "symbol"))
        moi_idx = _find_col_index(header, ("moi", "moi_title", "mode_of_inheritance"))
        disease_idx = _find_col_index(header, ("disease", "diseases", "phenotype", "disorder", "condition", "omim_disease"))

        if gene_idx is None:
            raise ValueError("Could not find a gene symbol column in the mapping file header.")
        if moi_idx is None and disease_idx is None:
            raise ValueError("Could not find any MOI or disease column in the mapping file header.")

        for row in reader:
            if not row or len(row) <= gene_idx:
                continue
            gene = (row[gene_idx] or "").strip()
            if not gene:
                continue
            gene_sym = _normalize_gene_token(gene).upper()

            moi_val = (row[moi_idx].strip() if (moi_idx is not None and moi_idx < len(row)) else "")
            dis_val = (row[disease_idx].strip() if (disease_idx is not None and disease_idx < len(row)) else "")

            if gene_sym not in mapping:
                mapping[gene_sym] = {"moi": [], "disease": []}
            if moi_val:
                mapping[gene_sym]["moi"].append(moi_val)
            if dis_val:
                mapping[gene_sym]["disease"].append(dis_val)

    # Deduplicate lists
    for gene_sym, dd in mapping.items():
        dd["moi"] = _dedup_order_preserving([x for x in dd["moi"] if x])
        dd["disease"] = _dedup_order_preserving([x for x in dd["disease"] if x])

    return mapping


# ----------------------------
# Table I/O
# ----------------------------

def read_table(path: str) -> Tuple[List[Row], Optional[Row]]:
    """
    Read a TSV table. Returns (rows, header).
    If a header is present, it is returned. Otherwise, header is None.
    """
    rows: List[Row] = []
    with open(path, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        # Peek first row
        first = next(reader, None)
        if first is None:
            return [], None

        # Heuristics: if the row contains any non-digit column names like "chrom" or "Ann_Genes",
        # consider it a header. Otherwise, treat it as data and no header.
        is_header = any(tok.strip().lower() in ("chrom", "chr", "ann_genes", "gene", "start", "end") for tok in first)

        header: Optional[Row] = first if is_header else None
        if not is_header:
            rows.append(first)

        for row in reader:
            if row:
                rows.append(row)

    return rows, header


def write_table(path: str, rows: List[Row], header: Optional[Row]) -> None:
    """
    Write a TSV table to the given path, including header if provided.
    """
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        if header:
            w.writerow(header)
        w.writerows(rows)


# ----------------------------
# Core annotation
# ----------------------------

def annotate_moi_rows(
    rows: List[Row],
    header: Optional[Row],
    mapping: Dict[str, Dict[str, List[str]]],
    gene_col_candidates: Tuple[str, ...] = ("Ann_Genes", "Gene", "Genes"),
    out_moi_col: str = "Ann_MOI",
    out_dis_col: str = "Ann_Diseases",
) -> Tuple[List[Row], Optional[Row]]:
    """
    For each row, read the gene column (semicolon-joined symbols),
    map to MOI and disease(s), and append two new columns.

    Args:
        rows:    Input table rows (without header).
        header:  Header row if exists, otherwise None.
        mapping: Dict of {GENE -> {"moi": [...], "disease": [...]} }.
        gene_col_candidates: Candidate names to locate the gene column.
        out_moi_col: Name of the output MOI column to append.
        out_dis_col: Name of the output Disease column to append.

    Returns:
        (out_rows, out_header)
    """
    gene_col_idx: Optional[int] = None
    out_header: Optional[Row] = None

    if header:
        gene_col_idx = _find_col_index(header, gene_col_candidates)
        if gene_col_idx is None:
            raise ValueError(
                f"Could not find a gene column among candidates: {', '.join(gene_col_candidates)}"
            )
        out_header = list(header) + [out_moi_col, out_dis_col]

    out_rows: List[Row] = []
    for row in rows:
        # If no header, assume genes are in the last four columns appended by overlap:
        # [..., Ann_Genes, Ann_RegionTypes, Ann_Strands, Ann_Overlap_bps]
        if gene_col_idx is None:
            if len(row) < 4:
                # Not enough columns to guess
                genes_cell = "."
            else:
                genes_cell = row[-4]
        else:
            if gene_col_idx >= len(row):
                genes_cell = "."
            else:
                genes_cell = row[gene_col_idx]

        genes = _split_genes(genes_cell)
        moi_vals: List[str] = []
        dis_vals: List[str] = []

        for g in genes:
            key = g.upper()
            if key in mapping:
                moi_vals.extend(mapping[key].get("moi", []))
                dis_vals.extend(mapping[key].get("disease", []))

        moi_joined = ";".join(_dedup_order_preserving([x for x in moi_vals if x])) if moi_vals else "."
        dis_joined = ";".join(_dedup_order_preserving([x for x in dis_vals if x])) if dis_vals else "."

        out_rows.append(list(row) + [moi_joined, dis_joined])

    return out_rows, out_header


def annotate_moi_paths(
    regions_path: str,
    mapping_path: str,
    output_path: str,
    gene_col_candidates: Tuple[str, ...] = ("Ann_Genes", "Gene", "Genes"),
    out_moi_col: str = "Ann_MOI",
    out_dis_col: str = "Ann_Diseases",
) -> None:
    """
    Library entry point: annotate a regions table with MOI and disease columns
    by looking up gene symbols using a mapping TSV.
    """
    rows, header = read_table(regions_path)
    mapping = load_moi_mapping(mapping_path)
    out_rows, out_header = annotate_moi_rows(
        rows, header, mapping,
        gene_col_candidates=gene_col_candidates,
        out_moi_col=out_moi_col,
        out_dis_col=out_dis_col,
    )
    write_table(output_path, out_rows, out_header)


# ----------------------------
# CLI entry
# ----------------------------

def annotate_moi_cmd(args) -> None:
    """
    CLI entry point compatible with other TRiAs subcommands.
    """
    annotate_moi_paths(
        regions_path=args.regions,
        mapping_path=args.mapping,
        output_path=args.output,
        gene_col_candidates=tuple(args.gene_col_candidates.split(",")) if args.gene_col_candidates else ("Ann_Genes", "Gene", "Genes"),
        out_moi_col=args.out_moi_col,
        out_dis_col=args.out_dis_col,
    )
