from __future__ import annotations
import argparse
import glob
import math
import os
import re
import sys
import shutil
from typing import List, Tuple, Dict, Any, Set
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

DEFAULT_FMT = ".2f"

# ---------- tiny utils ----------
def _is_missing(s):
    if s is None:
        return True
    s = str(s).strip()
    return s == "" or s == "." or s.lower() == "nan"

def _split_multi(s):
    if _is_missing(s):
        return []
    return str(s).strip().split("_")

def _to_float(token):
    if token is None:
        return math.nan
    t = str(token).strip()
    if t in ("", ".", "NaN", "N/A", "NA", "nan"):
        return math.nan
    try:
        return float(t)
    except Exception:
        return math.nan

def _fmt_num(x, fmt=DEFAULT_FMT):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NaN"
    if isinstance(x, int) or (isinstance(x, float) and x.is_integer()):
        return str(int(x))
    return format(x, fmt)

def _join_nums(vals: List[float], fmt=DEFAULT_FMT) -> str:
    return "_".join(_fmt_num(v, fmt) for v in vals) if vals else "-"

def parse_means(s):  # "1.2_3.4"
    return [_to_float(x) for x in _split_multi(s)]

def parse_counts(s):  # "2_0_11"
    return [_to_float(x) for x in _split_multi(s)]

def counts_str(s):
    return "-" if _is_missing(s) else str(s).strip()

def samples_to_set(samples_cell: str) -> Set[str]:
    if _is_missing(samples_cell):
        return set()
    toks = []
    for raw in str(samples_cell).split(","):
        t = raw.strip()
        if t.startswith("+") and t.endswith("+") and len(t) >= 3:
            t = t[1:-1]
        if t:
            toks.append(t)
    return set(toks)

def base_id(sample_token: str) -> str:
    t = str(sample_token).strip()
    if t.endswith(".1") or t.endswith(".2"):
        return t[:-2]
    return t

def hap_suffix(sample_token: str) -> str:
    t = str(sample_token).strip()
    if t.endswith(".1"): return ".1"
    if t.endswith(".2"): return ".2"
    return ""

def is_hap_token(sample_token: str) -> bool:
    t = str(sample_token).strip()
    return t.endswith(".1") or t.endswith(".2")

def chr_key(chr_value: str) -> int:
    v = str(chr_value).replace("chr", "").upper()
    if v == "X": return 23
    if v == "Y": return 24
    if v in ("M", "MT"): return 25
    try:
        return int(v)
    except Exception:
        return 1000

# ---------- MC ratio ----------
def mc_ratio_list(allele_vals: List[float], mean_vals: List[float], min_pop_mean=1.0) -> List[float]:
    L = max(len(allele_vals), len(mean_vals))
    out = []
    for i in range(L):
        av = allele_vals[i] if i < len(allele_vals) else math.nan
        mu = mean_vals[i] if i < len(mean_vals) else math.nan
        if math.isnan(av) or math.isnan(mu) or mu < min_pop_mean or mu == 0.0:
            out.append(math.nan)
        else:
            out.append((av - mu) / mu)
    return out

# ---------- Locus Step1: allowed indices ----------
def _range_max_list(range_str: str) -> List[float]:
    toks = _split_multi(range_str)
    out = []
    for t in toks:
        t = t.strip()
        if "-" in t:
            try:
                out.append(float(t.split("-")[-1]))
            except Exception:
                out.append(math.nan)
        else:
            out.append(_to_float(t))
    return out

def allowed_indices_from_locus(row, case_prefix: str, ctrl_prefix: str, min_pop_mean: float):
    """
    반환:
      allowed_total_idx: set[int]  where mu_total>=min_pop_mean & case_total_max>ctrl_total_max
      allowed_maxc_idx:  set[int]  where mu_maxc >=min_pop_mean & case_maxc_max>ctrl_maxc_max
    """
    # total
    case_tot = _range_max_list(row[f"{case_prefix}_total_repeat_dna_range"])
    ctrl_tot = _range_max_list(row[f"{ctrl_prefix}_total_repeat_dna_range"])
    mu_tot   = parse_means(row["Population_total_repeat_dna_mean"])
    L = max(len(case_tot), len(ctrl_tot), len(mu_tot))
    allowed_total = set()
    for i in range(L):
        ca = case_tot[i] if i < len(case_tot) else math.nan
        cb = ctrl_tot[i] if i < len(ctrl_tot) else math.nan
        mu = mu_tot[i]   if i < len(mu_tot)   else math.nan
        if (not math.isnan(ca)) and (not math.isnan(cb)) and (not math.isnan(mu)) and mu >= min_pop_mean and ca > cb:
            allowed_total.add(i)

    # max_consec
    case_maxc = _range_max_list(row[f"{case_prefix}_max_consec_dna_range"])
    ctrl_maxc = _range_max_list(row[f"{ctrl_prefix}_max_consec_dna_range"])
    mu_maxc   = parse_means(row["Population_max_consec_dna_mean"])
    L2 = max(len(case_maxc), len(ctrl_maxc), len(mu_maxc))
    allowed_maxc = set()
    for i in range(L2):
        ca = case_maxc[i] if i < len(case_maxc) else math.nan
        cb = ctrl_maxc[i] if i < len(ctrl_maxc) else math.nan
        mu = mu_maxc[i]   if i < len(mu_maxc)   else math.nan
        if (not math.isnan(ca)) and (not math.isnan(cb)) and (not math.isnan(mu)) and mu >= min_pop_mean and ca > cb:
            allowed_maxc.add(i)

    return allowed_total, allowed_maxc

# ---------- group helpers ----------
def any_group_present(row: pd.Series, group_cols: List[str]) -> bool:
    for col in group_cols:
        if col in row:
            try:
                if int(str(row[col] or "0")) > 0:
                    return True
            except Exception:
                pass
    return False

def verify_group_columns_exist(df: pd.DataFrame, cols: List[str]) -> List[str]:
    return [c for c in cols if c in df.columns]

# ---------- trio helpers ----------
def collect_counts_for_sample_in_group(group_df: pd.DataFrame, sample_base: str) -> Tuple[bool, str, str]:
    hap_order = [f"{sample_base}.1", f"{sample_base}.2"]
    per_hap_totals = {hap: [] for hap in hap_order}
    per_hap_maxs   = {hap: [] for hap in hap_order}
    any_found = False

    for _, r in group_df.iterrows():
        sset = samples_to_set(r.get("Samples", ""))
        for hap in hap_order:
            if hap in sset:
                any_found = True
                per_hap_totals[hap].append(counts_str(r["total_repeat_dna"]))
                per_hap_maxs[hap].append(counts_str(r["max_consec_dna"]))

    if not any_found:
        return False, "-", "-"

    totals_list, maxs_list = [], []
    for hap in hap_order:
        totals_list.extend(per_hap_totals[hap])
        maxs_list.extend(per_hap_maxs[hap])

    total_str = ",".join(totals_list) if totals_list else "-"
    max_str   = ",".join(maxs_list)   if maxs_list   else "-"
    return True, total_str, max_str

# ---------- chr parsing ----------
CHR_RE = re.compile(r"(?:^|[^A-Za-z0-9])chr(?P<chr>[0-9]+|X|Y|M|MT)(?:[^A-Za-z0-9]|$)", re.IGNORECASE)

def extract_chr_from_path(path: str) -> str:
    m = CHR_RE.search(os.path.basename(path))
    if not m:
        return ""
    raw = m.group("chr").upper()
    if raw == "MT":
        raw = "M"
    return "chr" + raw

def normalize_chr_value(val: str) -> str:
    if val is None:
        return ""
    t = str(val).strip()
    if t == "":
        return ""
    if t.lower().startswith("chr"):
        return "chr" + t[3:].upper().replace("MT", "M")
    return "chr" + t.upper().replace("MT", "M")

# ---------- per-chromosome worker ----------
def process_one_chr(chr_name: str,
                    loci_paths: List[str],
                    allele_paths: List[str],
                    args_dict: Dict[str, Any]) -> Tuple[str, int, int, int]:

    class A: pass
    args = A()
    for k, v in args_dict.items():
        setattr(args, k, v)

    loci_frames = []
    for p in loci_paths:
        try:
            df = pd.read_csv(p, sep="\t", dtype=str)
        except Exception as e:
            sys.stderr.write(f"[WARN] [{chr_name}] Failed to read loci file {p}: {e}\n")
            continue
        if "CHROM" not in df.columns:
            sys.stderr.write(f"[WARN] [{chr_name}] Missing CHROM in loci {p}\n")
            continue
        df = df[df["CHROM"].apply(lambda x: normalize_chr_value(x) == chr_name)]
        if not df.empty:
            loci_frames.append(df)

    if not loci_frames:
        return (chr_name, 0, 0, 0)
    loci = pd.concat(loci_frames, ignore_index=True)

    need_loci = [
        "CHROM","START","END","MOTIF",
        "Population_total_repeat_dna_mean","Population_max_consec_dna_mean",
        f"{args.loci_case_prefix}_total_repeat_dna_range",
        f"{args.loci_case_prefix}_max_consec_dna_range",
        f"{args.loci_ctrl_prefix}_total_repeat_dna_range",
        f"{args.loci_ctrl_prefix}_max_consec_dna_range",
    ]
    miss = [c for c in need_loci if c not in loci.columns]
    if miss:
        sys.stderr.write(f"[ERROR] [{chr_name}] Missing columns in loci files: {miss}\n")
        return (chr_name, 0, 0, 0)

    allele_frames = []
    for p in allele_paths:
        try:
            df = pd.read_csv(p, sep="\t", dtype=str).fillna("")
        except Exception as e:
            sys.stderr.write(f"[WARN] [{chr_name}] Skip allele file {p} (read error: {e})\n")
            continue
        if "CHROM" not in df.columns:
            sys.stderr.write(f"[WARN] [{chr_name}] Missing CHROM in allele {p}\n")
            continue
        df = df[df["CHROM"].apply(lambda x: normalize_chr_value(x) == chr_name)]
        if df.empty:
            continue
        need2 = ["CHROM","START","END","MOTIF","total_repeat_dna","max_consec_dna","Samples"]
        miss2 = [c for c in need2 if c not in df.columns]
        if miss2:
            sys.stderr.write(f"[WARN] [{chr_name}] Skip {p} (missing required allele columns: {miss2})\n")
            continue

        present_case = verify_group_columns_exist(df, args.case_cols)
        present_ctrl = verify_group_columns_exist(df, args.ctrl_cols)
        if not present_case and not present_ctrl:
            sys.stderr.write(f"[WARN] [{chr_name}] Skip {p} (no case/control columns found)\n")
            continue

        allele_frames.append(df)

    if not allele_frames:
        return (chr_name, 0, 0, 0)

    alleles_all = pd.concat(allele_frames, ignore_index=True)
    key_cols = ["CHROM","START","END","MOTIF"]
    alleles_all["_key"] = alleles_all[key_cols].agg("\t".join, axis=1)
    allele_groups = dict(tuple(alleles_all.groupby("_key", sort=False)))

    trio_rows: List[Dict[str,str]] = []
    if args.trio:
        try:
            tdf = pd.read_csv(args.trio, sep="\t", header=None, dtype=str)
        except Exception as e:
            sys.stderr.write(f"[ERROR] [{chr_name}] Cannot read trio: {e}\n")
            tdf = None
        if tdf is not None:
            if tdf.shape[1] < 3:
                sys.stderr.write(f"[ERROR] [{chr_name}] Trio TSV must have 3 columns: child<TAB>father<TAB>mother\n")
            else:
                tdf = tdf.iloc[:, :3]
                tdf.columns = ["child","father","mother"]
                trio_rows = tdf.to_dict(orient="records")

    out_dntre, out_loci, out_all = [], [], []

    for _, lr in loci.iterrows():
        chrom, start, end, motifs_str = str(lr["CHROM"]), str(lr["START"]), str(lr["END"]), str(lr["MOTIF"])
        if normalize_chr_value(chrom) != chr_name:
            continue
        key = "\t".join([chrom, start, end, motifs_str])

        allowed_total, allowed_maxc = allowed_indices_from_locus(
            lr, args.loci_case_prefix, args.loci_ctrl_prefix, args.min_pop_mean
        )
        if not allowed_total and not allowed_maxc:
            continue
        if key not in allele_groups:
            continue

        gdf = allele_groups[key].copy()
        mu_total = parse_means(lr["Population_total_repeat_dna_mean"])
        mu_maxc  = parse_means(lr["Population_max_consec_dna_mean"])

        cache = []
        present_case_cols = [c for c in args.case_cols if c in gdf.columns]
        present_ctrl_cols = [c for c in args.ctrl_cols if c in gdf.columns]

        for _, r in gdf.iterrows():
            if normalize_chr_value(r["CHROM"]) != chr_name:
                continue
            tot_list = parse_counts(r["total_repeat_dna"])
            max_list = parse_counts(r["max_consec_dna"])
            rat_total = mc_ratio_list(tot_list, mu_total, args.min_pop_mean)
            rat_max   = mc_ratio_list(max_list, mu_maxc,  args.min_pop_mean)

            case_present = any_group_present(r, present_case_cols)
            ctrl_present = any_group_present(r, present_ctrl_cols)

            pass_idx_total = sorted([i for i in allowed_total
                                     if i < len(rat_total) and (not math.isnan(rat_total[i])) and (rat_total[i] > args.threshold)])
            pass_idx_maxc  = sorted([i for i in allowed_maxc
                                     if i < len(rat_max) and (not math.isnan(rat_max[i])) and (rat_max[i] > args.threshold)])

            cache.append({
                "row": r,
                "samples_set": samples_to_set(r.get("Samples","")),
                "rat_total": rat_total, "rat_max": rat_max,
                "tot_str": counts_str(r["total_repeat_dna"]),
                "max_str": counts_str(r["max_consec_dna"]),
                "case_present": case_present, "ctrl_present": ctrl_present,
                "pass_idx_total": pass_idx_total,
                "pass_idx_maxc":  pass_idx_maxc,
            })

        any_case_over_total = any(c["case_present"] and len(c["pass_idx_total"])>0 for c in cache)
        any_case_over_max   = any(c["case_present"] and len(c["pass_idx_maxc"]) >0 for c in cache)

        ctrl_all_ok_total = True
        ctrl_all_ok_max   = True
        for c in cache:
            if c["ctrl_present"]:
                if any((i < len(c["rat_total"])) and (not math.isnan(c["rat_total"][i])) and (c["rat_total"][i] > args.threshold)
                       for i in allowed_total):
                    ctrl_all_ok_total = False
                if any((i < len(c["rat_max"])) and (not math.isnan(c["rat_max"][i])) and (c["rat_max"][i] > args.threshold)
                       for i in allowed_maxc):
                    ctrl_all_ok_max = False
                if (not ctrl_all_ok_total) and (not ctrl_all_ok_max):
                    break

        dn_total = "Y" if (any_case_over_total and ctrl_all_ok_total) else "N"
        dn_max   = "Y" if (any_case_over_max   and ctrl_all_ok_max)   else "N"
        if dn_total == "N" and dn_max == "N":
            continue

        trig_total_idx = sorted({i for c in cache if c["case_present"] for i in c["pass_idx_total"]})
        trig_max_idx   = sorted({i for c in cache if c["case_present"] for i in c["pass_idx_maxc"]})
        trig_total_idx_str = ",".join(map(str, trig_total_idx)) if trig_total_idx else "-"
        trig_max_idx_str   = ",".join(map(str, trig_max_idx))   if trig_max_idx   else "-"

        pop_total_mean_str = _join_nums(mu_total, args.fmt)
        pop_max_mean_str   = _join_nums(mu_maxc,  args.fmt)

        out_loci.append({
            "Chrom": chrom, "Start": start, "End": end, "Motifs": motifs_str,
            "dnTRE_total": dn_total, "dnTRE_max": dn_max,
            "Population_total_mean": pop_total_mean_str,
            "Population_max_mean":   pop_max_mean_str,
            "motif_index_total_dna": trig_total_idx_str,
            "motif_index_max_consec": trig_max_idx_str,
        })

        for c in cache:
            out_all.append({
                "Chrom": chrom, "Start": start, "End": end, "Motifs": motifs_str,
                "dnTRE_total": dn_total, "dnTRE_max": dn_max,
                "Allele_origin": ("case" if c["case_present"] else ("control" if c["ctrl_present"] else "none")),
                "TR_allele": c["row"].get("TR_allele",""),
                "Allele_total_MCratio": _join_nums(c["rat_total"], args.fmt),
                "Allele_max_MCratio":   _join_nums(c["rat_max"],   args.fmt),
                "Allele_total_MC": c["tot_str"],
                "Allele_max_MC":   c["max_str"],
                "Population_total_mean": pop_total_mean_str,
                "Population_max_mean":   pop_max_mean_str,
                "Case_present": 1 if c["case_present"] else 0,
                "Control_present": 1 if c["ctrl_present"] else 0,
                "Case_cols": ",".join(args.case_cols),
                "Control_cols": ",".join(args.ctrl_cols),
                "Samples": c["row"]["Samples"],
            })

        if args.trio:
            try:
                tdf = pd.read_csv(args.trio, sep="\t", header=None, dtype=str)
                tdf = tdf.iloc[:, :3]
                tdf.columns = ["child","father","mother"]
                trio_rows = tdf.to_dict(orient="records")
            except Exception:
                trio_rows = []

        if trio_rows:
            for trio in trio_rows:
                child = str(trio["child"]).strip()
                father = str(trio["father"]).strip()
                mother = str(trio["mother"]).strip()

                f_present, f_tot_s, f_max_s = collect_counts_for_sample_in_group(gdf, father)
                m_present, m_tot_s, m_max_s = collect_counts_for_sample_in_group(gdf, mother)
                c_present, c_tot_s, c_max_s = collect_counts_for_sample_in_group(gdf, child)

                father_ids_out = father if f_present else "-"
                mother_ids_out = mother if m_present else "-"

                for c in cache:
                    if not c["case_present"]:
                        continue
                    pass_total = len(c["pass_idx_total"]) > 0
                    pass_max   = len(c["pass_idx_maxc"])  > 0
                    if not (pass_total or pass_max):
                        continue

                    sset = c["samples_set"]
                    child_haps = []
                    if f"{child}.1" in sset: child_haps.append(f"{child}.1")
                    if f"{child}.2" in sset: child_haps.append(f"{child}.2")
                    if not child_haps:
                        continue

                    for ch in child_haps:
                        origin = args.hap1_label if ch.endswith(".1") else args.hap2_label
                        out_dntre.append({
                            "Chrom":  chrom, "Start": start, "End": end, "Motifs": motifs_str,
                            "dnTRE_total": "Y" if pass_total else "N",
                            "dnTRE_max":   "Y" if pass_max   else "N",
                            "dnTRE_total_MC": c["tot_str"],
                            "dnTRE_max_MC":   c["max_str"],
                            "Population_total_mean": pop_total_mean_str,
                            "Population_max_mean":   pop_max_mean_str,
                            "Allele_origin": origin,
                            "motif_index_total_dna": ",".join(map(str, c["pass_idx_total"])) if pass_total else "-",
                            "motif_index_max_consec": ",".join(map(str, c["pass_idx_maxc"]))  if pass_max   else "-",
                            "Father_IDs": father_ids_out,
                            "Mother_IDs": mother_ids_out,
                            "Child_ID":   child,
                            "Father_total_count": f_tot_s,
                            "Father_max_count":   f_max_s,
                            "Mother_total_count": m_tot_s,
                            "Mother_max_count":   m_max_s,
                            "Child_total_count":  c_tot_s if c_present else "-",
                            "Child_max_count":    c_max_s if c_present else "-",
                            "Child_total_MCratio": _join_nums(c["rat_total"], args.fmt),
                            "Child_max_MCratio":   _join_nums(c["rat_max"],   args.fmt),
                        })
        else:
            group_sample_cache: Dict[str, Tuple[bool, str, str]] = {}
            for c in cache:
                if not c["case_present"]:
                    continue
                pass_total = len(c["pass_idx_total"]) > 0
                pass_max   = len(c["pass_idx_maxc"])  > 0
                if not (pass_total or pass_max):
                    continue

                haps = sorted([t for t in c["samples_set"] if is_hap_token(t)])
                if not haps:
                    continue

                for tok in haps:
                    base = base_id(tok)
                    suf  = hap_suffix(tok)
                    origin = args.hap1_label if suf == ".1" else (args.hap2_label if suf == ".2" else "unknown")

                    if base not in group_sample_cache:
                        group_sample_cache[base] = collect_counts_for_sample_in_group(gdf, base)
                    present, all_tot_s, all_max_s = group_sample_cache[base]

                    out_dntre.append({
                        "Chrom":  chrom, "Start": start, "End": end, "Motifs": motifs_str,
                        "dnTRE_total": "Y" if pass_total else "N",
                        "dnTRE_max":   "Y" if pass_max   else "N",
                        "dnTRE_total_MC": c["tot_str"],
                        "dnTRE_max_MC":   c["max_str"],
                        "Population_total_mean": pop_total_mean_str,
                        "Population_max_mean":   pop_max_mean_str,
                        "Sample_origin": origin,
                        "motif_index_total_dna": ",".join(map(str, c["pass_idx_total"])) if pass_total else "-",
                        "motif_index_max_consec": ",".join(map(str, c["pass_idx_maxc"]))  if pass_max   else "-",
                        "Sample_ID": base,
                        "Sample_total_count":  all_tot_s if present else "-",
                        "Sample_max_count":    all_max_s if present else "-",
                        "Sample_total_MCratio": _join_nums(c["rat_total"], args.fmt),
                        "Sample_max_MCratio":   _join_nums(c["rat_max"],   args.fmt),
                    })

    os.makedirs(args.tmpdir, exist_ok=True)

    # dnTRE alleles
    if args.trio:
        cols_dntre = [
            "Chrom","Start","End","Motifs","dnTRE_total","dnTRE_max",
            "dnTRE_total_MC","dnTRE_max_MC",
            "Population_total_mean","Population_max_mean",
            "Allele_origin",
            "motif_index_total_dna","motif_index_max_consec",
            "Father_IDs","Mother_IDs","Child_ID",
            "Father_total_count","Father_max_count",
            "Mother_total_count","Mother_max_count",
            "Child_total_count","Child_max_count",
            "Child_total_MCratio","Child_max_MCratio",
        ]
    else:
        cols_dntre = [
            "Chrom","Start","End","Motifs","dnTRE_total","dnTRE_max",
            "dnTRE_total_MC","dnTRE_max_MC",
            "Population_total_mean","Population_max_mean",
            "Sample_origin",
            "motif_index_total_dna","motif_index_max_consec",
            "Sample_ID",
            "Sample_total_count","Sample_max_count",
            "Sample_total_MCratio","Sample_max_MCratio",
        ]
    pd.DataFrame(out_dntre, columns=cols_dntre).to_csv(
        os.path.join(args.tmpdir, f"dnTRE.{chr_name}.tmp.tsv"), sep="\t", index=False
    )

    cols_loci = [
        "Chrom","Start","End","Motifs","dnTRE_total","dnTRE_max",
        "Population_total_mean","Population_max_mean",
        "motif_index_total_dna","motif_index_max_consec",
    ]
    pd.DataFrame(out_loci, columns=cols_loci).to_csv(
        os.path.join(args.tmpdir, f"dnTRE_loci.{chr_name}.tmp.tsv"), sep="\t", index=False
    )

    cols_all = [
        "Chrom","Start","End","Motifs","dnTRE_total","dnTRE_max",
        "Allele_origin","TR_allele",
        "Allele_total_MC","Allele_max_MC",
        "Allele_total_MCratio","Allele_max_MCratio",
        "Population_total_mean","Population_max_mean",
        "Case_present","Control_present","Case_cols","Control_cols",
        "Samples",
    ]
    pd.DataFrame(out_all, columns=cols_all).to_csv(
        os.path.join(args.tmpdir, f"all_alleles_in_dnTRE_loci.{chr_name}.tmp.tsv"), sep="\t", index=False
    )

    return (chr_name, len(out_dntre), len(out_loci), len(out_all))

# ---------- orchestrator ----------
def default_jobs() -> int:
    try:
        return max(multiprocessing.cpu_count(), 1)
    except Exception:
        return 1

def run_dntre(args: argparse.Namespace) -> None:
    case_cols = [c.strip() for c in str(args.case).split(",") if c.strip()]
    ctrl_cols = [c.strip() for c in str(args.control).split(",") if c.strip()]
    if not case_cols or not ctrl_cols:
        sys.stderr.write("[ERROR] --case and --control must have at least one column each\n")
        sys.exit(1)

    loci_case_prefix = args.loci_case_prefix if args.loci_case_prefix else case_cols[0]
    loci_ctrl_prefix = args.loci_control_prefix if args.loci_control_prefix else ctrl_cols[0]

    loci_paths_all = sorted(glob.glob(args.loci)) if any(ch in args.loci for ch in "*?[]") else [args.loci]
    allele_paths_all = sorted(glob.glob(args.alleles)) if any(ch in args.alleles for ch in "*?[]") else [args.alleles]
    if not loci_paths_all:
        sys.stderr.write(f"[ERROR] No loci files matched: {args.loci}\n")
        sys.exit(1)
    if not allele_paths_all:
        sys.stderr.write(f"[ERROR] No allele files matched: {args.alleles}\n")
        sys.exit(1)

    loci_map: Dict[str, List[str]] = {}
    for p in loci_paths_all:
        c = extract_chr_from_path(p)
        if not c:
            sys.stderr.write(f"[WARN] Could not parse chr from loci file name: {p}\n")
            continue
        loci_map.setdefault(c, []).append(p)

    allele_map: Dict[str, List[str]] = {}
    for p in allele_paths_all:
        c = extract_chr_from_path(p)
        if not c:
            sys.stderr.write(f"[WARN] Could not parse chr from allele file name: {p}\n")
            continue
        allele_map.setdefault(c, []).append(p)

    chr_list = sorted(set(loci_map.keys()) & set(allele_map.keys()),
                      key=lambda c: chr_key(c))
    if not chr_list:
        sys.stderr.write(f"[ERROR] No overlapping chromosomes between loci and alleles by filename\n")
        sys.exit(1)

    args_dict = {
        "threshold": args.threshold,
        "min_pop_mean": args.min_pop_mean,
        "fmt": args.fmt,
        "trio": args.trio,
        "hap1_label": args.hap1_label,
        "hap2_label": args.hap2_label,
        "tmpdir": args.tmpdir,
        "case_cols": case_cols,
        "ctrl_cols": ctrl_cols,
        "loci_case_prefix": loci_case_prefix,
        "loci_ctrl_prefix": loci_ctrl_prefix,
    }

    os.makedirs(args.tmpdir, exist_ok=True)

    futures = []
    with ProcessPoolExecutor(max_workers=max(1, args.jobs)) as ex:
        for c in chr_list:
            futures.append(ex.submit(
                process_one_chr, c, loci_map.get(c, []), allele_map.get(c, []), args_dict
            ))
        for fut in as_completed(futures):
            chr_name, n_d, n_l, n_a = fut.result()
            sys.stderr.write(f"[INFO] {chr_name}: dnTRE={n_d}, loci={n_l}, all={n_a}\n")

    # Merge & sort
    def concat_sort_write(pattern: str, out_path: str):
        paths = sorted(glob.glob(os.path.join(args.tmpdir, pattern)))
        frames = []
        for p in paths:
            try:
                frames.append(pd.read_csv(p, sep="\t", dtype=str))
            except Exception as e:
                sys.stderr.write(f"[WARN] Failed to read tmp {p}: {e}\n")
        if not frames:
            pd.DataFrame(columns=[]).to_csv(out_path, sep="\t", index=False)
            return
        df = pd.concat(frames, ignore_index=True)

        cols = [c for c in ["Chrom","Start","End"] if c in df.columns]
        if cols:
            def keyfunc(col: pd.Series):
                if col.name == "Chrom":
                    return col.map(chr_key)
                if col.name in ("Start","End"):
                    return col.astype(int, errors="ignore")
                return col
            df = df.sort_values(by=cols, key=keyfunc)

        df.to_csv(out_path, sep="\t", index=False)

    concat_sort_write("dnTRE.*.tmp.tsv", args.out_dntre)
    concat_sort_write("dnTRE_loci.*.tmp.tsv", args.out_loci)
    concat_sort_write("all_alleles_in_dnTRE_loci.*.tmp.tsv", args.out_alleles)

    if not args.keep_temp:
        try:
            shutil.rmtree(args.tmpdir)
        except Exception as e:
            sys.stderr.write(f"[WARN] Failed to remove tmpdir {args.tmpdir}: {e}\n")
