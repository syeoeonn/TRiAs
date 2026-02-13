# src/trias/allele/report_dntre.py
from __future__ import annotations
import sys
import math
import glob
from typing import Set, Tuple, Dict, Any, List
import pandas as pd

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

def parse_means(s):
    toks = _split_multi(s)
    return [_to_float(x) for x in toks]

def parse_counts(s):
    toks = _split_multi(s)
    return [_to_float(x) for x in toks]

def counts_str(s):
    if _is_missing(s):
        return "-"
    return str(s).strip()

def parse_ratio_list(s):
    toks = _split_multi(s)
    return [_to_float(x) for x in toks]

def chr_key(chr_value: str) -> int:
    """Natural chromosome order: 1..22, X=23, Y=24, M/MT=25."""
    v = str(chr_value).replace("chr", "").upper()
    if v == "X": return 23
    if v == "Y": return 24
    if v in ("M", "MT"): return 25
    return int(v)

def mc_ratio_list(allele_vals, mean_vals, min_pop_mean=1.0):
    """
    Per-motif MC ratio: (allele_value - population_mean) / population_mean
    If pop_mean < min_pop_mean or NaN -> NaN
    """
    L = max(len(allele_vals), len(mean_vals))
    out = []
    for i in range(L):
        av = allele_vals[i] if i < len(allele_vals) else math.nan
        mu = mean_vals[i] if i < len(mean_vals) else math.nan
        if math.isnan(av) or math.isnan(mu) or mu < min_pop_mean:
            out.append(math.nan)
        else:
            out.append((av - mu) / mu)
    return out

def indices_ge_threshold(vals, indices, thr):
    out = []
    for i in indices:
        if 0 <= i < len(vals):
            v = vals[i]
            if not math.isnan(v) and v >= thr:
                out.append(i)
    return out

def clean_sample_token(tok: str) -> str:
    if tok is None:
        return ""
    t = tok.strip()
    if t.startswith("+") and t.endswith("+") and len(t) >= 3:
        t = t[1:-1]
    return t

def samples_to_set(samples_cell: str) -> Set[str]:
    if _is_missing(samples_cell):
        return set()
    toks = [clean_sample_token(x) for x in str(samples_cell).split(",")]
    return set([t for t in toks if t])

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

def collect_counts_for_sample_in_group(group_df: pd.DataFrame, sample_base: str) -> Tuple[bool, str, str]:
    """
    From all allele rows in this locus, gather total/max count strings for sample_base.
    Order: all .1 entries (row order), then all .2 entries (row order).
    Returns: present(bool), total_str, max_str
    """
    hap_order = [f"{sample_base}.1", f"{sample_base}.2"]
    per_hap_totals = {hap: [] for hap in hap_order}
    per_hap_maxs   = {hap: [] for hap in hap_order}
    any_found = False

    for _, r in group_df.iterrows():
        allele_samples = samples_to_set(r.get("Samples", ""))
        for hap in hap_order:
            if hap in allele_samples:
                any_found = True
                per_hap_totals[hap].append(counts_str(r["total_repeat_dna"]))
                per_hap_maxs[hap].append(counts_str(r["max_consec_dna"]))

    if not any_found:
        return False, "-", "-"

    totals_list, maxs_list = [], []
    for hap in hap_order:  # .1 first then .2
        totals_list.extend(per_hap_totals[hap])
        maxs_list.extend(per_hap_maxs[hap])

    total_str = ",".join(totals_list) if totals_list else "-"
    max_str   = ",".join(maxs_list)   if maxs_list   else "-"
    return True, total_str, max_str

def report_dntre_cmd(args):
    """
    CLI entrypoint for `trias report-dntre`.
    """
    # Trio (optional)
    has_trio = bool(args.trio)
    trio_rows: List[Dict[str, str]] = []
    if has_trio:
        trio_df = pd.read_csv(args.trio, sep="\t", header=None, dtype=str)
        if trio_df.shape[1] >= 3:
            trio_df = trio_df.iloc[:, :3]
            trio_df.columns = ["child", "father", "mother"]
            trio_rows = trio_df.to_dict(orient="records")
        else:
            sys.stderr.write("[ERROR] Trio TSV must have 3 columns: child<TAB>father<TAB>mother\n")
            return 1

    # file1 -> locus map
    df1 = pd.read_csv(args.mc, sep="\t", dtype=str)
    required_file1 = [
        "CHROM","START","END","MOTIF",
        "High_MCratio_total","High_MCratio_maxconsec",
        "Population_total_repeat_dna_mean","Population_max_consec_dna_mean",
        "Our_total_repeat_dna_MCratio","Our_max_consec_dna_MCratio",
        "Unaffected_total_repeat_dna_MCratio","Unaffected_max_consec_dna_MCratio",
    ]
    missing = [c for c in required_file1 if c not in df1.columns]
    if missing:
        sys.stderr.write(f"[ERROR] Missing columns in mc TSV: {missing}\n")
        return 1

    locus_map: Dict[Any, Any] = {}
    for _, r in df1[required_file1].dropna(subset=["CHROM","START","END","MOTIF"]).iterrows():
        key = (str(r["CHROM"]), str(r["START"]), str(r["END"]), str(r["MOTIF"]))
        check_total = str(r["High_MCratio_total"]).strip().upper() == "Y"
        check_maxc  = str(r["High_MCratio_maxconsec"]).strip().upper() == "Y"

        mu_total = parse_means(r["Population_total_repeat_dna_mean"])
        mu_maxc  = parse_means(r["Population_max_consec_dna_mean"])

        our_total = parse_ratio_list(r["Our_total_repeat_dna_MCratio"])
        una_total = parse_ratio_list(r["Unaffected_total_repeat_dna_MCratio"])
        our_maxc  = parse_ratio_list(r["Our_max_consec_dna_MCratio"])
        una_maxc  = parse_ratio_list(r["Unaffected_max_consec_dna_MCratio"])

        allowed_idx_total = set()
        if check_total:
            L = max(len(our_total), len(una_total))
            for i in range(L):
                cr = our_total[i] if i < len(our_total) else math.nan
                ur = una_total[i] if i < len(una_total) else math.nan
                if (not math.isnan(cr)) and (not math.isnan(ur)) and (cr >= args.threshold) and (ur <= args.threshold):
                    allowed_idx_total.add(i)

        allowed_idx_maxc = set()
        if check_maxc:
            L = max(len(our_maxc), len(una_maxc))
            for i in range(L):
                cr = our_maxc[i] if i < len(our_maxc) else math.nan
                ur = una_maxc[i] if i < len(una_maxc) else math.nan
                if (not math.isnan(cr)) and (not math.isnan(ur)) and (cr >= args.threshold) and (ur <= args.threshold):
                    allowed_idx_maxc.add(i)

        locus_map[key] = {
            "check_total": check_total,
            "check_maxc":  check_maxc,
            "mu_total":    mu_total,
            "mu_maxc":     mu_maxc,
            "allowed_idx_total": allowed_idx_total,
            "allowed_idx_maxc":  allowed_idx_maxc,
        }

    keys_df = pd.DataFrame([(k[0],k[1],k[2],k[3]) for k in locus_map.keys()],
                           columns=["CHROM","START","END","MOTIF"])

    # input2 files
    file2_list = sorted(glob.glob(args.population))
    if not file2_list:
        sys.stderr.write(f"[ERROR] No files matched pattern: {args.population}\n")
        return 1

    out_is_stdout = (args.output == "-" or args.output.lower() == "stdout")
    all_batches: List[pd.DataFrame] = []

    # column sets
    shared_leading = [
        "Chrom","Start","End","Motifs",
        "dnTRE_total","dnTRE_max",
        "dnTRE_total_MC","dnTRE_max_MC",
        "Population_total_mean","Population_max_mean",
    ]
    motif_cols = ["motif_index_total_dna","motif_index_max_consec"]

    if has_trio:
        out_cols = shared_leading + [
            "Allele_origin",
            *motif_cols,
            "Father_IDs","Mother_IDs","Child_ID",
            "Father_total_count","Father_max_count",
            "Mother_total_count","Mother_max_count",
            "Child_total_count","Child_max_count",
            "Child_total_MCratio","Child_max_MCratio",
        ]
    else:
        out_cols = shared_leading + [
            "Sample_origin",
            *motif_cols,
            "Sample_ID",
            "Sample_total_count","Sample_max_count",
            "Sample_total_MCratio","Sample_max_MCratio",
        ]

    for f in file2_list:
        df2 = pd.read_csv(f, sep="\t", dtype=str)
        req2 = ["CHROM","START","END","MOTIF","total_repeat_dna","max_consec_dna","Samples"]
        miss2 = [c for c in req2 if c not in df2.columns]
        if miss2:
            sys.stderr.write(f"[WARN] Skip {f} (missing columns: {miss2})\n")
            continue

        if keys_df.empty:
            break

        merged = df2.merge(keys_df, on=["CHROM","START","END","MOTIF"], how="inner")
        if merged.empty:
            continue

        # group by locus
        out_rows: List[Dict[str, Any]] = []
        for (g_chrom, g_start, g_end, g_motif), gdf in merged.groupby(["CHROM","START","END","MOTIF"], sort=False):
            key = (g_chrom, g_start, g_end, g_motif)
            info = locus_map.get(key)
            if not info:
                continue

            mu_total_list = info["mu_total"]
            mu_maxc_list  = info["mu_maxc"]
            allowed_total = sorted(info["allowed_idx_total"])
            allowed_maxc  = sorted(info["allowed_idx_maxc"])

            # per-allele cache
            allele_cache = []
            for _, r in gdf.iterrows():
                allele_samples = samples_to_set(r.get("Samples", ""))

                allele_total_counts_list = parse_counts(r["total_repeat_dna"])
                allele_maxc_counts_list  = parse_counts(r["max_consec_dna"])
                ratios_total = mc_ratio_list(allele_total_counts_list, mu_total_list, args.min_pop_mean)
                ratios_maxc  = mc_ratio_list(allele_maxc_counts_list,  mu_maxc_list,  args.min_pop_mean)

                pass_idx_total = indices_ge_threshold(ratios_total, allowed_total, args.threshold) if info["check_total"] else []
                pass_idx_maxc  = indices_ge_threshold(ratios_maxc,  allowed_maxc,  args.threshold) if info["check_maxc"]  else []
                pass_total = len(pass_idx_total) > 0
                pass_max   = len(pass_idx_maxc)  > 0

                allele_cache.append({
                    "row": r,
                    "samples_set": allele_samples,
                    "ratios_total": ratios_total,
                    "ratios_maxc":  ratios_maxc,
                    "pass_idx_total": pass_idx_total,
                    "pass_idx_maxc":  pass_idx_maxc,
                    "pass_total": pass_total,
                    "pass_max":   pass_max,
                    "total_count_str": counts_str(r["total_repeat_dna"]),
                    "max_count_str":   counts_str(r["max_consec_dna"]),
                })

            pop_total_mean_str = "_".join(_fmt_num(x, args.fmt) for x in mu_total_list) if mu_total_list else "-"
            pop_max_mean_str   = "_".join(_fmt_num(x, args.fmt) for x in mu_maxc_list)  if mu_maxc_list  else "-"

            if has_trio:
                # Trio mode: per child-hap rows on dnTRE allele rows
                for trio in trio_rows:
                    child_base  = str(trio["child"]).strip()
                    father_base = str(trio["father"]).strip()
                    mother_base = str(trio["mother"]).strip()

                    # locus-wide counts
                    f_present, f_tot_s, f_max_s = collect_counts_for_sample_in_group(gdf, father_base)
                    m_present, m_tot_s, m_max_s = collect_counts_for_sample_in_group(gdf, mother_base)
                    c_present, c_tot_s, c_max_s = collect_counts_for_sample_in_group(gdf, child_base)

                    father_ids_out = father_base if f_present else "-"
                    mother_ids_out = mother_base if m_present else "-"

                    for ac in allele_cache:
                        if not (ac["pass_total"] or ac["pass_max"]):
                            continue
                        sset = ac["samples_set"]
                        # child hap tokens inside this allele row
                        ch_haps = []
                        h1 = f"{child_base}.1"
                        h2 = f"{child_base}.2"
                        if h1 in sset: ch_haps.append(h1)
                        if h2 in sset: ch_haps.append(h2)
                        if not ch_haps:
                            continue

                        for ch in ch_haps:
                            allele_origin = args.hap1_label if ch.endswith(".1") else args.hap2_label
                            out_rows.append({
                                "Chrom":  g_chrom,
                                "Start":  g_start,
                                "End":    g_end,
                                "Motifs": g_motif,
                                "dnTRE_total": "Y" if ac["pass_total"] else "N",
                                "dnTRE_max":   "Y" if ac["pass_max"]   else "N",
                                "dnTRE_total_MC": ac["total_count_str"],
                                "dnTRE_max_MC":   ac["max_count_str"],
                                "Population_total_mean": pop_total_mean_str,
                                "Population_max_mean":   pop_max_mean_str,
                                "Allele_origin": allele_origin,
                                "motif_index_total_dna": ",".join(map(str, ac["pass_idx_total"])) if ac["pass_idx_total"] else "-",
                                "motif_index_max_consec": ",".join(map(str, ac["pass_idx_maxc"])) if ac["pass_idx_maxc"] else "-",
                                "Father_IDs": father_ids_out,
                                "Mother_IDs": mother_ids_out,
                                "Child_ID":   child_base,
                                "Father_total_count": f_tot_s,
                                "Father_max_count":   f_max_s,
                                "Mother_total_count": m_tot_s,
                                "Mother_max_count":   m_max_s,
                                "Child_total_count":  c_tot_s if c_present else "-",
                                "Child_max_count":    c_max_s if c_present else "-",
                                "Child_total_MCratio": "_".join(_fmt_num(x, args.fmt) for x in ac["ratios_total"]) if ac["ratios_total"] else "-",
                                "Child_max_MCratio":   "_".join(_fmt_num(x, args.fmt) for x in ac["ratios_maxc"])  if ac["ratios_maxc"]  else "-",
                            })
            else:
                # Auto-sample mode: per sample-hap rows for dnTRE allele rows
                for ac in allele_cache:
                    if not (ac["pass_total"] or ac["pass_max"]):
                        continue

                    sset = ac["samples_set"]
                    hap_tokens = sorted([t for t in sset if is_hap_token(t)])
                    if not hap_tokens:
                        continue

                    group_sample_cache: Dict[str, Tuple[bool, str, str]] = {}

                    for tok in hap_tokens:
                        base = base_id(tok)
                        suffix = hap_suffix(tok)
                        if suffix == ".1":
                            origin = args.hap1_label
                        elif suffix == ".2":
                            origin = args.hap2_label
                        else:
                            continue

                        if base not in group_sample_cache:
                            group_sample_cache[base] = collect_counts_for_sample_in_group(gdf, base)
                        present, all_tot_s, all_max_s = group_sample_cache[base]

                        out_rows.append({
                            "Chrom":  g_chrom,
                            "Start":  g_start,
                            "End":    g_end,
                            "Motifs": g_motif,
                            "dnTRE_total": "Y" if ac["pass_total"] else "N",
                            "dnTRE_max":   "Y" if ac["pass_max"]   else "N",
                            "dnTRE_total_MC": ac["total_count_str"],
                            "dnTRE_max_MC":   ac["max_count_str"],
                            "Population_total_mean": pop_total_mean_str,
                            "Population_max_mean":   pop_max_mean_str,
                            "Sample_origin": origin,
                            "motif_index_total_dna": ",".join(map(str, ac["pass_idx_total"])) if ac["pass_idx_total"] else "-",
                            "motif_index_max_consec": ",".join(map(str, ac["pass_idx_maxc"])) if ac["pass_idx_maxc"] else "-",
                            "Sample_ID": base,
                            "Sample_total_count":  all_tot_s if present else "-",
                            "Sample_max_count":    all_max_s if present else "-",
                            "Sample_total_MCratio": "_".join(_fmt_num(x, args.fmt) for x in ac["ratios_total"]) if ac["ratios_total"] else "-",
                            "Sample_max_MCratio":   "_".join(_fmt_num(x, args.fmt) for x in ac["ratios_maxc"])  if ac["ratios_maxc"]  else "-",
                        })

        if not out_rows:
            continue

        out_df = pd.DataFrame(out_rows, columns=out_cols)
        all_batches.append(out_df)

    # write out
    if all_batches:
        final_df = pd.concat(all_batches, ignore_index=True)
        sort_cols = [c for c in ["Chrom","Start","End"] if c in final_df.columns]
        final_df = final_df.sort_values(
            by=sort_cols,
            key=lambda col: (
                col.map(chr_key) if col.name == "Chrom"
                else col.astype(int)
            )
        )
        if out_is_stdout:
            final_df.to_csv(sys.stdout, sep="\t", index=False)
        else:
            final_df.to_csv(args.output, sep="\t", index=False)
    else:
        if out_is_stdout:
            pd.DataFrame(columns=out_cols).to_csv(sys.stdout, sep="\t", index=False)
        else:
            pd.DataFrame(columns=out_cols).to_csv(args.output, sep="\t", index=False)

    return 0
