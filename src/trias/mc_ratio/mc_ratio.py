# src/trias/mc_ratio.py
import argparse
import math
import sys
import pandas as pd

DEFAULT_FMT = ".2f"

# ---------- Utility functions ----------
def _is_missing(s):
    if s is None:
        return True
    s = str(s).strip()
    return s == "" or s == "." or s.lower() == "nan"

def _split_multi(s):
    if _is_missing(s):
        return []
    return str(s).strip().split("_")

def _range_to_max(token):
    if token is None:
        return math.nan
    t = str(token).strip()
    if t == "" or t == ".":
        return math.nan
    if "-" in t:
        try:
            return float(t.split("-")[-1])
        except Exception:
            return math.nan
    try:
        return float(t)
    except Exception:
        return math.nan

def _to_float(token):
    if token is None:
        return math.nan
    t = str(token).strip()
    if t == "" or t == ".":
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

def compute_mc_ratio_list(range_str, pop_mean_str, min_pop_mean=1.0):
    range_tokens = _split_multi(range_str)
    pop_tokens = _split_multi(pop_mean_str)
    L = max(len(range_tokens), len(pop_tokens))
    ratios = []
    for i in range(L):
        range_max = _range_to_max(range_tokens[i] if i < len(range_tokens) else None)
        pop_mean = _to_float(pop_tokens[i] if i < len(pop_tokens) else None)
        if (pop_mean is None or math.isnan(pop_mean) or pop_mean < min_pop_mean
            or math.isnan(range_max) or pop_mean == 0.0):
            ratios.append(math.nan)
        else:
            ratios.append((range_max - pop_mean) / pop_mean)
    return ratios

def compute_max_list(range_str):
    toks = _split_multi(range_str)
    return [_range_to_max(t) for t in toks]

def compute_diff_list(case_range_str, control_range_str):
    case_max_list = compute_max_list(case_range_str)
    ctrl_max_list = compute_max_list(control_range_str)
    L = max(len(case_max_list), len(ctrl_max_list))
    diffs = []
    for i in range(L):
        cm = case_max_list[i] if i < len(case_max_list) else math.nan
        um = ctrl_max_list[i] if i < len(ctrl_max_list) else math.nan
        if math.isnan(cm) or math.isnan(um):
            diffs.append(math.nan)
        else:
            diffs.append(cm - um)
    return diffs

def join_ratios(ratios, fmt=DEFAULT_FMT):
    return "_".join(_fmt_num(x, fmt) for x in ratios)

def any_gt_threshold(ratios, thr=2.0):
    return any((not math.isnan(x)) and x > thr for x in ratios)

def any_dual_filter_condition(case_ratios, ctrl_ratios,
                              case_max_list=None, ctrl_max_list=None,
                              thr=2.0, use_max_condition=False):
    L = max(len(case_ratios), len(ctrl_ratios))
    if use_max_condition and (case_max_list is not None and ctrl_max_list is not None):
        L = max(L, len(case_max_list), len(ctrl_max_list))

    for i in range(L):
        cr = case_ratios[i] if i < len(case_ratios) else math.nan
        ur = ctrl_ratios[i] if i < len(ctrl_ratios) else math.nan
        if not math.isnan(cr) and not math.isnan(ur):
            dual_condition = (cr > thr) and (ur < thr)
            if use_max_condition and dual_condition:
                cm = case_max_list[i] if i < len(case_max_list) else math.nan
                um = ctrl_max_list[i] if i < len(ctrl_max_list) else math.nan
                if not math.isnan(cm) and not math.isnan(um):
                    if cm > um:
                        return True
                elif dual_condition:  # max 값 결측이면 ratio 조건만
                    return True
            elif dual_condition:
                return True
    return False

def any_ratio_and_max_condition(ratios, case_max_list, ctrl_max_list, thr=2.0):
    L = max(len(ratios), len(case_max_list), len(ctrl_max_list))
    for i in range(L):
        r = ratios[i] if i < len(ratios) else math.nan
        cm = case_max_list[i] if i < len(case_max_list) else math.nan
        um = ctrl_max_list[i] if i < len(ctrl_max_list) else math.nan
        if not math.isnan(r) and not math.isnan(cm) and not math.isnan(um):
            if r > thr and cm > um:
                return True
    return False

# ---------- Main function ----------
def mc_ratio_cmd(args):
    df = pd.read_csv(args.input, sep="\t", dtype=str)

    required = [
        f"{args.case}_total_repeat_dna_range",
        f"{args.control}_total_repeat_dna_range",
        f"{args.case}_max_consec_dna_range",
        f"{args.control}_max_consec_dna_range",
        "Population_total_repeat_dna_mean",
        "Population_max_consec_dna_mean",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    case_total_ratio_str = []
    case_maxc_ratio_str = []
    ctrl_total_ratio_str = []
    ctrl_maxc_ratio_str = []
    flag_total = []
    flag_maxc = []
    diff_total_str, diff_maxc_str = [], []

    for _, row in df.iterrows():
        # Case ratios
        r_total = compute_mc_ratio_list(row[f"{args.case}_total_repeat_dna_range"],
                                        row["Population_total_repeat_dna_mean"],
                                        args.min_pop_mean)
        r_maxc = compute_mc_ratio_list(row[f"{args.case}_max_consec_dna_range"],
                                       row["Population_max_consec_dna_mean"],
                                       args.min_pop_mean)
        case_total_ratio_str.append(join_ratios(r_total, args.fmt))
        case_maxc_ratio_str.append(join_ratios(r_maxc, args.fmt))

        # Control + flags
        if args.dual_filter:
            ur_total = compute_mc_ratio_list(row[f"{args.control}_total_repeat_dna_range"],
                                             row["Population_total_repeat_dna_mean"],
                                             args.min_pop_mean)
            ur_maxc = compute_mc_ratio_list(row[f"{args.control}_max_consec_dna_range"],
                                            row["Population_max_consec_dna_mean"],
                                            args.min_pop_mean)
            ctrl_total_ratio_str.append(join_ratios(ur_total, args.fmt))
            ctrl_maxc_ratio_str.append(join_ratios(ur_maxc, args.fmt))

            if args.max:
                case_total_max_list = compute_max_list(row[f"{args.case}_total_repeat_dna_range"])
                ctrl_total_max_list = compute_max_list(row[f"{args.control}_total_repeat_dna_range"])
                case_maxc_max_list = compute_max_list(row[f"{args.case}_max_consec_dna_range"])
                ctrl_maxc_max_list = compute_max_list(row[f"{args.control}_max_consec_dna_range"])

                flag_total.append("Y" if any_dual_filter_condition(
                    r_total, ur_total,
                    case_total_max_list, ctrl_total_max_list,
                    args.threshold, use_max_condition=True) else "N")
                flag_maxc.append("Y" if any_dual_filter_condition(
                    r_maxc, ur_maxc,
                    case_maxc_max_list, ctrl_maxc_max_list,
                    args.threshold, use_max_condition=True) else "N")
            else:
                flag_total.append("Y" if any_dual_filter_condition(r_total, ur_total, thr=args.threshold) else "N")
                flag_maxc.append("Y" if any_dual_filter_condition(r_maxc, ur_maxc, thr=args.threshold) else "N")

        else:  # not dual-filter
            ctrl_total_ratio_str.append("N/A")
            ctrl_maxc_ratio_str.append("N/A")

            if args.max:
                case_total_max_list = compute_max_list(row[f"{args.case}_total_repeat_dna_range"])
                ctrl_total_max_list = compute_max_list(row[f"{args.control}_total_repeat_dna_range"])
                case_maxc_max_list = compute_max_list(row[f"{args.case}_max_consec_dna_range"])
                ctrl_maxc_max_list = compute_max_list(row[f"{args.control}_max_consec_dna_range"])

                flag_total.append("Y" if any_ratio_and_max_condition(
                    r_total, case_total_max_list, ctrl_total_max_list, args.threshold) else "N")
                flag_maxc.append("Y" if any_ratio_and_max_condition(
                    r_maxc, case_maxc_max_list, ctrl_maxc_max_list, args.threshold) else "N")
            else:
                flag_total.append("Y" if any_gt_threshold(r_total, args.threshold) else "N")
                flag_maxc.append("Y" if any_gt_threshold(r_maxc, args.threshold) else "N")

        # Diff
        if args.diff:
            total_diff = compute_diff_list(row[f"{args.case}_total_repeat_dna_range"],
                                           row[f"{args.control}_total_repeat_dna_range"])
            maxc_diff = compute_diff_list(row[f"{args.case}_max_consec_dna_range"],
                                          row[f"{args.control}_max_consec_dna_range"])
            diff_total_str.append(join_ratios(total_diff, args.fmt))
            diff_maxc_str.append(join_ratios(maxc_diff, args.fmt))

    # Add columns
    df.insert(0, "High_MCratio_total", flag_total)
    df.insert(1, "High_MCratio_maxconsec", flag_maxc)
    df[f"{args.case}_total_repeat_dna_MCratio"] = case_total_ratio_str
    df[f"{args.case}_max_consec_dna_MCratio"] = case_maxc_ratio_str
    df[f"{args.control}_total_repeat_dna_MCratio"] = ctrl_total_ratio_str
    df[f"{args.control}_max_consec_dna_MCratio"] = ctrl_maxc_ratio_str
    if args.diff:
        df[f"{args.case}_vs_{args.control}_total_diff"] = diff_total_str
        df[f"{args.case}_vs_{args.control}_maxconsec_diff"] = diff_maxc_str

    # Write output
    if args.output == "-" or args.output.lower() == "stdout":
        df.to_csv(sys.stdout, sep="\t", index=False)
    else:
        df.to_csv(args.output, sep="\t", index=False)
