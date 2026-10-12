"""v7 scoring with patient-level replication + per-non-target-cell-type contamination.

For each target cell_type, output per peak:
  Replication (patient-level primary, cohort-level secondary):
    n_target_patients      — count of target patients where this peak was detected
    n_target_cohorts       — count of distinct cohorts among those patients
    cohort_balance         — fraction of target signal from top-contributing cohort
                              (close to 1/N_cohorts = healthy; 1 = cohort-dominated)
  Specificity (traditional):
    target_mean_cpm        — mean CPM across target patient columns
    bg_mean_cpm            — mean CPM across all non-target patient columns
    daigle_z               — (target_mean - bg_mean) / bg_std across non-target cols
  AAV off-target safety (NEW):
    offtarget_max_cpm      — max per-cell-type mean CPM among non-target cell types
    offtarget_top_cell_type— which non-target cell type is the worst
    offtarget_pooled_mean  — mean across non-target cell-type means (uniform avg)
    n_offtarget_ct_above_1 — number of non-target cell types with mean CPM > 1

Pre-rank filter: target_mean_cpm >= CPM_FLOOR to clip mathematical-outlier Z.
"""
import sys
import time
import re
from pathlib import Path

import numpy as np
import pandas as pd

# v1.1 CNV-aware annotation + arm-adjusted Z
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cnv_patch import (annotate_cnv_class_and_arm, arm_adjusted_z,
                       patient_fraction_elevated)

ROOT = Path("/data/projects/atacseq")
PILOT = ROOT / "pilot_v5"
REQ = PILOT / "requant"
RPT = PILOT / "reports_v7"
RPT.mkdir(exist_ok=True, parents=True)

LOG = PILOT / "score_v7.log"
CPM_FLOOR = 0.5


def log(m):
    s = f"[score-v7] {m}"
    print(s, flush=True)
    with open(LOG, "a") as f:
        f.write(s + "\n")


COHORTS = ["greenwald_2024", "sundaram_2024", "sojka_2025", "zhao_2024",
           "mathewson_2021", "guo_2020", "guilhamon_2021", "nagaraja_2021"]


def parse_col(col):
    """<cohort>__<patient>__<label> → (cohort, patient, label)."""
    m = re.match(r"^([^_]+(?:_\d{4}))__(.+)__([^_]+(?:_[^_]+)*)$", col)
    if m:
        return m.group(1), m.group(2), m.group(3)
    # fallback: split on __
    parts = col.split("__", 2)
    if len(parts) == 3:
        return parts[0], parts[1], parts[2]
    return None, None, None


def load_v7():
    peaks = pd.read_parquet(REQ / "v5_peak_index.parquet")
    log(f"peaks: {len(peaks):,}")
    pbs = []
    for c in COHORTS:
        p = REQ / f"pseudobulk_v7_{c}.parquet"
        if p.exists():
            df = pd.read_parquet(p)
            log(f"  {c}: {df.shape}")
            pbs.append(df)
    big_raw = pd.concat(pbs, axis=1)
    log(f"combined raw: {big_raw.shape}")
    # Parse column meta
    col_meta = {}
    for col in big_raw.columns:
        cohort, patient, label = parse_col(col)
        col_meta[col] = {"cohort": cohort, "patient": patient, "label": label}
    return peaks, big_raw, col_meta


def cpm_normalize(df):
    out = df.copy().astype(float)
    for col in out.columns:
        s = out[col].sum()
        if s > 0:
            out[col] = out[col] * 1e6 / s
    return out


def score_one_target(target_ct, big_raw, big_cpm, col_meta, peaks):
    target_cols = [c for c, m in col_meta.items() if m["label"] == target_ct]
    bg_cols = [c for c in big_cpm.columns if c not in target_cols]
    if len(target_cols) == 0 or len(bg_cols) < 2:
        log(f"  skip {target_ct}: targets={len(target_cols)} bg={len(bg_cols)}")
        return None

    # ---- target metrics
    target_mean = big_cpm[target_cols].mean(axis=1)
    # patient-level replication: a column has signal if raw count > 0
    target_detected = (big_raw[target_cols] > 0)
    n_target_patients = target_detected.sum(axis=1)

    # cohort of each target column
    target_cohorts = pd.Series({c: col_meta[c]["cohort"] for c in target_cols})

    # n_target_cohorts: number of distinct cohorts where ≥1 target patient has signal
    # cohort_balance: fraction of total target CPM from top-contributing cohort
    cohort_cpm_sums = {}
    for cohort in target_cohorts.unique():
        cohort_target_cols = [c for c in target_cols if col_meta[c]["cohort"] == cohort]
        cohort_cpm_sums[cohort] = big_cpm[cohort_target_cols].sum(axis=1)
    cohort_cpm_df = pd.DataFrame(cohort_cpm_sums)
    row_sum = cohort_cpm_df.sum(axis=1).replace(0, np.nan)
    top_cohort_frac = cohort_cpm_df.max(axis=1) / row_sum
    n_target_cohorts = (cohort_cpm_df > 0).sum(axis=1)

    # ---- bg (specificity) metrics
    bg_mean = big_cpm[bg_cols].mean(axis=1)
    bg_std = big_cpm[bg_cols].std(axis=1).replace(0, np.nan)
    daigle_z = ((target_mean - bg_mean) / bg_std).fillna(0)

    # ---- per-non-target-cell-type contamination
    other_labels = {m["label"] for c, m in col_meta.items() if m["label"] != target_ct and m["label"] not in ("", "nan")}
    per_ct_mean = {}
    for lab in other_labels:
        lab_cols = [c for c, m in col_meta.items() if m["label"] == lab]
        if lab_cols:
            per_ct_mean[lab] = big_cpm[lab_cols].mean(axis=1)
    per_ct_df = pd.DataFrame(per_ct_mean)
    offtarget_max_cpm = per_ct_df.max(axis=1)
    offtarget_top_cell_type = per_ct_df.idxmax(axis=1)
    offtarget_pooled_mean = per_ct_df.mean(axis=1)
    n_offtarget_ct_above_1 = (per_ct_df > 1.0).sum(axis=1)

    # ---- assemble + filter
    df = peaks[["chrom", "start", "end", "peak_id"]].copy()
    df["target_mean_cpm"] = target_mean.values
    df["bg_mean_cpm"] = bg_mean.values
    df["daigle_z_raw"] = daigle_z.values
    df["n_target_patients"] = n_target_patients.values
    df["n_target_cohorts"] = n_target_cohorts.values
    df["cohort_balance"] = top_cohort_frac.fillna(0).values
    df["offtarget_max_cpm"] = offtarget_max_cpm.values
    df["offtarget_top_cell_type"] = offtarget_top_cell_type.values
    df["offtarget_pooled_mean"] = offtarget_pooled_mean.values
    df["n_offtarget_ct_above_1"] = n_offtarget_ct_above_1.values
    df["aav_safety_ratio"] = df["target_mean_cpm"] / (df["offtarget_max_cpm"] + 0.1)

    # v1.1: CNV class + arm + arm-conditional Z
    df = annotate_cnv_class_and_arm(df)
    df["daigle_z_arm_adjusted"] = arm_adjusted_z(df, "daigle_z_raw").values

    # patient_fraction_chr7_elevated — only meaningful when target cols are patient-level
    if "malignant" in target_ct.lower() or target_ct == "malignant":
        frac, pat_ratios = patient_fraction_elevated(target_cols, col_meta, big_cpm, df)
        df["patient_fraction_chr7_elevated"] = frac.values
    else:
        df["patient_fraction_chr7_elevated"] = np.nan

    n_pre = len(df)
    df = df[df["target_mean_cpm"] >= CPM_FLOOR].copy()
    df = df.sort_values("daigle_z_arm_adjusted", ascending=False).reset_index(drop=True)
    log(f"  {target_ct}: {len(target_cols)} patient-cols; {n_pre:,}→{len(df):,} after CPM floor; "
        f"top Z_adj={df['daigle_z_arm_adjusted'].max():.1f} raw={df['daigle_z_raw'].max():.1f}")
    return df


def score_merged_target(label, target_cols, big_raw, big_cpm, col_meta, peaks):
    """For pan_myeloid etc. — targets are multiple labels pooled."""
    bg_cols = [c for c in big_cpm.columns if c not in target_cols]
    if len(target_cols) == 0 or len(bg_cols) < 2:
        return None

    target_mean = big_cpm[target_cols].mean(axis=1)
    target_detected = (big_raw[target_cols] > 0)
    n_target_patients = target_detected.sum(axis=1)
    target_cohorts_per_col = {c: col_meta[c]["cohort"] for c in target_cols}
    unique_target_cohorts = set(target_cohorts_per_col.values())

    cohort_cpm_sums = {}
    for cohort in unique_target_cohorts:
        cohort_target_cols = [c for c in target_cols if col_meta[c]["cohort"] == cohort]
        cohort_cpm_sums[cohort] = big_cpm[cohort_target_cols].sum(axis=1)
    cohort_cpm_df = pd.DataFrame(cohort_cpm_sums)
    row_sum = cohort_cpm_df.sum(axis=1).replace(0, np.nan)
    top_cohort_frac = cohort_cpm_df.max(axis=1) / row_sum
    n_target_cohorts = (cohort_cpm_df > 0).sum(axis=1)

    bg_mean = big_cpm[bg_cols].mean(axis=1)
    bg_std = big_cpm[bg_cols].std(axis=1).replace(0, np.nan)
    daigle_z = ((target_mean - bg_mean) / bg_std).fillna(0)

    # per-cell-type contamination — the "target" labels here (myeloid or malignant)
    # should be excluded from the off-target computation
    target_labels = {col_meta[c]["label"] for c in target_cols}
    other_labels = {m["label"] for c, m in col_meta.items()
                    if m["label"] not in target_labels and m["label"] not in ("", "nan")}
    per_ct_mean = {}
    for lab in other_labels:
        lab_cols = [c for c, m in col_meta.items() if m["label"] == lab]
        if lab_cols:
            per_ct_mean[lab] = big_cpm[lab_cols].mean(axis=1)
    per_ct_df = pd.DataFrame(per_ct_mean)
    offtarget_max_cpm = per_ct_df.max(axis=1)
    offtarget_top_cell_type = per_ct_df.idxmax(axis=1)
    offtarget_pooled_mean = per_ct_df.mean(axis=1)
    n_offtarget_ct_above_1 = (per_ct_df > 1.0).sum(axis=1)

    df = peaks[["chrom", "start", "end", "peak_id"]].copy()
    df["target_mean_cpm"] = target_mean.values
    df["bg_mean_cpm"] = bg_mean.values
    df["daigle_z_raw"] = daigle_z.values
    df["n_target_patients"] = n_target_patients.values
    df["n_target_cohorts"] = n_target_cohorts.values
    df["cohort_balance"] = top_cohort_frac.fillna(0).values
    df["offtarget_max_cpm"] = offtarget_max_cpm.values
    df["offtarget_top_cell_type"] = offtarget_top_cell_type.values
    df["offtarget_pooled_mean"] = offtarget_pooled_mean.values
    df["n_offtarget_ct_above_1"] = n_offtarget_ct_above_1.values
    df["aav_safety_ratio"] = df["target_mean_cpm"] / (df["offtarget_max_cpm"] + 0.1)

    # v1.1: CNV class + arm + arm-conditional Z
    df = annotate_cnv_class_and_arm(df)
    df["daigle_z_arm_adjusted"] = arm_adjusted_z(df, "daigle_z_raw").values

    if "malignant" in label.lower():
        frac, _ = patient_fraction_elevated(target_cols, col_meta, big_cpm, df)
        df["patient_fraction_chr7_elevated"] = frac.values
    else:
        df["patient_fraction_chr7_elevated"] = np.nan

    n_pre = len(df)
    df = df[df["target_mean_cpm"] >= CPM_FLOOR].copy()
    df = df.sort_values("daigle_z_arm_adjusted", ascending=False).reset_index(drop=True)
    log(f"  {label}: {len(target_cols)} target cols; {n_pre:,}→{len(df):,} after CPM floor; "
        f"top Z_adj={df['daigle_z_arm_adjusted'].max():.1f} raw={df['daigle_z_raw'].max():.1f}")
    return df


def main():
    log(f"=== score v7 start {time.strftime('%Y-%m-%d %H:%M:%S')} CPM_FLOOR={CPM_FLOOR} ===")
    peaks, big_raw, col_meta = load_v7()

    big_cpm = cpm_normalize(big_raw)
    log("CPM-normalized")

    all_labels = sorted({m["label"] for m in col_meta.values() if m["label"] not in ("", "nan")})
    log(f"labels across cohorts: {all_labels}")

    for lab in all_labels:
        out_label = "pan_malignant" if lab == "malignant" else f"per_celltype_{lab}"
        df = score_one_target(lab, big_raw, big_cpm, col_meta, peaks)
        if df is not None:
            df.to_parquet(RPT / f"{out_label}.parquet", index=False)

    # pan_myeloid: microglia + TAM + macrophage
    mye_cols = [c for c, m in col_meta.items()
                if any(k in m["label"].lower() for k in ["microglia", "tam", "macrophage"])]
    if mye_cols:
        df = score_merged_target("pan_myeloid", mye_cols, big_raw, big_cpm, col_meta, peaks)
        if df is not None:
            df.to_parquet(RPT / "pan_myeloid.parquet", index=False)

    log(f"=== score v7 end {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
