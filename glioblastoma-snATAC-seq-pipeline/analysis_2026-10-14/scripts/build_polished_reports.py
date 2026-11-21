"""Polished top-100 distal reports.

Data dependency: expects the annotated per-report CSVs (produced by
scripts/build_reports.py from the Zenodo-deposited parquets:
https://doi.org/10.5281/zenodo.23252324). Run build_reports.py first.

Reads the annotated per-report CSVs from reports/detail/ and reports/detail/per_celltype/,
applies the "polished" filters:
  - distance_category != "promoter"      (distal / proximal / far-distal only)
  - target_mean_cpm >= 1.0               (minimum signal floor)
  - n_target_patients >= 2               (patient-level replication, inherited)
then takes the top-100 by daigle_z per report.

Outputs to reports/ root (not detail/):
  <report>_top100.csv       — per report
  top100_distal_combined.csv/xlsx/parquet — all reports combined, one tab per report
"""
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
REPORTS = ROOT / "reports"
DETAIL = REPORTS / "detail"
PER_CT = DETAIL / "per_celltype"

CPM_FLOOR = 1.0
TOP_N = 100


def label_from_path(p):
    name = p.stem.replace("_annotated", "")
    return name.replace("per_celltype_", "")


# Exclude Neftel-state subtypes from the main combined/polished reports — they
# have per_celltype_<state>.parquet in detail/ for drill-down, but pan_malignant
# (which merges all malignant_* columns) is what goes in the top-level HTML.
NEFTEL_SUBTYPE_LABELS = {"malignant_AC", "malignant_MES", "malignant_OPC",
                        "malignant_NPC", "malignant"}


def main():
    all_per_ct = sorted(PER_CT.glob("per_celltype_*_annotated.csv"))
    # Filter out Neftel subtypes from the main combined — they're available
    # in detail/ but not in the top-level polished xlsx/html.
    per_ct_for_combined = [p for p in all_per_ct
                           if label_from_path(p) not in NEFTEL_SUBTYPE_LABELS]
    sources = per_ct_for_combined + sorted(DETAIL.glob("pan_*_annotated.csv"))
    rows = []
    for p in sources:
        label = label_from_path(p)
        df = pd.read_csv(p)
        n_pre = len(df)
        df = df[(df["distance_category"] != "promoter") & (df["target_mean_cpm"] >= CPM_FLOOR)].copy()
        n_post = len(df)
        df = df.sort_values("daigle_z_arm_adjusted", ascending=False).head(TOP_N).reset_index(drop=True)
        df.insert(0, "rank", range(1, len(df) + 1))
        df.insert(0, "report", label)
        # Per-report polished CSV
        df.drop(columns=["report"]).to_csv(REPORTS / f"{label}_top100.csv", index=False)
        rows.append(df)
        print(f"  {label}: {n_pre:,} → {n_post:,} after distal+CPM≥{CPM_FLOOR} → top {len(df)}")

    combined = pd.concat(rows, ignore_index=True)
    combined.to_parquet(REPORTS / "top100_distal_combined.parquet", index=False)
    combined.to_csv(REPORTS / "top100_distal_combined.csv", index=False)

    # v1.1 — also build top-10 per chromosome arm (per report) for CNV-diversity view
    if "chromosome_arm" in combined.columns:
        rows_arm = []
        for p in sources:
            label = label_from_path(p)
            df = pd.read_csv(p)
            df = df[(df["distance_category"] != "promoter") & (df["target_mean_cpm"] >= CPM_FLOOR)].copy()
            df = df.sort_values("daigle_z_arm_adjusted", ascending=False)
            per_arm = df.groupby("chromosome_arm").head(10).copy()
            per_arm.insert(0, "report", label)
            rows_arm.append(per_arm)
        arm_combined = pd.concat(rows_arm, ignore_index=True)
        arm_combined.to_csv(REPORTS / "top10_per_arm_combined.csv", index=False)
        with pd.ExcelWriter(REPORTS / "top10_per_arm_combined.xlsx", engine="openpyxl") as xw:
            for lab in arm_combined["report"].unique():
                arm_combined[arm_combined["report"] == lab].to_excel(xw, sheet_name=lab[:31], index=False)
        print(f"wrote top10_per_arm_combined ({len(arm_combined)} rows)")
    with pd.ExcelWriter(REPORTS / "top100_distal_combined.xlsx", engine="openpyxl") as xw:
        summary = combined.groupby("report").agg(
            n=("rank", "count"),
            median_z=("daigle_z_arm_adjusted", "median"),
            median_target_cpm=("target_mean_cpm", "median"),
            median_n_patients=("n_target_patients", "median"),
            median_cohort_balance=("cohort_balance", "median"),
            median_aav_safety=("aav_safety_ratio", "median"),
            median_offtarget_max=("offtarget_max_cpm", "median"),
        ).reset_index()
        summary.to_excel(xw, sheet_name="summary", index=False)
        for label in combined["report"].unique():
            combined[combined["report"] == label].to_excel(xw, sheet_name=label[:31], index=False)
    print(f"\nwrote combined {combined.shape}")
    print(f"  {REPORTS / 'top100_distal_combined.csv'}")
    print(f"  {REPORTS / 'top100_distal_combined.xlsx'}")
    print(f"  {REPORTS / 'top100_distal_combined.parquet'}")
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
