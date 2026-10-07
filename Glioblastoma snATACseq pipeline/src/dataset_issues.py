"""
Audit the matrix for coverage gaps that would force prospective data generation.
Output verdict: "existing public data suffices" vs "we need our own cohort".

Checks:
    - min cohorts per cell type (default 3)
    - min patients per cell type (default 15)
    - min cells per cell type per patient (default 20)
    - Multiome availability per cell type (needed for enhancer→gene linkage)
    - malignant Neftel state coverage (NPC/OPC/AC/MES-like each ≥10 patients)
"""

from __future__ import annotations
import argparse, sys, yaml
from pathlib import Path


def log(msg): print(f"[issues] {msg}", flush=True)


def main() -> int:
    import polars as pl

    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    with open("config/pipeline.yaml") as f:
        checks = yaml.safe_load(f)["dataset_issues_checks"]

    cfg = {}
    for check in checks:
        if isinstance(check, dict):
            cfg.update(check)
    min_cohorts = cfg.get("min_cohorts_for_tme_call", 3)
    min_patients = cfg.get("min_patients_per_celltype", 15)
    min_cells = cfg.get("min_cells_per_celltype_per_patient", 20)
    needs_multiome = cfg.get("require_multiome_for_enhancer_to_gene_linkage", True)
    neftel_targets = cfg.get("acceptable_malignant_subtype_coverage", {})

    log(f"Loading matrix: {args.src}")
    m = pl.read_parquet(args.src)

    issues = []

    log("Checking cohort/patient coverage per cell type")
    per_ct = (m.group_by("cell_type")
                .agg([
                    pl.col("cohort").n_unique().alias("n_cohorts"),
                    pl.col("patient_id").n_unique().alias("n_patients"),
                    pl.col("modality").is_in(["multiome"]).any().alias("has_multiome"),
                ]))
    for row in per_ct.iter_rows(named=True):
        ct = row["cell_type"]
        if row["n_cohorts"] < min_cohorts:
            issues.append({"cell_type": ct, "issue": "insufficient_cohorts",
                           "observed": row["n_cohorts"], "required": min_cohorts, "severity": "high"})
        if row["n_patients"] < min_patients:
            issues.append({"cell_type": ct, "issue": "insufficient_patients",
                           "observed": row["n_patients"], "required": min_patients, "severity": "high"})
        if needs_multiome and not row["has_multiome"]:
            issues.append({"cell_type": ct, "issue": "no_multiome_coverage",
                           "observed": 0, "required": 1, "severity": "medium"})

    log("Checking Neftel malignant state coverage")
    for state, min_pat in neftel_targets.items():
        obs = (m.filter(pl.col("cell_type") == state)
                 .select(pl.col("patient_id").n_unique())
                 .item() or 0)
        if obs < min_pat:
            issues.append({"cell_type": state, "issue": "insufficient_malignant_state_coverage",
                           "observed": obs, "required": min_pat, "severity": "high"})

    log("Checking cells-per-celltype-per-patient")
    per_pt = (m.group_by(["cell_type", "patient_id", "cohort"])
                .agg(pl.col("n_cells").max().alias("cells_in_bin")))
    weak = per_pt.filter(pl.col("cells_in_bin") < min_cells)
    weak_ct = (weak.group_by("cell_type")
                   .agg(pl.col("patient_id").n_unique().alias("n_underpowered_patients")))
    for row in weak_ct.iter_rows(named=True):
        if row["n_underpowered_patients"] >= 3:
            issues.append({"cell_type": row["cell_type"], "issue": "chronic_undercapture",
                           "observed": row["n_underpowered_patients"], "required": 0, "severity": "medium"})

    df = pl.DataFrame(issues) if issues else pl.DataFrame(schema={
        "cell_type": pl.Utf8, "issue": pl.Utf8, "observed": pl.Int64, "required": pl.Int64, "severity": pl.Utf8
    })

    high_count = df.filter(pl.col("severity") == "high").height if df.height > 0 else 0
    log("")
    log("=== VERDICT ===")
    if high_count == 0:
        log("  Existing public data suffices. R01 can be written for VALIDATION, not discovery.")
        log("  → Aim 3 pitches as prospective validation of candidate enhancers in slice-culture AAV screen.")
    else:
        log(f"  {high_count} HIGH-SEVERITY coverage issues detected.")
        log("  Recommend prospective n=24 Multiome cohort to close the gaps.")
        log("  → Aim 3 pitches as DISCOVERY: generate cohort, identify candidates, validate in slice-culture AAV screen.")
        for row in df.filter(pl.col("severity") == "high").iter_rows(named=True):
            log(f"    - {row['cell_type']:30s} {row['issue']:35s} obs={row['observed']} < req={row['required']}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(args.out, compression="zstd")
    log(f"Wrote → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
