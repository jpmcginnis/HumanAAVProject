"""Filter tests — assert the replication filter + distal-selective filter
produce the counts the reports claim.

Run: pytest tests/test_filters.py
"""
from pathlib import Path
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent.parent
REPORTS = HERE / "reports" / "detail"
PER_CT = REPORTS / "per_celltype"


@pytest.fixture(scope="module")
def pan_malignant():
    return pd.read_parquet(REPORTS / "pan_malignant.parquet")


def test_patient_level_replication(pan_malignant):
    """Primary replication is patient-level in v1. Every scored peak passes
    n_target_patients >= 2 by construction of score_v7.py's MIN_PATIENTS gate."""
    assert (pan_malignant["n_target_patients"] >= 2).all(), (
        f"found {(pan_malignant['n_target_patients'] < 2).sum()} peaks with <2 target patients"
    )
    # n_target_cohorts can be 1 if all target patients come from one cohort —
    # that's a batch-dominated peak, flagged by cohort_balance ≈ 1.
    single_cohort = (pan_malignant["n_target_cohorts"] == 1).sum()
    multi_cohort = (pan_malignant["n_target_cohorts"] >= 2).sum()
    assert multi_cohort > 0, "no multi-cohort peaks — likely a scoring bug"
    assert single_cohort < len(pan_malignant), "every peak looks single-cohort — scoring bug"


def test_distal_top30_bed_excludes_promoters():
    # Pick one per-celltype top30 distal BED and verify every row's name has a kb distance >= 2kb
    import re
    bed = PER_CT / "per_celltype_microglia_top30_distal.bed"
    assert bed.exists()
    df = pd.read_csv(bed, sep="\t", header=None, names=["chrom", "start", "end", "name", "score", "strand"])
    assert len(df) == 30, f"expected 30 peaks, got {len(df)}"
    # name format: GENE_±XkB_peakid  (e.g., GAD2_+1.2kb_chr10:...)
    bad = []
    for n in df["name"]:
        m = re.search(r"([+-]?[\d.]+)kb_chr", n)
        if m:
            kb = float(m.group(1))
            if abs(kb) < 2:
                bad.append(n)
    assert not bad, f"top30 distal contains promoter-adjacent peaks: {bad[:3]}"


def test_top50_combined_has_800_rows():
    df = pd.read_parquet(REPORTS / "top50_combined.parquet")
    # 16 reports × 50 = 800
    #   10 non-malignant per_celltype: GABA_neuron, OPC, TAM, T_cell, astrocyte,
    #       endothelial, microglia, neuron, oligodendrocyte, unassigned
    #   4 Neftel-state per_celltype (detail-only): malignant_AC/MES/OPC/NPC
    #   2 pooled: pan_malignant, pan_myeloid
    assert len(df) == 800, f"expected 800, got {len(df)}"
    counts = df.groupby("report").size()
    assert (counts == 50).all(), f"not all reports have 50 rows: {counts.to_dict()}"
    # Replication filter applied at patient level, not cohort level.
    assert (df["n_target_patients"] >= 2).all(), (
        f"found {(df['n_target_patients'] < 2).sum()} rows below patient replication threshold"
    )


def test_annotated_distance_categories_consistent():
    # Load one annotated CSV, verify distance_category matches |dist_to_tss_signed|
    import sys
    sys.path.insert(0, str(HERE / "scripts"))
    from build_reports import PROMOTER_WINDOW, PROXIMAL_WINDOW, DISTAL_WINDOW

    df = pd.read_csv(PER_CT / "per_celltype_microglia_annotated.csv", nrows=1000)
    abs_d = df["dist_to_tss_signed"].abs().values
    expected = []
    for d in abs_d:
        if d <= PROMOTER_WINDOW: expected.append("promoter")
        elif d <= PROXIMAL_WINDOW: expected.append("proximal")
        elif d <= DISTAL_WINDOW: expected.append("distal")
        else: expected.append("far-distal")
    mismatches = (df["distance_category"].values != expected).sum()
    assert mismatches == 0, f"{mismatches} rows have wrong distance_category"
