"""Filter tests — assert the AAV-clean pool filters produce sane counts
against the live parquet data.

Run: pytest tests/test_filters.py
"""
import polars as pl
import pandas as pd
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent.parent
MATRIX_DIR = HERE / "data" / "matrix"


@pytest.fixture(scope="module")
def pan_scores():
    return pl.read_parquet(MATRIX_DIR / "pan_malignant_scores_v4.parquet").to_pandas()


@pytest.fixture(scope="module")
def candidate_scores():
    return pl.read_parquet(MATRIX_DIR / "candidate_scores.parquet").to_pandas()


def test_pan_malignant_pool_size(pan_scores):
    # n_cohorts >= 4 is the primary pool filter. Documented in README: 170,673 peaks.
    pool = pan_scores[pan_scores["n_cohorts"] >= 4]
    assert pool.shape[0] == 170_673, f"pool size changed: {pool.shape[0]} (expected 170,673)"


def test_pan_malignant_filters_match_report(pan_scores, candidate_scores):
    # Replicate the v3 clean-pool filter:
    # distal (|TSS|>=2kb, no upper cap) + selective (>=2x non-mal) + non-chr7 + non-housekeeping.
    # Skip the nearest-gene annotation step (which requires refgene_hg38.bed and is tested in test_report_smoke).
    pool = pan_scores[pan_scores["n_cohorts"] >= 4].copy()

    # Compute selectivity like the builder does
    ct_nonmal = candidate_scores[~candidate_scores["cell_type"].isin(["malignant_unresolved", "unassigned"])]
    nonmal_strength = ct_nonmal.groupby("peak_id")["strength"].max().rename("nonmal_max").reset_index()
    pool = pool.merge(nonmal_strength, on="peak_id", how="left")
    pool["nonmal_max"] = pool["nonmal_max"].fillna(0.0)
    pool["selectivity"] = pool["strength_mean"] / pool["nonmal_max"].clip(lower=0.01)

    # Filter: selectivity >= 2.0, non-chr7
    # (not applying distal filter here because that needs nearest-gene annotation
    #  — tested separately in test_report_smoke)
    filtered = pool[(pool["selectivity"] >= 2.0) & (pool["chrom"] != "chr7")]
    # Should be in the hundreds to low thousands range. Verify not empty + not everything.
    assert 100 < filtered.shape[0] < 10_000, (
        f"selectivity+non-chr7 pool size: {filtered.shape[0]} (expected hundreds-to-thousands)"
    )


def test_pan_malignant_consistency_range(pan_scores):
    # consistency is a fraction in [0, 1]
    pool = pan_scores[pan_scores["n_cohorts"] >= 4]
    assert pool["consistency"].min() >= 0.0
    assert pool["consistency"].max() <= 1.0


def test_pan_malignant_n_cohorts_bounded(pan_scores):
    assert pan_scores["n_cohorts"].min() >= 0
    assert pan_scores["n_cohorts"].max() <= 8   # 8 cohorts total


def test_pan_malignant_n_patients_bounded(pan_scores):
    # n_patients_detected cannot exceed 45 (post-fix pool size)
    assert pan_scores["n_patients_detected"].max() <= 45
    # n_patients_accessible cannot exceed n_patients_detected
    assert (pan_scores["n_patients_accessible"] <= pan_scores["n_patients_detected"]).all()


def test_strength_range(pan_scores):
    # strength_mean is frac_accessible in [0, 1]
    assert pan_scores["strength_mean"].min() >= 0.0
    assert pan_scores["strength_mean"].max() <= 1.0
    assert pan_scores["strength_median"].min() >= 0.0
    assert pan_scores["strength_median"].max() <= 1.0


def test_candidate_scores_cell_types(candidate_scores):
    # Expected 11 cell types in the atlas
    expected = {
        "astrocyte", "endothelial", "GABA_neuron", "malignant_unresolved",
        "microglia", "neuron", "oligodendrocyte", "OPC", "T_cell", "TAM",
        "unassigned",
    }
    actual = set(candidate_scores["cell_type"].unique())
    assert actual == expected, f"cell_types diverged: {actual ^ expected}"


def test_pass2_was_run_count(candidate_scores):
    # Pass-2 posterior was run on top 2000 per cell type -> 22,000 rows
    pass2 = candidate_scores[candidate_scores["pass2_was_run"]]
    assert pass2.shape[0] == 22_000
    # Every pass2 row should have a posterior_mean
    assert pass2["posterior_mean"].notna().all()
