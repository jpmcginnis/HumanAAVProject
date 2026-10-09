"""Schema tests — assert every deposited parquet has the exact columns + dtypes
documented in SCHEMA.md.

Run: pytest tests/test_schema.py
"""
import polars as pl
from pathlib import Path
import pytest

HERE = Path(__file__).resolve().parent.parent
MATRIX_DIR = HERE / "data" / "matrix"


# ---- Expected schemas (frozen in SCHEMA.md) ----

PAN_MALIGNANT_SCORES_SCHEMA = {
    "peak_id": pl.String,
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "strength_mean": pl.Float32,
    "strength_median": pl.Float32,
    "n_patients_accessible": pl.UInt32,
    "n_patients_detected": pl.UInt32,
    "n_cohorts": pl.UInt32,
    "consistency": pl.Float64,
}

PAN_MALIGNANT_MATRIX_SCHEMA = {
    "peak_id": pl.String,
    "chrom": pl.String,
    "start": pl.Int64,
    "end": pl.Int64,
    "cohort": pl.String,
    "patient_id": pl.String,
    "modality": pl.String,
    "n_cells_malignant": pl.Int64,
    "n_accessible": pl.Int64,
    "frac_accessible": pl.Float32,
    "mean_counts": pl.Float32,
}

CANDIDATE_SCORES_SCHEMA = {
    "peak_id": pl.String,
    "cell_type": pl.String,
    "strength": pl.Float32,
    "n_patients": pl.UInt32,
    "consistency": pl.Float64,
    "n_cohorts": pl.UInt32,
    "selectivity": pl.Float32,
    "posterior_mean": pl.Float64,
    "posterior_ci_lo": pl.Float64,
    "posterior_ci_hi": pl.Float64,
    "posterior_p_specific": pl.Float64,
    "pass2_was_run": pl.Boolean,
    "composite_score": pl.Float64,
}


def _assert_schema(parquet_path: Path, expected: dict):
    assert parquet_path.exists(), f"missing: {parquet_path}"
    df = pl.read_parquet(parquet_path, n_rows=1)
    for col, dtype in expected.items():
        assert col in df.columns, f"{parquet_path.name}: missing column '{col}'"
        assert df.schema[col] == dtype, (
            f"{parquet_path.name}: column '{col}' is {df.schema[col]}, "
            f"expected {dtype}"
        )
    extra = set(df.columns) - set(expected.keys())
    assert not extra, f"{parquet_path.name}: unexpected columns {extra}"


def test_pan_malignant_scores_schema():
    _assert_schema(MATRIX_DIR / "pan_malignant_scores_v4.parquet", PAN_MALIGNANT_SCORES_SCHEMA)


def test_pan_malignant_matrix_schema():
    _assert_schema(MATRIX_DIR / "pan_malignant_matrix_v4.parquet", PAN_MALIGNANT_MATRIX_SCHEMA)


def test_candidate_scores_schema():
    _assert_schema(MATRIX_DIR / "candidate_scores.parquet", CANDIDATE_SCORES_SCHEMA)


def test_pan_malignant_scores_row_count():
    df = pl.read_parquet(MATRIX_DIR / "pan_malignant_scores_v4.parquet")
    # CATLAS atlas is 544,735 peaks; v4 should score every one
    assert df.height == 544_735


def test_pan_malignant_matrix_patient_count():
    df = pl.read_parquet(MATRIX_DIR / "pan_malignant_matrix_v4.parquet", columns=["cohort","patient_id"])
    # Post-CNV-fix: 45 distinct patients in 8 cohorts
    n_patients = df.unique().height
    assert n_patients == 45, f"expected 45 (cohort, patient_id), got {n_patients}"
    n_cohorts = df["cohort"].n_unique()
    assert n_cohorts == 8, f"expected 8 cohorts, got {n_cohorts}"


def test_cohorts_match_attribution_csv():
    # Every cohort in the matrix must have a row in cohort_attribution.csv
    import csv
    with open(HERE / "cohort_attribution.csv") as f:
        attributed = {row["cohort_slug"] for row in csv.DictReader(f)}
    matrix_cohorts = set(
        pl.read_parquet(MATRIX_DIR / "pan_malignant_matrix_v4.parquet", columns=["cohort"])["cohort"].unique().to_list()
    )
    missing = matrix_cohorts - attributed
    assert not missing, f"cohorts in matrix but not in cohort_attribution.csv: {missing}"
