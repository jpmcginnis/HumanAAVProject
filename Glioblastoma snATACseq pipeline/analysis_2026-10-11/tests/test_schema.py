"""Schema tests — assert the deposited parquet/BED files have the columns
and dtypes documented in README.md.

Run: pytest tests/test_schema.py
"""
from pathlib import Path
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent.parent
MATRIX_DIR = HERE / "data" / "matrix"
REPORTS = HERE / "reports" / "detail"


# ---- Expected column sets ----

UNION_PROVENANCE_COLS = {"peak_id", "chrom", "start", "end", "source", "n_sources_total"}

ENHANCER_MATRIX_FIRST_COLS = {"peak_id"}  # followed by many <slug>__<cell_type> columns

REPORT_SCORE_COLS = {
    "chrom", "start", "end", "peak_id",
    "target_mean_cpm", "bg_mean_cpm",
    "daigle_z_raw", "daigle_z_arm_adjusted",
    "n_target_patients", "n_target_pseudobulks", "n_target_cohorts",
    "cohort_balance", "greenwald_fraction",
    "offtarget_max_cpm", "offtarget_top_cell_type",
    "offtarget_pooled_mean", "n_offtarget_ct_above_1", "aav_safety_ratio",
    "chromosome_arm", "chromosome_cnv_class",
}

ANNOTATED_CSV_COLS = REPORT_SCORE_COLS | {"nearest_gene", "dist_to_tss_signed", "distance_category"}


# ---- Fixtures ----

@pytest.fixture(scope="module")
def peak_provenance():
    return pd.read_parquet(MATRIX_DIR / "peak_provenance.parquet")


@pytest.fixture(scope="module")
def enhancer_matrix():
    return pd.read_parquet(MATRIX_DIR / "enhancer_candidate_matrix.parquet")


@pytest.fixture(scope="module")
def pan_malignant():
    return pd.read_parquet(REPORTS / "pan_malignant.parquet")


# ---- Schema tests ----

def test_peak_provenance_columns(peak_provenance):
    assert UNION_PROVENANCE_COLS <= set(peak_provenance.columns), (
        f"missing provenance cols: {UNION_PROVENANCE_COLS - set(peak_provenance.columns)}"
    )


def test_peak_provenance_sources_match_canonical_slugs(peak_provenance):
    slugs = set(peak_provenance["source"].unique())
    expected = {
        "greenwald_2024", "sojka_2025", "sundaram_2024", "zhao_2024",
        "mathewson_2021", "guo_2020", "guilhamon_2021", "nagaraja_2021",
        "catlas",
    }
    # All source values must be in the expected set
    unexpected = slugs - expected
    assert not unexpected, f"unexpected source slugs in provenance: {unexpected}"


def test_enhancer_matrix_has_peak_id(enhancer_matrix):
    assert "peak_id" in enhancer_matrix.columns
    assert enhancer_matrix["peak_id"].str.contains(":").all(), (
        "peak_id should look like 'chrN:start-end'"
    )


def test_enhancer_matrix_has_pseudobulk_columns(enhancer_matrix):
    pb_cols = [c for c in enhancer_matrix.columns if "__" in c]
    assert len(pb_cols) >= 60, f"expected 60+ cohort__celltype columns, got {len(pb_cols)}"
    # Every pseudobulk column name must be <slug>__<cell_type>
    for c in pb_cols:
        slug, celltype = c.split("__", 1)
        assert slug and celltype, f"malformed pseudobulk col: {c}"


def test_pan_malignant_report_schema(pan_malignant):
    missing = REPORT_SCORE_COLS - set(pan_malignant.columns)
    assert not missing, f"pan_malignant missing cols: {missing}"


def test_annotated_csvs_have_gene_cols():
    # pick one per-celltype annotated CSV and check columns
    one = HERE / "reports" / "detail" / "per_celltype" / "per_celltype_microglia_annotated.csv"
    assert one.exists(), f"missing {one}"
    df = pd.read_csv(one, nrows=5)
    missing = ANNOTATED_CSV_COLS - set(df.columns)
    assert not missing, f"annotated CSV missing cols: {missing}"


def test_n_target_patients_within_atlas_bounds(pan_malignant, enhancer_matrix):
    """Regression guard: n_target_patients must never exceed the atlas total
    (unique (cohort, patient_id) across the matrix). The reviewer caught a
    previous version where this field was actually counting pseudobulk columns
    and reported 53 against a 48-patient atlas."""
    pb_cols = [c for c in enhancer_matrix.columns if "__" in c]
    atlas_total = len({(c.split("__", 2)[0], c.split("__", 2)[1]) for c in pb_cols})
    assert pan_malignant["n_target_patients"].max() <= atlas_total, (
        f"n_target_patients {pan_malignant['n_target_patients'].max()} > atlas {atlas_total}"
    )
    # n_target_pseudobulks may legitimately exceed n_target_patients for pooled reports
    assert (pan_malignant["n_target_patients"] <= pan_malignant["n_target_pseudobulks"]).all(), (
        "n_target_patients should never exceed n_target_pseudobulks"
    )


def test_greenwald_fraction_in_range(pan_malignant):
    """greenwald_fraction must be in [0, 1]; also disclose the known-high
    median for pan_malignant so a reviewer isn't surprised."""
    gf = pan_malignant["greenwald_fraction"]
    assert gf.min() >= 0 and gf.max() <= 1
    # pan_malignant is Greenwald-weighted by construction (5 Neftel labels)
    assert gf.median() > 0.5, f"pan_malignant greenwald_fraction median {gf.median():.2f} < 0.5 — unexpected, double-check pool rule"


def test_union_bed_well_formed():
    bed = MATRIX_DIR / "union_peak_reference.bed"
    assert bed.exists()
    df = pd.read_csv(bed, sep="\t", header=None, nrows=1000,
                     names=["chrom", "start", "end", "sources"])
    assert (df["chrom"].str.startswith("chr")).all()
    assert (df["end"] > df["start"]).all()
    # widths should be reasonable (no chromosome-arm chaining)
    widths = df["end"] - df["start"]
    assert widths.max() < 50000, (
        f"peak width {widths.max()} suggests transitive chaining; mean={widths.mean():.0f}"
    )
