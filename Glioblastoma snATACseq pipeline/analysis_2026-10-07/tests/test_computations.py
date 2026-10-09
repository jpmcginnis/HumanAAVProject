"""Unit tests for the computational helpers in the report builders:
distance_category, Daigle Z-score, selectivity, aav_score.

These test the inline helper functions by replicating their logic against
hand-crafted toy inputs with expected outputs.

Run: pytest tests/test_computations.py
"""
import math
import pandas as pd
import pytest


# ---- distance_category (copied from the three report builders) ----
# All three builders use the identical categorization.

def _dist_cat(kb):
    if pd.isna(kb): return "unknown"
    if kb < 2:    return "promoter"
    if kb < 10:   return "near"
    if kb < 100:  return "distal"
    if kb < 500:  return "far-distal"
    return "gene-desert"


def test_distance_category_boundaries():
    # Edge cases around each boundary
    assert _dist_cat(0) == "promoter"
    assert _dist_cat(1.9) == "promoter"
    assert _dist_cat(2.0) == "near"          # 2 kb floor is excluded by the filter, labeled "near"
    assert _dist_cat(9.9) == "near"
    assert _dist_cat(10.0) == "distal"
    assert _dist_cat(99.9) == "distal"
    assert _dist_cat(100.0) == "far-distal"  # old upper cap, now just an annotation boundary
    assert _dist_cat(499.9) == "far-distal"
    assert _dist_cat(500.0) == "gene-desert"
    assert _dist_cat(5000.0) == "gene-desert"


def test_distance_category_nan():
    assert _dist_cat(float("nan")) == "unknown"


# ---- Daigle Z-score ----
# Z = (strength_target - mean_other) / sd_other, with sd clipped at 0.001

def _daigle_z(target, others_mean, others_sd):
    return (target - others_mean) / max(others_sd, 0.001)


def test_daigle_z_basic():
    # Target 0.5, others mean 0.1 sd 0.1 -> Z = 4
    assert _daigle_z(0.5, 0.1, 0.1) == pytest.approx(4.0)


def test_daigle_z_passes_at_2():
    # Target 0.3 vs others mean 0.1 sd 0.1 -> Z = 2.0 exactly
    assert _daigle_z(0.3, 0.1, 0.1) == pytest.approx(2.0)


def test_daigle_z_negative_target():
    # Target below others mean -> negative Z
    assert _daigle_z(0.05, 0.1, 0.1) == pytest.approx(-0.5)


def test_daigle_z_sd_floor_prevents_inf():
    # Zero sd should not blow up — clipped to 0.001
    z = _daigle_z(0.5, 0.1, 0.0)
    assert math.isfinite(z)
    assert z == pytest.approx(0.4 / 0.001)


def test_daigle_z_criterion_threshold():
    # Daigle criterion: Z >= 2 is the Armamentarium specificity threshold.
    # (0.3 - 0.1) / 0.1 is 1.9999...998 in IEEE-754 so we use pytest.approx for the pass case.
    assert _daigle_z(0.3, 0.1, 0.1) == pytest.approx(2.0)
    assert _daigle_z(0.29, 0.1, 0.1) < 2.0   # fail
    # Clear pass: 3 sigma above the mean
    assert _daigle_z(0.4, 0.1, 0.1) >= 2.0


# ---- Selectivity (pan-malignant: strength_mean / strength_max_nonmalignant) ----

def _selectivity(target, max_other):
    return target / max(max_other, 0.01)


def test_selectivity_above_1():
    assert _selectivity(0.4, 0.2) == pytest.approx(2.0)


def test_selectivity_below_1_tumor_abundant_type():
    # Astrocyte in GBM has selectivity ~0.8 because reactive glia overlap with
    # malignant signal
    assert _selectivity(0.08, 0.1) == pytest.approx(0.8)


def test_selectivity_denominator_floor():
    # Peak absent from non-target cell types -> denominator floored at 0.01
    assert _selectivity(0.5, 0.0) == pytest.approx(50.0)
    assert _selectivity(0.5, 0.005) == pytest.approx(50.0)


# ---- aav_score (pan-malignant composite) ----
# aav_score = 0.35*consistency + 0.25*(n_cohorts/8) + 0.15*(n_patients/45 clipped <=1)
#           + 0.15*(selectivity/10 clipped <=1) + 0.10*(strength/0.5 clipped <=1)

def _aav_score(consistency, n_cohorts, n_patients, selectivity, strength,
               n_patients_denom=45):
    return (
        consistency * 0.35
        + float(n_cohorts) / 8.0 * 0.25
        + min(1.0, float(n_patients) / n_patients_denom) * 0.15
        + min(1.0, min(selectivity, 10) / 10) * 0.15
        + min(1.0, min(strength, 0.5) / 0.5) * 0.10
    )


def test_aav_score_bounded_in_zero_one():
    assert _aav_score(0.0, 0, 0, 0.0, 0.0) == 0.0
    perfect = _aav_score(1.0, 8, 45, 10.0, 0.5)
    assert perfect == pytest.approx(1.0)
    # Clipping prevents scores >1
    over = _aav_score(1.0, 8, 100, 50.0, 1.0)
    assert over == pytest.approx(1.0)


def test_aav_score_weights_consistency_most():
    # Peak with high consistency only
    high_consistency = _aav_score(1.0, 4, 20, 1.0, 0.1)
    # Peak with high strength only
    high_strength = _aav_score(0.0, 4, 20, 1.0, 0.5)
    # Consistency weighted 0.35 vs strength 0.10 => consistency-only should win
    assert high_consistency > high_strength


def test_aav_score_approximate_realistic():
    # From v3 report top-10: CPNE4 ~0.488 ballpark
    # consistency 0.27, n_cohorts 7, n_patients 37, selectivity 2.17, strength 0.124
    score = _aav_score(0.27, 7, 37, 2.17, 0.124)
    # Just sanity: score should be in a reasonable range
    assert 0.3 < score < 0.6
