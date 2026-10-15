"""Computation tests — distance_category boundaries, nearest-TSS correctness,
Daigle Z formula.

Run: pytest tests/test_computations.py
"""
import sys
import numpy as np
import pandas as pd
import pytest

HERE_SCRIPTS = __import__("pathlib").Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(HERE_SCRIPTS))

# Import the real functions used by the report builder so tests catch drift.
from build_reports import nearest_tss, PROMOTER_WINDOW, PROXIMAL_WINDOW, DISTAL_WINDOW


def _category(abs_d):
    if abs_d <= PROMOTER_WINDOW:
        return "promoter"
    if abs_d <= PROXIMAL_WINDOW:
        return "proximal"
    if abs_d <= DISTAL_WINDOW:
        return "distal"
    return "far-distal"


# ---- distance_category boundaries ----

def test_distance_category_boundaries():
    assert _category(0) == "promoter"
    assert _category(PROMOTER_WINDOW) == "promoter"          # inclusive upper
    assert _category(PROMOTER_WINDOW + 1) == "proximal"
    assert _category(PROXIMAL_WINDOW) == "proximal"
    assert _category(PROXIMAL_WINDOW + 1) == "distal"
    assert _category(DISTAL_WINDOW) == "distal"
    assert _category(DISTAL_WINDOW + 1) == "far-distal"
    assert _category(1_000_000) == "far-distal"


# ---- nearest_tss correctness ----

def test_nearest_tss_picks_closest():
    # Build a tiny tss_by_chrom with two genes on chr1
    tss_df = pd.DataFrame({"chrom": ["chr1", "chr1"], "tss": [1000, 10000], "gene": ["GENE_A", "GENE_B"]})
    tss_by_chrom = {"chrom": tss_df, "chr1": tss_df.sort_values("tss").reset_index(drop=True)}
    # Peak 1: center 1200 → nearest is GENE_A (TSS 1000), signed distance +200
    # Peak 2: center 5000 → nearest is GENE_A (TSS 1000, dist 4000) vs GENE_B (TSS 10000, dist -5000) → A
    # Peak 3: center 8000 → nearest is GENE_B (TSS 10000, signed -2000)
    peaks = pd.DataFrame({
        "chrom": ["chr1", "chr1", "chr1"],
        "start": [1000, 4750, 7750],
        "end":   [1400, 5250, 8250],
    })
    annotated = nearest_tss(peaks, tss_by_chrom)
    assert annotated.loc[0, "nearest_gene"] == "GENE_A"
    assert annotated.loc[0, "dist_to_tss_signed"] == 200
    assert annotated.loc[1, "nearest_gene"] == "GENE_A"
    assert annotated.loc[1, "dist_to_tss_signed"] == 4000
    assert annotated.loc[2, "nearest_gene"] == "GENE_B"
    assert annotated.loc[2, "dist_to_tss_signed"] == -2000


def test_nearest_tss_handles_missing_chrom():
    tss_by_chrom = {"chr1": pd.DataFrame({"chrom": ["chr1"], "tss": [1000], "gene": ["GENE_A"]})}
    peaks = pd.DataFrame({"chrom": ["chrX"], "start": [500], "end": [1000]})
    out = nearest_tss(peaks, tss_by_chrom)
    # No TSS data for chrX → gene stays "", dist stays 0
    assert out.loc[0, "nearest_gene"] == ""
    assert out.loc[0, "dist_to_tss_signed"] == 0


# ---- Daigle Z formula ----

def test_daigle_z_formula():
    # target >> background → large positive Z
    target = pd.Series([10.0, 2.0, 0.5])
    bg = pd.DataFrame({"a": [1.0, 1.0, 1.0], "b": [2.0, 2.0, 0.4], "c": [0.5, 0.5, 0.6]})
    mu = bg.mean(axis=1)
    sd = bg.std(axis=1)
    z = (target - mu) / sd
    # First peak: target=10, bg~1.17 ± 0.76 → z ≈ 11.6
    assert z.iloc[0] > 10
    # Last peak: target 0.5, bg ~0.5 → z near 0
    assert abs(z.iloc[2]) < 2
