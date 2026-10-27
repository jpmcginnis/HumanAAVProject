#!/usr/bin/env python3
"""Patch replicability_major_v2.4_fresh_counts.parquet so strength_wn_gbm matches
GBM atlas semantics EXACTLY:

  strength_wn_gbm = (sum_frac_at_this_peak_across_evaluable_donors / n_evaluable_class)
                    × (500 / peak_width)

Donors evaluable for a class but with 0 cells accessing the peak contribute 0 to the
numerator (not NaN). Previously the denominator was n_evaluable_fresh (donors with
signal at the peak), which inflated scores on rare peaks.

Also updates mean_frac_fresh, strength_wn_gbm_fresh with the same correction.
"""
import polars as pl
from pathlib import Path

AT = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0")
IN = AT / "replicability_major_v2.4_fresh_counts.parquet"
OUT = AT / "replicability_major_v2.5_gbm_strength.parquet"
REF_WIDTH = 500

print(f"loading {IN}", flush=True)
d = pl.read_parquet(IN)
print(f"  rows: {d.shape[0]:,}", flush=True)

# `sum_frac` was dropped from spec recompute output but we have:
#   mean_frac_fresh = sum_frac / n_evaluable_fresh (old divisor)
# So: sum_frac = mean_frac_fresh * n_evaluable_fresh
# New mean = sum_frac / n_evaluable_class
# New strength = new_mean * (500/peak_width)
d = d.with_columns([
    (pl.col("mean_frac_fresh").cast(pl.Float64) * pl.col("n_evaluable_fresh").cast(pl.Float64)).alias("_sum_frac"),
])
d = d.with_columns([
    pl.when(pl.col("n_evaluable_class") > 0)
      .then(pl.col("_sum_frac") / pl.col("n_evaluable_class").cast(pl.Float64))
      .otherwise(0.0)
      .alias("mean_frac_fresh_gbm"),
])
d = d.with_columns([
    pl.when(pl.col("peak_width") > 0)
      .then(pl.col("mean_frac_fresh_gbm") * (REF_WIDTH / pl.col("peak_width").cast(pl.Float64)))
      .otherwise(0.0)
      .alias("strength_wn_gbm_new"),
])
# Replace the stored columns
d = d.with_columns([
    pl.col("mean_frac_fresh_gbm").alias("mean_frac_fresh"),
    pl.col("strength_wn_gbm_new").alias("strength_wn_gbm"),
    pl.col("strength_wn_gbm_new").alias("strength_wn_gbm_fresh"),
    pl.col("mean_frac_fresh_gbm").alias("mean_frac"),  # mirror for HTML
])
d = d.drop(["_sum_frac", "mean_frac_fresh_gbm", "strength_wn_gbm_new"])

# Sanity: top Microglia by new strength
print("\nTop 10 Microglia by GBM-fixed strength:", flush=True)
print(d.filter(pl.col("major_class")=="Microglia")
       .sort("strength_wn_gbm", descending=True)
       .head(10)
       .select(["peak_id","n_donors","n_evaluable_class","mean_frac","strength_wn_gbm","patient_detection_rate_gbm"]))

d.write_parquet(OUT)
print(f"\nwrote {OUT} ({OUT.stat().st_size/1e9:.2f} GB)", flush=True)
