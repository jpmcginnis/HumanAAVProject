#!/usr/bin/env python3
"""Run the atlas scoring chain against ONE synthetic cohort parquet. Mirrors the
real atlas pipeline (nb_specificity_recompute → nb_counts_recompute → fix_strength_to_gbm)
but inlined and self-contained so the smoke test doesn't need the full union-peak
bed or the real multi-cohort scored parquet.

usage: python3 run_scoring_on_synth.py IN.parquet OUT.parquet
"""
import sys
import polars as pl
import numpy as np

MIN_CELLS = 100
EPS = 1e-6
REF_WIDTH = 500

in_pq, out_pq = sys.argv[1], sys.argv[2]
d = pl.read_parquet(in_pq)

# All synthetic donors are PASS. Evaluability gate: n_cells_total_direct ≥ MIN_CELLS.
d_eval = d.filter(pl.col("n_cells_total_direct") >= MIN_CELLS)
donors_per_class = (d_eval.group_by("major_class")
                    .agg(pl.col("donor_id").n_unique().alias("n_evaluable_class")))

# Per-(peak, class) sum_frac + n_donors with signal
peak_class = (d_eval.filter(pl.col("fraction_open_direct") > 0)
                    .group_by(["peak_num","major_class"])
                    .agg([
                        pl.col("fraction_open_direct").sum().alias("sum_frac"),
                        pl.col("donor_id").n_unique().alias("n_donors"),
                    ])
                    .join(donors_per_class, on="major_class", how="left"))
# mean_frac uses n_evaluable_class denominator (GBM-matched)
peak_class = peak_class.with_columns([
    (pl.col("sum_frac") / pl.col("n_evaluable_class").cast(pl.Float64)).alias("mean_frac"),
])

# Pivot to compute per-peak log2FC vs best off-target class
pivot = peak_class.pivot(values="mean_frac", index="peak_num", on="major_class").fill_null(0.0)
class_cols = [c for c in pivot.columns if c != "peak_num"]
arr = pivot.select(class_cols).to_numpy()
n_peaks, n_classes = arr.shape

log2fc = np.zeros_like(arr)
for j in range(n_classes):
    masked = arr.copy()
    masked[:, j] = -np.inf
    off = masked.max(axis=1)
    off = np.where(np.isinf(off), 0.0, off)
    log2fc[:, j] = np.log2((arr[:, j] + EPS) / (off + EPS))

peak_nums = pivot["peak_num"].to_numpy()
spec_rows = []
for j, cls in enumerate(class_cols):
    for i, p in enumerate(peak_nums):
        spec_rows.append({"peak_num": int(p), "major_class": cls, "log2FC_fresh": float(log2fc[i, j])})
spec_df = pl.DataFrame(spec_rows)

out = peak_class.join(spec_df, on=["peak_num","major_class"], how="left")

# Synthetic peak_width (not stored in synth parquet); use a fixed value so strength is deterministic
out = out.with_columns([
    pl.lit(500).alias("peak_width"),
])
# strength_wn_gbm = mean_frac × (REF_WIDTH / peak_width)
out = out.with_columns(
    (pl.col("mean_frac") * (REF_WIDTH / pl.col("peak_width").cast(pl.Float64))).alias("strength_wn_gbm")
)
out = out.with_columns([
    (pl.col("n_donors").cast(pl.Float64) / pl.col("n_evaluable_class").cast(pl.Float64))
      .alias("patient_detection_rate_gbm"),
])

out.write_parquet(out_pq)
print(f"wrote {out_pq} ({out.shape[0]} rows)")
