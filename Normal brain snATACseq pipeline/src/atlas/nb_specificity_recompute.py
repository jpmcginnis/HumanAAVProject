#!/usr/bin/env python3
"""v2.3 specificity recompute on the 340-donor PASS pool from metadata_FINAL_v7.

Fixes two bugs the audit surfaced:
  1. stale v2.0 specificity columns (never recomputed after nnz filter)
  2. disease filter never applied — cases contaminated the 'normal brain' column means

Formula (GBM-matched):
  for each (peak, class):
    evaluable_donors[class] = donors with >=100 cells of class  (same as strength denominator)
    mean_frac_fresh[peak, class] = mean(fraction_open_direct) across evaluable donors
                                   from PASS/PASS_adjacent pool only
  for each peak:
    top_class            = argmax(mean_frac_fresh)
    max_off_target_frac  = max(mean_frac_fresh) over classes != top_class
    specificity_ratio    = mean_frac_fresh[top] / (max_off_target_frac + eps)
    log2FC_fresh         = log2((top + eps) / (off + eps))
"""
import polars as pl
import numpy as np
from pathlib import Path

ROOT = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas 10-13-2026/atlas_v2.0.0")
DIRECT_DIR = ROOT / "h5ads_v2.2_nnz_rebuild"
SCORED = ROOT / "replicability_major_v2.2_nnz_rebuild_extended.parquet"
META = ROOT / "audit_2026-10-25/metadata_FINAL_v17.parquet"
OUT = ROOT / "replicability_major_v2.3_fresh_specificity.parquet"

MIN_CELLS = 100
EPS = 1e-6

meta = pl.read_parquet(META)
pass_donors = set(meta.filter(pl.col('atlas_usability_tier').is_in(['PASS','PASS_adjacent']))['donor_id'].to_list())
print(f"[1/3] PASS pool: {len(pass_donors)} donors", flush=True)

files = sorted(DIRECT_DIR.glob("*_direct_nnz1k.parquet"))
print(f"      cohort parquets: {len(files)}", flush=True)

# Pass 1: per-(peak, class) aggregation restricted to PASS donors
# Chunk each cohort by donor to keep peak RAM small on huge cohorts (hmba-bg has 204 donors × ~20M rows).
import tempfile, gc, shutil
tmp_dir = Path(tempfile.mkdtemp(prefix="nb_spec_v23_"))
print(f"      tmp dir: {tmp_dir}", flush=True)
pass_donors_list = list(pass_donors)

for i, f in enumerate(files):
    # Get donors in this file
    cohort_donors = (pl.scan_parquet(f).select('donor_id').unique()
                     .collect(engine='streaming')['donor_id'].to_list())
    cohort_pass_donors = [d for d in cohort_donors if d in pass_donors]
    print(f"  [{i+1}/{len(files)}] {f.name}  ({len(cohort_pass_donors)} PASS donors)", flush=True)
    if not cohort_pass_donors:
        continue

    # Chunk donors into batches of 50 to cap per-chunk RAM
    BATCH = 50
    for batch_start in range(0, len(cohort_pass_donors), BATCH):
        batch = cohort_pass_donors[batch_start:batch_start+BATCH]
        d = pl.scan_parquet(f)
        d = d.filter(pl.col('donor_id').is_in(batch))
        d = d.filter(pl.col('n_cells_total_direct') >= MIN_CELLS)
        d = d.with_columns(pl.col('peak_num').cast(pl.Int64))
        agg = d.group_by(['peak_num','major_class']).agg([
            pl.col('fraction_open_direct').sum().alias('sum_frac'),
            pl.col('donor_id').n_unique().alias('n_evaluable_fresh'),
        ]).collect(engine='streaming')
        if agg.shape[0] > 0:
            agg.write_parquet(tmp_dir / f"{i:03d}_{batch_start:04d}.parquet")
        del agg
        gc.collect()

# Final merge: streaming scan over all tmp parquets, group + sink
print(f"  merging {len(files)} tmp parquets lazily", flush=True)
final_lazy = pl.scan_parquet(str(tmp_dir / "*.parquet")).group_by(['peak_num','major_class']).agg([
    pl.col('sum_frac').sum(),
    pl.col('n_evaluable_fresh').sum(),
]).with_columns(
    (pl.col('sum_frac') / pl.col('n_evaluable_fresh').cast(pl.Float64)).alias('mean_frac_fresh')
)
# Streaming sink to another tmp parquet
final_pq = tmp_dir / "FINAL.parquet"
final_lazy.sink_parquet(final_pq)
final_agg = pl.read_parquet(final_pq)
print(f"[2/3] (peak, class) rows: {final_agg.shape[0]:,}", flush=True)
# Cleanup intermediate per-cohort tmp files but keep FINAL
import glob as _glob
for p in _glob.glob(str(tmp_dir / "*_*.parquet")):
    Path(p).unlink(missing_ok=True)

# Pass 2: per-peak top vs off-target
# Pivot to peak x class matrix of mean_frac_fresh
pivot = final_agg.pivot(
    values='mean_frac_fresh',
    index='peak_num',
    on='major_class',
).fill_null(0.0)

class_cols = [c for c in pivot.columns if c != 'peak_num']
print(f"      classes: {class_cols}", flush=True)

# Numpy for speed — compute PER-(peak, class) log2FC: log2(this_class / max(other_classes))
arr = pivot.select(class_cols).to_numpy()  # (n_peaks, n_classes)
top_idx = arr.argmax(axis=1)
top_val = arr.max(axis=1)
top_classes = np.array(class_cols)[top_idx]

# Per-(peak, class) off-target max: for each column j, max over all other columns
# Trick: find top-2 values per row; off_other = top if class != top_class, else second
# More general: for each class j, mask column j and take row-max
n_peaks, n_classes = arr.shape
log2FC_per_class = np.zeros_like(arr)
off_per_class = np.zeros_like(arr)
for j in range(n_classes):
    masked = arr.copy()
    masked[:, j] = -np.inf
    off = masked.max(axis=1)
    off = np.where(np.isinf(off), 0.0, off)
    off_per_class[:, j] = off
    log2FC_per_class[:, j] = np.log2((arr[:, j] + EPS) / (off + EPS))

# Melt to long form (peak, class, log2FC_fresh, max_off_target_frac_fresh, specificity_ratio_fresh)
peak_nums = pivot['peak_num'].to_numpy().repeat(n_classes)
class_long = np.tile(np.array(class_cols), n_peaks)
log2FC_long = log2FC_per_class.flatten()
off_long = off_per_class.flatten()
# this_class_mean_frac
this_frac_long = arr.flatten()
spec_df = pl.DataFrame({
    'peak_num': pl.Series(peak_nums, dtype=pl.Int64),
    'major_class': class_long,
    'mean_frac_top_fresh': this_frac_long,       # this class's mean_frac at this peak
    'max_off_target_frac_fresh': off_long,
    'specificity_ratio_fresh': (this_frac_long + EPS) / (off_long + EPS),
    'log2FC_fresh': log2FC_long,
})

# Also build the per-peak global top_class_fresh (constant per peak, derived from top class across classes)
per_peak_top = pl.DataFrame({
    'peak_num': pivot['peak_num'].cast(pl.Int64),
    'top_class_fresh': top_classes,
})

# Melt final_agg to long form keyed by (peak, class)
long = final_agg.select(['peak_num','major_class','mean_frac_fresh','n_evaluable_fresh']).with_columns(pl.col('peak_num').cast(pl.Int64))

# Join into scored parquet
print(f"[3/3] joining to scored parquet + overwriting stale columns", flush=True)
scored = pl.read_parquet(SCORED)
scored = scored.with_columns(
    pl.col('peak_id').str.strip_prefix('peak_').cast(pl.Int64).alias('peak_num')
)

out = scored.join(long, on=['peak_num','major_class'], how='left')
out = out.join(spec_df, on=['peak_num','major_class'], how='left')
out = out.join(per_peak_top, on='peak_num', how='left')
# Null-fill for peaks absent from direct quantification
out = out.with_columns([
    pl.col('mean_frac_fresh').fill_null(0.0),
    pl.col('n_evaluable_fresh').fill_null(0),
    pl.col('specificity_ratio_fresh').fill_null(1.0),
    pl.col('log2FC_fresh').fill_null(0.0),
])
# ALSO recompute strength_wn_gbm from mean_frac_fresh — same formula as nb_rescore_gbm_style.py
# This replays GBM-style scoring for new cell types that didn't exist when the Oct 25 rescore ran.
REF_WIDTH = 500
out = out.with_columns(
    pl.when(pl.col('peak_width') > 0)
      .then(pl.col('mean_frac_fresh') * (REF_WIDTH / pl.col('peak_width').cast(pl.Float64)))
      .otherwise(0.0)
      .alias('strength_wn_gbm_fresh')
)
# Overwrite the old strength_wn_gbm (which is correct for pre-existing classes but zero for new ones)
# Keep the old value where fresh is zero (= no data in pseudobulks); use fresh otherwise.
out = out.with_columns(
    pl.when(pl.col('strength_wn_gbm_fresh') > 0)
      .then(pl.col('strength_wn_gbm_fresh'))
      .otherwise(pl.col('strength_wn_gbm'))
      .alias('strength_wn_gbm')
)

# OVERWRITE the stale columns (keep old names for HTML compatibility)
# Old column: specificity_ratio -> write specificity_ratio_fresh into it
replace_map = {
    'specificity_ratio': 'specificity_ratio_fresh',
    'max_off_target_frac': 'max_off_target_frac_fresh',
    'mean_frac': 'mean_frac_fresh',
    'log2FC_top_vs_next': 'log2FC_fresh',
}
# Only replace columns that exist in scored
exprs = []
for old, new in replace_map.items():
    if old in out.columns and new in out.columns:
        exprs.append(pl.col(new).alias(old))
out = out.with_columns(exprs)

out.write_parquet(OUT)
print(f"  wrote {OUT.name} ({out.shape[0]:,} rows)", flush=True)
print(f"  cohort pool size: {len(pass_donors)} control donors (metadata_v7 tier=PASS/PASS_adjacent)")
