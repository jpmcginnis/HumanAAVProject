#!/usr/bin/env python3
"""v2.4: recompute n_donors + n_studies per (peak, class) with GBM evaluability gate.

Previously the HTML displayed n_donors / n_studies read from pre-nnz-filter columns,
which could overcount by including donors that failed the ≥100-cells evaluability gate
or whose cells were dropped by the nnz≥1k per-cell filter.

Fresh definition (matches strength denominator exactly):
  per (peak, class):
    n_donors_fresh  = # unique donors where
                        - donor is in metadata_v7 PASS/PASS_adjacent pool
                        - donor has ≥ MIN_CELLS cells of class  (evaluable)
                        - fraction_open_direct > 0            (peak has signal)
    n_studies_fresh = # unique study_id among those detected donors

Writes to replicability_major_v2.4_fresh_counts.parquet with n_donors and n_studies
columns overwritten.
"""
import polars as pl
from pathlib import Path

ROOT = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas 10-13-2026/atlas_v2.0.0")
DIRECT_DIR = ROOT / "h5ads_v2.2_nnz_rebuild"
IN  = ROOT / "replicability_major_v2.3_fresh_specificity.parquet"
OUT = ROOT / "replicability_major_v2.4_fresh_counts.parquet"
META = ROOT / "audit_2026-10-25/metadata_FINAL_v17.parquet"
MIN_CELLS = 100

meta = pl.read_parquet(META)
pass_donors = set(meta.filter(pl.col('atlas_usability_tier').is_in(['PASS','PASS_adjacent']))['donor_id'].to_list())
print(f"PASS pool: {len(pass_donors)} donors", flush=True)

files = sorted(DIRECT_DIR.glob("*_direct_nnz1k.parquet"))

# -- Pass 1: compute per-class n_evaluable (unique PASS donors with >=MIN_CELLS of class)
eval_parts = []
for f in files:
    d = pl.scan_parquet(f)
    d = d.filter(pl.col('donor_id').is_in(list(pass_donors)))
    d = d.filter(pl.col('n_cells_total_direct') >= MIN_CELLS)
    agg = d.with_columns(
        pl.concat_str([pl.col('donor_id'), pl.lit('|'), pl.col('study_id')]).alias('donor_study_key')
    ).group_by('major_class').agg(pl.col('donor_study_key').n_unique().alias('n_eval_cohort')).collect(engine='streaming')
    eval_parts.append(agg)
n_eval_class = pl.concat(eval_parts).group_by('major_class').agg(
    pl.col('n_eval_cohort').sum().alias('n_evaluable_class_fresh')
)
print(f"per-class evaluable PASS donors:\n{n_eval_class.sort('n_evaluable_class_fresh', descending=True)}")

# -- Pass 2: per-(peak, class) n_donors with signal + n_studies
# OOM-safe: chunk per cohort into donor batches, dump each batch's aggregate to tmp parquets,
# then streaming scan + group lazily for final.
import tempfile, gc, glob as _glob
from pathlib import Path as _Path
tmp_dir = _Path(tempfile.mkdtemp(prefix="nb_counts_v24_"))
print(f"  tmp dir: {tmp_dir}", flush=True)
pass_donors_list = list(pass_donors)
for i, f in enumerate(files):
    cohort_donors = (pl.scan_parquet(f).select('donor_id').unique()
                     .collect(engine='streaming')['donor_id'].to_list())
    cohort_pass_donors = [dn for dn in cohort_donors if dn in pass_donors]
    print(f"  [{i+1}/{len(files)}] {f.name} ({len(cohort_pass_donors)} PASS donors)", flush=True)
    if not cohort_pass_donors:
        continue
    BATCH = 50
    # Within one cohort, study_id is constant → skip n_unique on study_id (always 1 per cohort).
    # donor_study_key equals donor_id for same-cohort rows → use plain donor_id n_unique.
    # We also tag each batch with cohort_idx so we can n_unique(cohort_idx) at final merge for n_studies.
    for batch_start in range(0, len(cohort_pass_donors), BATCH):
        batch = cohort_pass_donors[batch_start:batch_start+BATCH]
        d = pl.scan_parquet(f)
        d = d.filter(pl.col('donor_id').is_in(batch))
        d = d.filter(pl.col('n_cells_total_direct') >= MIN_CELLS)
        d = d.with_columns(pl.col('peak_num').cast(pl.Int64))
        agg = d.group_by(['peak_num','major_class']).agg([
            pl.col('donor_id').n_unique().alias('n_donors_cohort'),
        ]).with_columns(pl.lit(i, dtype=pl.Int32).alias('cohort_idx')).collect(engine='streaming')
        if agg.shape[0] > 0:
            agg.write_parquet(tmp_dir / f"{i:03d}_{batch_start:04d}.parquet")
        del agg
        gc.collect()

final_pq = tmp_dir / "FINAL.parquet"
pl.scan_parquet(str(tmp_dir / "*_*.parquet")).group_by(['peak_num','major_class']).agg([
    pl.col('n_donors_cohort').sum().alias('n_donors_fresh'),
    pl.col('cohort_idx').n_unique().alias('n_studies_fresh'),
]).sink_parquet(final_pq)
final = pl.read_parquet(final_pq)
for p in _glob.glob(str(tmp_dir / "*_*.parquet")):
    _Path(p).unlink(missing_ok=True)
# Attach fresh evaluable denominator
final = final.join(n_eval_class, on='major_class', how='left')
# Fresh patient detection rate
final = final.with_columns(
    pl.when(pl.col('n_evaluable_class_fresh') > 0)
      .then(pl.col('n_donors_fresh').cast(pl.Float64) / pl.col('n_evaluable_class_fresh').cast(pl.Float64))
      .otherwise(0.0).alias('patient_detection_rate_fresh')
)
print(f"(peak, class) rows: {final.shape[0]:,}", flush=True)

scored = pl.read_parquet(IN)
if 'peak_num' not in scored.columns:
    scored = scored.with_columns(pl.col('peak_id').str.strip_prefix('peak_').cast(pl.Int64).alias('peak_num'))
scored = scored.with_columns(pl.col('peak_num').cast(pl.Int64))
out = scored.join(final, on=['peak_num','major_class'], how='left')
# Also fill n_evaluable_class_fresh per-class for peaks absent from pseudobulks
out = out.join(n_eval_class, on='major_class', how='left', suffix='_cls')
# Prefer the per-peak value if present, else per-class value
if 'n_evaluable_class_fresh_cls' in out.columns:
    out = out.with_columns(
        pl.coalesce(['n_evaluable_class_fresh','n_evaluable_class_fresh_cls']).alias('n_evaluable_class_fresh')
    ).drop('n_evaluable_class_fresh_cls')

out = out.with_columns([
    pl.col('n_donors_fresh').fill_null(0),
    pl.col('n_studies_fresh').fill_null(0),
    pl.col('n_evaluable_class_fresh').fill_null(0),
    pl.col('patient_detection_rate_fresh').fill_null(0.0),
])

# OVERWRITE the stale columns the HTML reads from
out = out.with_columns([
    pl.col('n_donors_fresh').alias('n_donors'),
    pl.col('n_studies_fresh').alias('n_studies'),
    pl.col('n_donors_fresh').alias('n_detected'),                    # events-count → unique-donors-count
    pl.col('n_evaluable_class_fresh').alias('n_evaluable_class'),    # dedupe + PASS filter
    pl.col('patient_detection_rate_fresh').alias('patient_detection_rate_gbm'),
])

out.write_parquet(OUT)
print(f"wrote {OUT.name} ({out.shape[0]:,} rows)", flush=True)

# Sanity
sanity = (out.filter(pl.col('major_class')=='Microglia')
             .sort('strength_wn_gbm', descending=True)
             .head(5)
             .select(['peak_id','n_donors','n_studies','n_detected','n_evaluable_class','patient_detection_rate_gbm']))
print(sanity)
