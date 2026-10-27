#!/usr/bin/env python3
"""Append peak × class rows to the base scored parquet for cell types that only
exist in new pseudobulks (e.g. D1-MSN, D2-MSN, Granule, Purkinje from Menon/Corces).

Workflow: for each new class, duplicate the per-peak metadata rows from an existing
class (e.g. Excitatory) with scoring columns zeroed. Downstream recompute scripts
(nb_specificity_recompute_v23 + nb_counts_recompute_v24) will then fill those with
fresh values from the pseudobulks.
"""
import polars as pl
from pathlib import Path

ROOT = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas 10-13-2026/atlas_v2.0.0")
BASE = ROOT / "replicability_major_v2.2_nnz_rebuild.parquet"
OUT  = ROOT / "replicability_major_v2.2_nnz_rebuild_extended.parquet"
DIRECT_DIR = ROOT / "h5ads_v2.2_nnz_rebuild"

# Discover classes in pseudobulks
pb_classes = set()
for f in sorted(DIRECT_DIR.glob("*_direct_nnz1k.parquet")):
    c = pl.scan_parquet(f).select('major_class').unique().collect(engine='streaming')['major_class'].to_list()
    pb_classes.update(c)
print(f"classes across all pseudobulks: {sorted(pb_classes)}")

base = pl.read_parquet(BASE)
base_classes = set(base['major_class'].unique().to_list())
print(f"\nclasses in base scored parquet: {sorted(base_classes)}")

new_classes = sorted(pb_classes - base_classes)
print(f"\nnew classes to add: {new_classes}")
if not new_classes:
    print("nothing to do"); import sys; sys.exit(0)

# Template rows: take one existing class, duplicate per peak, swap class name
template_class = 'Excitatory'  # has rows for most peaks
template = base.filter(pl.col('major_class') == template_class)
print(f"\ntemplate class {template_class}: {template.shape[0]:,} rows")

# Identify scoring columns (to zero) vs peak-metadata columns (to keep)
peak_cols = {'peak_id','peak_num','peak_width','chrom','start','end','peak_type',
             'nearest_gene','gene_type','distance_to_tss'}
scoring_cols = [c for c in base.columns if c not in peak_cols and c != 'major_class']
print(f"peak-metadata cols kept: {sorted(c for c in template.columns if c in peak_cols)}")
print(f"scoring cols to zero/null: {len(scoring_cols)} cols")

new_frames = []
for cls in new_classes:
    df = template.clone()
    df = df.with_columns(pl.lit(cls).alias('major_class'))
    # Zero out scoring columns
    for c in scoring_cols:
        if c in df.columns:
            dtype = df[c].dtype
            if dtype in (pl.Float64, pl.Float32):
                df = df.with_columns(pl.lit(0.0).cast(dtype).alias(c))
            elif dtype in (pl.Int32, pl.Int64, pl.UInt32, pl.UInt64):
                df = df.with_columns(pl.lit(0).cast(dtype).alias(c))
            elif dtype == pl.Boolean:
                df = df.with_columns(pl.lit(False).alias(c))
            else:
                df = df.with_columns(pl.lit(None).cast(dtype).alias(c))
    new_frames.append(df)
    print(f"  prepared {cls}: {df.shape[0]:,} rows")

extended = pl.concat([base] + new_frames)
print(f"\nextended total: {extended.shape[0]:,} rows (base {base.shape[0]:,} + new {sum(f.shape[0] for f in new_frames):,})")
extended.write_parquet(OUT)
print(f"wrote {OUT.name}")
