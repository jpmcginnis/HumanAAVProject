# Normal brain snATAC-seq enhancer atlas — pipeline

Companion to the [GBM snATAC-seq pipeline](../Glioblastoma%20snATACseq%20pipeline/) in the same repository. Builds the per-cell-type enhancer-candidate atlas served live at [data.mcginnislab.org/normal-brain-atlas/top-1000](https://data.mcginnislab.org/normal-brain-atlas/top-1000). Downloadable per-sample metadata at [`/normal-brain-atlas/metadata.csv`](https://data.mcginnislab.org/normal-brain-atlas/metadata.csv).

## Scope

Every publicly-accessible human brain snATAC-seq / snMultiome deposit — healthy, developmental, surgical margin, and disease (cases are excluded from scoring but are tracked in the metadata). Covers hippocampus, cortex, basal ganglia, cerebellum, brainstem, spinal cord. **33 cohorts, ~953 unique donors, 713 atlas-eligible (PASS + PASS-adjacent)**.

See [`docs/COHORTS.md`](docs/COHORTS.md) for the per-cohort ingest route.

## Pipeline overview

```
PER-COHORT INGEST (various routes per deposit — see docs/COHORTS.md)
       │  one *_direct_nnz1k.parquet per cohort, columns:
       │  peak_num, major_class, donor_id, n_cells_access_direct,
       │  n_cells_total_direct, fraction_open_direct, study_id
       ▼
ATLAS SCORING (src/atlas/)
   nb_add_new_classes.py          # extend base scored parquet with new cell-type classes
   nb_specificity_recompute.py    # per-(peak, class) log2FC vs best off-target class
   nb_counts_recompute.py         # per-(peak, class) n_donors with signal, n_studies, patient_detection_rate
   fix_strength_to_gbm.py         # strength_mean_evaluable_wn = sum_frac / n_evaluable_class × (500/peak_width)
   build_nb_atlas_html.py         # single-file self-contained HTML atlas with DataTables + row-expand UCSC links
       ▼
METADATA (src/metadata/)
   build_cohort_urls.py           # per-cohort data URL + PMID/DOI backfill (NCBI eutils + hand-curated)
       ▼
DEPLOY (manual: cp HTML + CSV → mcginnislab-data repo, git push)
```

## Peak set

Scored against the 990K-peak "v1.0.0 pre-Design-B union peak set" (median width 895 bp, mean 2,434 bp). This set includes wide peaks that absorb neighboring signal. For AAV cloning, the row expansion in the HTML surfaces the narrow Design B sub-peaks (2.06M peaks, median 518 bp) that fall inside each wide peak — those are the actual cloneable fragments. See `narrow_lookup.py` for the mapping.

## Scoring details

- **Quantification**: binary accessibility per cell (fraction_open_direct = cells accessing peak / total cells of that class per donor).
- **QC floor**: `nnz ≥ 1000` union peaks per cell. Applied at ingest. Required because some deposits (4DN hippocampus, GSE322785) shipped permissive barcode calling with median nnz 3–4.
- **Evaluability gate**: ≥ 100 cells of that class per donor.
- **Strength (GBM-matched)**: `strength_mean_evaluable_wn = (sum_frac / n_evaluable_class) × (500 / peak_width)`. Donors evaluable for the class but without peak signal count as 0 in the mean. Same convention as the GBM atlas.
- **Specificity**: `log2FC_fresh = log2(mean_frac_this_class / max_mean_frac_other_classes)` with epsilon handling.
- **Primary ranking**: `primary_score = (spec_pct + strength_pct) / 2`, where each is the percentile rank of its component **within distal peaks (|dist to TSS| ≥ 2 kb) of that cell type**. Promoter-proximal peaks are filtered out before top-N selection because they aren't useful AAV enhancer scaffolds.
- **Known limitation**: `primary_score` biases toward rare+specific peaks. For prevalence-filtered views, filter rows on `patient_detection_rate_gbm ≥ 0.5`. See `docs/PIPELINE.md`.

## Running the atlas-build chain

```bash
cd src/atlas
python3 nb_add_new_classes.py           # extend base parquet with new classes if pseudobulks have any
python3 nb_specificity_recompute.py     # ~20 min; writes replicability_major_v2.3_fresh_specificity.parquet
python3 nb_counts_recompute.py          # ~20 min; writes replicability_major_v2.4_fresh_counts.parquet
python3 fix_strength_to_gbm.py          # ~2 min; writes replicability_major_v2.5_gbm_strength.parquet
python3 build_nb_atlas_html.py          # ~5 min; writes top-1000.html (self-contained, ~25 MB)
```

All file paths at the top of each script.

## Testing

```bash
cd tests
bash smoke_test.sh
```

Runs the full atlas-build chain against a 2-donor synthetic cohort; verifies the output HTML has the expected top-1 peak + a non-empty specificity bar for the dominant class. ~1 min.

## Requirements

See [`requirements.txt`](requirements.txt). Key dependencies: `polars`, `pandas`, `numpy`, `scipy`, `pyranges`, `h5py` (snapatac2 cohort ingest only).

## Attribution

See [`docs/COHORTS.md`](docs/COHORTS.md) for per-cohort citation. Live CSV with every donor's deposit URL + PMID/DOI: [data.mcginnislab.org/normal-brain-atlas/metadata.csv](https://data.mcginnislab.org/normal-brain-atlas/metadata.csv).
