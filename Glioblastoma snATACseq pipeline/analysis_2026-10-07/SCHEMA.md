# Frozen schema — GBM Enhancer Atlas v1.0.0 (2026-10-08)

**This schema is frozen for deposit.** Zenodo deposits are versioned but not editable,
so column names and types below are the contract. Any future renames or additions
require a new Zenodo version (DOI increments on major schema change).

Column names use **snake_case**; no column will ever be renamed in-place
without a new deposit version. Dtypes are the exact Polars types shipped.

## `data/matrix/pan_malignant_scores_v4.parquet` (544,735 rows)

Per-peak summary across the 45-patient pan-malignant pool.

| column | dtype | meaning |
|---|---|---|
| `peak_id` | String | CATLAS peak id, format `chrN:start-end` (GRCh38) |
| `chrom` | String | GRCh38 chromosome (`chr1`-`chr22`, `chrX`) |
| `start` | Int64 | GRCh38 start (0-based) |
| `end` | Int64 | GRCh38 end (exclusive) |
| `strength_mean` | Float32 | Mean `frac_accessible` across patients where peak is detected |
| `strength_median` | Float32 | Median `frac_accessible` across patients |
| `n_patients_accessible` | UInt32 | Count of patients where ≥10% of their malignant cells access the peak |
| `n_patients_detected` | UInt32 | Count of patients where peak had ≥1 fragment in ≥1 malignant cell |
| `n_cohorts` | UInt32 | Count of the 8 cohorts with ≥1 patient detecting the peak |
| `consistency` | Float64 | Fraction of detecting patients where accessibility ≥ 10% |

## `data/matrix/pan_malignant_matrix_v4.parquet` (10.4M rows, 10.0 GB uncompressed)

Long-form (peak × patient × cohort) matrix for the pan-malignant pool.

| column | dtype | meaning |
|---|---|---|
| `peak_id` | String | CATLAS peak id |
| `chrom` | String | GRCh38 chromosome |
| `start` | Int64 | start |
| `end` | Int64 | end |
| `cohort` | String | Source cohort slug (see `cohort_attribution.csv`) |
| `patient_id` | String | Standardized patient id (unique within cohort) |
| `modality` | String | `snATAC` or `multiome` |
| `n_cells_malignant` | Int64 | Count of pan-malignant cells (CNV | label) for this patient |
| `n_accessible` | Int64 | Pan-malignant cells with ≥1 fragment at this peak |
| `frac_accessible` | Float32 | `n_accessible / n_cells_malignant` |
| `mean_counts` | Float32 | Mean fragment count per pan-malignant cell |

## `data/matrix/candidate_scores.parquet` (5.6M rows)

Per-peak × per-cell-type scores across the full 11-cell-type atlas.

| column | dtype | meaning |
|---|---|---|
| `peak_id` | String | CATLAS peak id |
| `cell_type` | String | One of: astrocyte, endothelial, GABA_neuron, malignant_unresolved, microglia, neuron, oligodendrocyte, OPC, T_cell, TAM, unassigned |
| `strength` | Float32 | Mean `frac_accessible` for this (peak, cell_type) across patients |
| `n_patients` | UInt32 | Count of patients with this cell type scoring this peak |
| `consistency` | Float64 | Fraction of detecting patients with accessibility ≥ 10% |
| `n_cohorts` | UInt32 | Cross-cohort replication count |
| `selectivity` | Float32 | `strength[this cell_type] / max(strength[other cell_types])` |
| `posterior_mean` | Float64 | Pass-2 hierarchical Bayesian posterior mean (top 2000 per cell type; null otherwise) |
| `posterior_ci_lo` | Float64 | Pass-2 posterior 95% CI lower (null for pass-1-only peaks) |
| `posterior_ci_hi` | Float64 | Pass-2 posterior 95% CI upper (null otherwise) |
| `posterior_p_specific` | Float64 | Pass-2 posterior P(peak is cell-type-specific) (null otherwise) |
| `pass2_was_run` | Boolean | True if hierarchical pass-2 was run on this (peak, cell_type) |
| `composite_score` | Float64 | Pass-1 closed-form composite = strength × selectivity × consistency × (n_cohorts/8) |

## `data/matrix/enhancer_candidate_matrix.parquet` (209M rows, 3.5 GB)

Long-form (peak × patient × cohort × cell_type) atlas workhorse. Schema in
`data/matrix/MATRIX_SCHEMA.md`. **Frozen.**

## `data/matrix/surprises.parquet` (1.1 GB)

Pattern-detector outputs. **Some columns are 100% null** (upstream artifact,
tracked as QC #15; see `qc/FINDINGS.md`). Use only populated columns: `peak_id`,
`detector`, `cell_type`, `pattern`.

## Report CSV schemas (shortlist deliverables)

### `reports/pan_malignant_top50_distal_selective_v3.csv` and `..._top100_all_v3.csv`

| column | meaning |
|---|---|
| `peak_id`, `chrom`, `start`, `end` | GRCh38 locus |
| `nearest_gene` | **Annotation only, not a filter.** Nearest TSS gene from GENCODE basic v45 |
| `dist_to_tss_signed` | Signed distance in bp; positive = downstream of TSS |
| `distance_category` | Categorical annotation (near 2-10 kb / distal 10-100 kb / far-distal 100-500 kb / gene-desert >500 kb). **Not a filter.** For AAV cloning, the enhancer is extracted from genomic context — native distance irrelevant. |
| `gbm_tags` | Gene-set membership (glioma driver, GSC/lineage TF, Neftel state, chromatin remodeler, TME myeloid, housekeeping) |
| `n_cohorts` | 0-8 |
| `n_patients_detected` | Count detecting out of 45 |
| `n_patients_accessible` | Count with ≥10% accessibility |
| `consistency` | 0.0-1.0 |
| `strength_mean` | 0.0-1.0 |
| `strength_max_nonmalignant` | 0.0-1.0 |
| `selectivity_vs_nonmal` | `strength_mean / strength_max_nonmalignant`, floor 0.01 on denominator |
| `nonmal_mean` | Mean strength across non-malignant cell types |
| `nonmal_sd` | SD strength across non-malignant cell types |
| `z_daigle` | `(strength_mean − nonmal_mean) / nonmal_sd` — Allen Institute Armamentarium specificity criterion |
| `passes_daigle_z2` | Boolean: `z_daigle >= 2.0` |
| `aav_score` | 0-1 weighted composite (consistency 35% + n_cohorts 25% + patient breadth 15% + selectivity 15% + strength 10%) |

### `reports/per_celltype/<cell_type>_top50.csv`

Same core columns, with `cell_type` identified by file name and `strength` + `selectivity`
in that cell type's reference frame. Includes `z_daigle` + `passes_daigle_z2`.

### `reports/pan_myeloid_top50_shared.csv` + `pan_myeloid_top50_distal_selective.csv`

Pan-myeloid union (TAM ∪ microglia) specific columns:

| column | meaning |
|---|---|
| `strength_TAM`, `strength_microglia` | Per-source strengths |
| `strength_myeloid` | `max(strength_TAM, strength_microglia)` |
| `strength_other_max` | Max strength across non-myeloid cell types |
| `selectivity_myeloid` | `strength_myeloid / strength_other_max` |
| `other_mean`, `other_sd` | For Daigle Z against non-myeloid |
| `z_daigle` | `(strength_myeloid − other_mean) / other_sd` |
| `passes_daigle_z2` | Boolean |
| `dominant_source` | Which of TAM/microglia has higher strength |
| `shared_strength_mean` | `(strength_TAM + strength_microglia) / 2` (SHARED table only) |
| `shared_rank_score` | SHARED-table ranking metric |

## BED files

GRCh38 BEDs for each shortlist:
- Column 1-3: `chrom`, `start`, `end`
- Column 4: `<nearest_gene>_<±kb>_<peak_id>`
- Column 5: Score × 1000 (BED convention: int 0-1000)
- Column 6: `.` (strand not meaningful)
