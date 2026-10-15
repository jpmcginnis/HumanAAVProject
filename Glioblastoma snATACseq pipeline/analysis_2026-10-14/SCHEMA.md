# Frozen schema — Harmonized human glioblastoma chromatin accessibility atlas v2026-10-14 (2026-10-14)

**This schema is frozen for deposit.** Zenodo deposits are versioned but not editable,
so column names and types below are the contract. Any future renames or additions
require a new Zenodo version (DOI increments on major schema change).

Column names use **snake_case**; no column will ever be renamed in-place without a
new deposit version.

Genome build: **GRCh38 (hg38)** throughout. All intervals are 0-based half-open BED convention.

---

## `data/matrix/union_peak_reference.bed` (1,958,753 rows)

4-column BED of the union peak set across 8 GBM cohorts.

| col | name | meaning |
|---|---|---|
| 1 | `chrom` | hg38 chromosome (`chr1`–`chr22`, `chrX`, `chrY`, `chrM`) |
| 2 | `start` | 0-based start |
| 3 | `end` | exclusive end |
| 4 | `sources` | comma-separated list of contributing cohort slugs (e.g. `greenwald_2024,catlas,sundaram_2024`) |

Mean width 594 bp, max 6.7 kb. No transitive-chaining artifacts (merges are at single-cohort peak-call precision).

## `data/matrix/peak_provenance.parquet` (long-form, one row per `(union_peak, source)`)

| column | dtype | meaning |
|---|---|---|
| `peak_id` | String | `chrN:start-end` of the union peak |
| `chrom` | String | hg38 chromosome |
| `start` | Int64 | union peak start |
| `end` | Int64 | union peak end |
| `source` | String | contributing cohort slug (one of 8 cohorts or `catlas`) |
| `source_chrom` | String | the source peak's chromosome |
| `source_start` | Int64 | the source peak's start |
| `source_end` | Int64 | the source peak's end |
| `n_sources_total` | Int32 | total contributing sources for this union peak |

## `data/matrix/enhancer_candidate_matrix.parquet` (1,958,753 × 338)

Wide-form patient × cell-type pseudobulk matrix. 1 row per union peak, 337 pseudobulk columns.

| column | dtype | meaning |
|---|---|---|
| `peak_id` | String | `chrN:start-end` matching union_peak_reference.bed |
| `<cohort>__<patient_id>__<cell_type>` | Int64 (raw counts) | integer fragment counts for that (patient × cell_type) in that cohort |

**Cohort coverage** (337 total pseudobulk columns):
- `greenwald_2024`: 144 columns (12 patients × up to 14 coarse labels, incl. 4 Neftel states + `malignant`)
- `sundaram_2024`: 62 columns
- `zhao_2024`: 50 columns
- `sojka_2025`: 33 columns
- `mathewson_2021`: 20 columns
- `guo_2020`: 15 columns
- `guilhamon_2021`: 10 columns
- `raviram_2023`: 3 columns

**Pool rule.** A cell with `malignant_cnv == 1` OR `cell_type == "malignant_unresolved"` is bucketed as `malignant`, not as its gbmap cell_type. This fixes "hidden malignant" contamination of non-malignant pseudobulks seen in earlier versions.

**MIN_CELLS = 50** per (patient × label) group — tinier groups skipped as noisy.

## `data/matrix/greenwald_fine_pseudobulk_backup.parquet`

Preserved fine-grained Greenwald pseudobulk (99 labels × 12 patients) using the full annotation granularity from Greenwald et al. (AC-progenitor-like-1-4, OPC-NPC-like-1-2, pro-inflammatory TAMs, astrocyte-like TAMs, CD8+ cytotoxic T cells, Exc L4-5 IT, etc.). Column format: `greenwald_2024__<patient>__<fine_label>`. Preserved so future releases can re-remap and generate super-fine reports without re-reading the 10.5B-nnz Greenwald ATAC h5ad.

## `data/matrix/slug_canonical.tsv`

Two-column TSV mapping legacy cohort slugs → canonical first-author_year slugs. Notably: `gse276177_khan_astro` → `sojka_2025` (misattribution corrected in v2026-10-11; the Sojka et al. 2025 Nat Cell Biol paper is the correct authorship).

---

## `reports/detail/<report>.parquet` (per-report score tables)

One parquet per report (16 reports total: 10 non-malignant per-celltype + 4 Neftel-state per-celltype + `pan_malignant` + `pan_myeloid`). Each row = one union peak scored against that report's target.

| column | dtype | meaning |
|---|---|---|
| `peak_id` | String | `chrN:start-end` |
| `chrom`, `start`, `end` | String, Int64, Int64 | hg38 locus |
| `target_mean_cpm` | Float64 | mean CPM across target patients (patient-level pseudobulks, patient-averaged) |
| `bg_mean_cpm` | Float64 | mean CPM across non-target columns genome-wide |
| `daigle_z_raw` | Float64 | `(target_mean − μ_bg) / σ_bg` across non-target cols, genome-wide null |
| `daigle_z_arm_adjusted` | Float64 | **primary ranking column.** Same formula but `μ_arm, σ_arm` are over non-target cols on the same chromosome arm as the peak. **Standardized per chromosome arm — not directly comparable across arms or to `daigle_z_raw`.** |
| `n_target_patients` | Int32 | unique (cohort, patient_id) tuples with nonzero signal in at least one target-label column. For single-label reports this equals the number of pseudobulk columns with signal; for `pan_malignant` / `pan_myeloid` where a patient can appear under multiple labels, it collapses them to one. |
| `n_target_pseudobulks` | Int32 | distinct patient × label pseudobulk columns with nonzero signal. Only differs from `n_target_patients` for the two pooled reports (`pan_malignant`, `pan_myeloid`), where it is strictly ≥ `n_target_patients`. |
| `n_target_cohorts` | Int32 | distinct cohorts those target patients span |
| `cohort_balance` | Float64 | fraction of target signal from the top-contributing cohort (1.0 = single-cohort peak) |
| `greenwald_fraction` | Float64 | fraction of total target CPM contributed by greenwald_2024 columns. **For `pan_malignant` the median is 0.78** (Greenwald is the only cohort with Neftel-split malignant labels, so it provides 60 of ~88 target columns); for `pan_myeloid` the median is 0.25 (cross-cohort balanced). Filter on this column to isolate strictly cross-cohort-replicating peaks. |
| `offtarget_max_cpm` | Float64 | max mean CPM across any non-target cell type |
| `offtarget_top_cell_type` | String | which non-target cell type is highest |
| `offtarget_pooled_mean` | Float64 | mean across non-target cell-type means |
| `n_offtarget_ct_above_1` | Int32 | count of non-target cell types with mean CPM > 1 |
| `aav_safety_ratio` | Float64 | `target_mean_cpm / offtarget_max_cpm` (clipped at a large value when offtarg_max ≈ 0) |
| `chromosome_arm` | String | hg38 centromere-split arm (`1p`, `1q`, …, `Xq`) |
| `chromosome_cnv_class` | String | `chr7_gain`, `chr10_loss`, or `none` — known-malignant arm-level CNVs |
| `patient_fraction_chr7_elevated` | Float64 (nullable) | malignant targets only: fraction of target patients where this peak falls on an arm with per-patient elevated signal (surfaces CNV-driven apparent specificity). Null for non-malignant reports. |

**Filters applied before writing:** `MIN_PATIENTS = 2` (n_target_patients ≥ 2). Ranked by `daigle_z_arm_adjusted` descending.

## `reports/detail/<report>_annotated.csv`

Same columns as the parquet, plus gene annotations:

| column | meaning |
|---|---|
| `nearest_gene` | RefGene symbol of nearest TSS |
| `dist_to_tss_signed` | signed distance (bp) from peak center to nearest TSS; negative = upstream |
| `distance_category` | `promoter` (≤2 kb), `proximal` (2–10 kb), `distal` (10–100 kb), `far-distal` (>100 kb) |

## `reports/detail/<report>_top100.csv` and `<report>_top50_distal.csv`

Ranked slice of the annotated CSV. `top50_distal` additionally filters `distance_category ≠ promoter`.

## `reports/detail/<report>_top30_distal.bed`

BED6 of the top-30 distal peaks for easy UCSC / IGV loading.

| col | content |
|---|---|
| 1-3 | `chrom`, `start`, `end` |
| 4 | name: `<nearest_gene>_<±kb>kb_<peak_id>` |
| 5 | score × 1000 (BED convention 0–1000) |
| 6 | strand placeholder `.` |

## `reports/<report>_top100.csv` (12 files, top-level polished)

Same schema as `reports/detail/<report>_top100.csv` but filtered with the **polished** criteria: `distance_category ≠ promoter`, `target_mean_cpm ≥ 1.0`, `n_target_patients ≥ 2`, then top-100 by `daigle_z_arm_adjusted`.

Only the 12 headline reports get a top-level polished CSV: 10 non-malignant per-celltype + `pan_malignant` + `pan_myeloid`. Neftel-state malignant subtypes are intentionally excluded from the top-level — they remain available in `reports/detail/per_celltype/`.

## `reports/top100_distal_combined.{csv,parquet,xlsx,html}` (1200 rows)

The 12 headline reports concatenated with a `report` and `rank` column prepended. The `.html` is interactive (sortable headers, text/select/number filters, hover tooltips on every column header + filter label explaining the equation and purpose).

## `reports/top10_per_arm_combined.{csv,xlsx}`

For each report × chromosome_arm pair, the top-10 peaks by `daigle_z_arm_adjusted`. CNV-diversity view — surfaces strong candidates on every arm, not just the globally top-ranked ones which may concentrate on a few arms.

## `reports/detail/top50_combined.{csv,xlsx,parquet}` and `top50_distal_combined.{...}` (800 rows)

The 16 reports (10 non-malignant + 4 Neftel + 2 pan) × 50 top rows. `top50_distal_combined` additionally filters `distance_category ≠ promoter`.

---

## Scored cohorts (8)

All 8 published GBM scATAC/snMultiome cohorts contribute to both the peak reference AND scoring:

| canonical slug | citation | accession | platform | n_patients | n_malignant_cells |
|---|---|---|---|---|---|
| `greenwald_2024` | Greenwald et al. 2024 Cell | E-MTAB-17183 | 10x snATAC + snMultiome | 12 | 438,994 |
| `sojka_2025` | Sojka et al. 2025 Nat Cell Biol | GSE276177 | 10x snATAC | 3 | 73,413 |
| `sundaram_2024` | Sundaram et al. 2024 Cell Rep | GDC/dbGaP phs001287 | 10x snATAC | 9 | 53,632 |
| `zhao_2024` | Zhao et al. 2024 Cell Rep Med | HRA004942 | snATAC | 7 | 591 |
| `mathewson_2021` | Mathewson 2021 Cell + Wang 2024 Sci Adv | GSE163654 + DNBelab | 10x + DNBelab | 5 | 2,202 |
| `guo_2020` | Guo et al. 2020 Nat Commun | GSE138794 | 10x snATAC | 3 | 1,817 |
| `guilhamon_2021` | Guilhamon et al. 2021 eLife | GSE139136 | 10x snATAC | 4 | 2,332 |
| `raviram_2023` | Raviram et al. 2023 PNAS 120:e2210991120 | GSE165037 | sci-ATAC | 2 | 2,990 |

See `cohort_attribution.csv` for data-use conditions and `CHANGELOG.md` for the slug rename rationale.
