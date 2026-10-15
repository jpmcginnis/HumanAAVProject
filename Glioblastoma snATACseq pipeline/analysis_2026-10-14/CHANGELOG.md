# CHANGELOG — Harmonized human glioblastoma chromatin accessibility atlas (snATAC-seq)

Semantic versioning: major = schema change, minor = new data / new reports, patch = fixes.

Everything numbered below v2026-10-14 was a working iteration — the artifact set at the root of this folder is **v2026-10-14 (current)**, with no version suffixes in filenames.

## v2026-10-14 — 2026-10-14 (CURRENT, released) — Raviram-attribution patch

Non-breaking metadata correction. No change to sample set, scoring, or scoring metrics.

- **GSE165037 attribution corrected.** The cohort previously labeled `nagaraja_2021` (and attributed to "Nagaraja et al. 2021 / Michelle Monje") is actually from the Raviram, Chen, Ren group at UCSD and published as Raviram et al. 2023 PNAS 120:e2210991120 ("Integrated analysis of single-cell chromatin state and transcriptome identified common vulnerability despite glioblastoma heterogeneity", NFIA/NFIB finding). Verified directly against the GEO series matrix: series title "Single Cell Analysis of Chromatin Accessibility Reveals Genetic and Regulatory Heterogeneity in Glioblastomas", contributors Raviram / Chen / Ren at UCSD, sci-ATAC combinatorial-indexing platform, 6 snATAC samples (3 IDH-WT GBMs: GBM4 / GBM9 / GBM12; 2 IDH-mutant astrocytomas: GBM1 / GBM11; non-tumor brain).
- Cohort slug renamed: **`nagaraja_2021` → `raviram_2023`** throughout the deposit:
  - `cohort_attribution.csv`: full row rewritten with corrected citation, PI (Bing Ren, UCSD), and expanded notes.
  - `data/matrix/enhancer_candidate_matrix.parquet`: 3 pseudobulk column names renamed (`nagaraja_2021__GBM4__malignant` → `raviram_2023__GBM4__malignant`, etc.).
  - `data/matrix/peak_provenance.parquet`: 78,432 rows' `source` column rewritten from `nagaraja_2021` to `raviram_2023`.
  - `data/matrix/union_peak_reference.bed`: 78,432 line-level `sources` column replacements.
  - `data/matrix/slug_canonical.tsv`: canonical slug updated; legacy `nagaraja_2021` mapping preserved as backwards-compat entry so downstream code using the old slug still resolves.
  - `SCHEMA.md`, `README.md`, `CHANGELOG.md`, `zenodo.json`, `scripts/pileup_region.py`, `tests/test_schema.py`, `workflow/REPLICATE_FOR_NORMAL_CORTEX.md` all updated.
- **Patient-count audit clarified.** Deposit records contain 50 patients across 8 cohorts (12 greenwald + 3 sojka + 9 sundaram + 11 zhao + 5 mathewson + 3 guo + 4 guilhamon + 3 raviram); 2 were excluded by QC — sundaram X008 (snapatac2 Rust panic) and raviram GBM9 (MIN_CELLS=50 floor). The headline 48-patient-pseudobulked total is unchanged.
- **`scripts/zenodo_upload.py` reliability improvement.** Retries extended from 3 → 10 attempts, with explicit handling of `ReadTimeout` alongside `ConnectionError` / `ChunkedEncodingError`. Backoff tuned to a 30s–5min window (vs the old exponential curve that could push a single file to 2.5 hours). Observed-helpful on multi-hundred-MB uploads over flaky links.
- **No change to atlas content, metrics, or scoring.** The 1.96M union peaks, 337 pseudobulk columns, Daigle-Z-arm-adjusted ranking, per-cell-type reports, and 12 polished top-100 deliverables are byte-identical except where the string `nagaraja_2021` has been substituted with `raviram_2023`. The 17/17 pytest regression suite passes.

## v2026-10-11-rc2 — 2026-10-11 (reviewer-driven pre-push revisions)

Diff vs. earlier internally-circulated v2026-10-11-rc1 (same date). If you were handed
a v2026-10-11 draft before this entry existed, you have rc1 and the numbers below are
what changed. If you cloned the GitHub repo, you have rc2 (now stamped v2026-10-11).

- **`n_target_patients` semantic fix.** Previously counted pseudobulk COLUMNS
  (patient × label groups) with signal, so `pan_myeloid` reported values up to
  55 against a 48-patient atlas. Fixed in-place to count unique
  `(cohort, patient_id)` tuples. The old quantity preserved as a new column
  `n_target_pseudobulks` (strictly ≥ `n_target_patients` for the two pooled
  reports; equal for single-label reports). Patch applied to all 16 report
  parquets with a vectorized full-matrix recomputation.
- **Patient-count correction.** Headline "51 patients" was wrong across README,
  CHANGELOG header, and zenodo.json; cohort_attribution.csv had `zhao_2024`
  transcribed as 7 patients when the h5ad holds 11 unique donor IDs — a
  4-patient attribution error. Corrected everywhere; `n_patients_included`
  per cohort now ties out to 48 (= 49 deposit records − 1 for X008 exclusion).
  See §"Patient-count provenance (disclosure)" below.
- **`greenwald_fraction` column added** to every report parquet (CPM-normalized
  share of target signal from greenwald_2024). Median for `pan_malignant` is
  0.78 (Greenwald has 5 Neftel labels × 12 patients = 60 of 88 target cols);
  for `pan_myeloid` 0.25 (cross-cohort balanced). Filterable in the HTML.
- **SCHEMA.md rewrite.** Previous SCHEMA described v3/v4 files that no longer
  exist (`pan_malignant_scores_v4.parquet`, `candidate_scores.parquet`,
  `surprises.parquet`). Now a frozen column dictionary for every artifact
  shipped in this release, with the `n_target_patients` vs `n_target_pseudobulks`
  distinction documented inline.
- **cohort_attribution.csv new columns:** `raw_fragments_released`
  (yes / no (h5ad only) — 2 of 8 cohorts shipped raw fragments publicly),
  `accession_version` (specific submission ids + which processed file was
  ingested), `atlas_build_date` (`2026-10-11`).
- **README gains a Known-Limitations entry** for (a) `pan_malignant` being
  Greenwald-weighted and (b) the TAM top-5 chrX pile-up driven by Greenwald's
  snMultiome label-transfer purity advantage over the other 6 label-transferred
  cohorts. Also an explicit Methods-style statement that X008 is excluded.
- **Arm-adjusted Z scale note.** README scoring-metrics table clarifies that
  `daigle_z_arm_adjusted` is standardized per chromosome arm and not directly
  comparable across arms or to `daigle_z_raw`.
- **Small-arm null stability check** written up: top-100 of every report has at
  most 3 peaks on chr21q and ≤2 on chr22q, zero on chr15p / chr21p / chr22p /
  chrY — so the per-arm σ estimates are not driving ranking. Documented as a
  Known-Limitations entry (confirmed non-issue).
- **+2 regression tests:** `test_n_target_patients_within_atlas_bounds` guards
  the arithmetic invariant that broke in rc1; `test_greenwald_fraction_in_range`
  guards the new column's `[0, 1]` bound and the known-high `pan_malignant` median.
  17/17 total passing — see `REGEN_LOG.md` at repo root.

## v2026-10-11 — 2026-10-11 (released)

**Headline:** First release of the atlas as a usable deliverable.
**48 patients** across 8 GBM cohorts (unique `(cohort, patient_id)` pseudobulks
with ≥50 cells in at least one label), 1,958,753 union peaks (74% de-novo vs
CATLAS), scored per cell type with patient-level replication, cohort batch
sanity, per-non-target-cell-type contamination for AAV safety, and arm-adjusted
Daigle Z for CNV correction.

### Patient-count provenance (disclosure)

Earlier working drafts of this release stated **"51 patients"** across the atlas.
That number was wrong on multiple fronts, caught by reviewer audit against the
pseudobulk matrix. Corrected in v2026-10-11:

- **`zhao_2024` was transcribed as 7 patients in prior `cohort_attribution.csv`.**
  The h5ad actually contains **11** unique donor IDs (P62, P64, P77, P78, P79,
  P80, P83, P84, P92, P98, P101). This is a data-provenance error, not a
  copyedit — the attribution file was wrong by 4 patients in one row. Corrected
  in-place; the matrix counts (11 columns per Zhao label) were right all along.
- **`sundaram_2024` X008 skip** was reflected in the matrix (8/9 patients
  pseudobulked) but `n_patients_included` in `cohort_attribution.csv` still
  showed 9. Corrected to 8; the exclusion rationale (snapatac2 Rust panic on
  this sample's obs structure) is now in the file's `notes` column, not only
  in README's Limitations section.
- **The old `n_target_patients` scored column was counting pseudobulk COLUMNS,
  not unique patients** — for `pan_myeloid` where target = TAM ∪ microglia, a
  patient contributing both labels was counted twice. Max was 55 (36 TAM cols
  + 19 microglia cols). Fixed in-place: the field now counts unique
  `(cohort, patient_id)` tuples with signal (max 38 for pan_myeloid, 40 for
  pan_malignant); the old quantity moved to `n_target_pseudobulks`.
  A regression test (`test_n_target_patients_within_atlas_bounds`) now guards
  the invariant `n_target_patients ≤ atlas_total ≤ n_target_pseudobulks`.
- Reconciled totals (post-fix): deposit records contain **50** patients
  (12 greenwald + 3 sojka + 9 sundaram + 11 zhao + 5 mathewson + 3 guo +
  4 guilhamon + 3 raviram). Of these, sundaram X008 was excluded (snapatac2 panic)
  and raviram GBM9 was excluded (MIN_CELLS=50 floor), leaving **48** patients
  pseudobulked. This is the denominator used
  everywhere in v2026-10-11 (README, SCHEMA, zenodo.json, tooltips, tests,
  cohort_attribution.csv `n_patients_included` column, and the HTML
  denominator-table tooltip on n_target_patients / n_target_cohorts).

### Peak reference

- `data/matrix/union_peak_reference.bed` — 1,958,753 peaks (mean width 594 bp,
  max 6.7 kb; no transitive-chaining artifacts).
- Composition: catlas 25.8%, greenwald 19.9%, sundaram 45.8%, sojka 36.5%,
  guo 16.7%, zhao 9.5%, mathewson 8.4%, guilhamon 9.5%, raviram 4.0%.
- **CATLAS-only (would be lost in a CATLAS-only pipeline): 132,960 (6.8%).**
- **De-novo (not in CATLAS): 1,453,136 (74.2%).**

### Peak-level provenance

- `data/matrix/peak_provenance.parquet` — one row per `(union_peak, source)`.
  Any merged interval traceable back to contributing cohort coordinates.

### Pseudobulk matrix (patient-level, with proper pool rule)

- `data/matrix/enhancer_candidate_matrix.parquet` — 1.96M peaks × 338 columns (337 `<cohort>__<patient_id>__<label>` pseudobulks + `peak_id`).
- Columns are `<cohort>__<patient_id>__<label>` so patient is the primary unit
  of replication and cohort is recoverable for batch sanity.
- Pool rule: a cell with `malignant_cnv == 1` OR `cell_type == "malignant_unresolved"`
  is bucketed as `malignant`, not as its gbmap cell_type. Fixes the "hidden
  malignant" contamination of astrocyte/OPC/neuron pseudobulks that affected
  earlier versions.
- MIN_CELLS = 50 per (patient × label) group — tiny groups skipped as noisy.
- **Scored cohorts (8, all GBM snATAC/snMultiome sources):** sundaram_2024,
  sojka_2025, zhao_2024, mathewson_2021, guo_2020, guilhamon_2021, raviram_2023,
  and **greenwald_2024** (10.5B-nnz CSR scored via an 8-worker parallel chunked
  decompress across the ATAC h5ad; labels transferred from the shared-barcode
  Multiome snRNA and remapped to the coarse schema, with 4 Neftel states
  — `malignant_AC` / `malignant_MES` / `malignant_OPC` / `malignant_NPC` —
  retained as separate detail-only reports; the fine-grained 99-label
  pseudobulk is preserved at `data/matrix/greenwald_fine_pseudobulk_backup.parquet`
  for future re-remapping).

### Cohort slugs — canonical first-author naming

Legacy slugs replaced site-wide:

| Legacy | v1 canonical | First author |
|---|---|---|
| gbm_space | greenwald_2024 | Greenwald et al. Cell |
| gse276177_khan_astro | **sojka_2025** | Sojka et al. Nat Cell Biol (fixed misattribution: NOT Khan) |
| tcga_scatac | sundaram_2024 | Sundaram et al. Cell Rep |
| gbm_tme_atlas_hra004942 | zhao_2024 | Zhao et al. Cell Rep Med |
| mathewson_lupien | mathewson_2021 | Mathewson Cell + Wang 2024 Sci Adv |
| gse138794_guo | guo_2020 | Guo et al. Nat Comm |
| guilhamon | guilhamon_2021 | Guilhamon et al. eLife |
| gse165037 | raviram_2023 | Raviram et al. 2023 PNAS 120:e2210991120 (Ren lab, UCSD — attribution corrected from Nagaraja/Monje) |

Mapping persisted at `data/matrix/slug_canonical.tsv`.

### Scoring

For each target cell_type (`pan_malignant`, `pan_myeloid`, 10 non-malignant
cell types, and 4 Neftel-state detail-only reports):

- `daigle_z_raw` — (target_mean − bg_mean) / bg_std across non-target columns,
  with CPM floor 0.5 applied before ranking (kills mathematical outliers where
  near-zero signal in both target and bg produce artifactually huge Z).
- `daigle_z_arm_adjusted` — **primary ranking column.** Same formula but the
  null distribution is computed per chromosome arm (centromere-split hg38),
  so an arm-wide gain (e.g. chr7 in GBM malignant) no longer inflates every
  peak on it against a genome-wide null. Used as the sort key in all reports.
- `chromosome_arm`, `chromosome_cnv_class` — p-arm/q-arm assignment and an
  indicator for the known-malignant chr7-gain / chr10-loss arms.
- `patient_fraction_chr7_elevated` — for malignant targets, the fraction of
  target patients where this peak falls on a chr7 arm with elevated signal;
  surfaces peaks whose malignant-specificity is actually driven by the
  patient-level CNV rather than cell-type biology.
- **Patient-level replication (primary):**
  - `n_target_patients` — count of target patients with nonzero signal
  - `n_target_cohorts` — distinct cohorts those patients span
  - `cohort_balance` — fraction of target signal from top-contributing cohort
    (close to 1 = cohort-dominated artifact)
- **AAV off-target safety:**
  - `offtarget_max_cpm` — worst non-target cell type's mean CPM
  - `offtarget_top_cell_type` — which non-target cell type leaks most
  - `offtarget_pooled_mean` — mean across non-target cell-type means
  - `n_offtarget_ct_above_1` — count of non-target cell types with CPM > 1
  - `aav_safety_ratio` — target / (offtarget_max + 0.1), directly actionable

### Reports folder layout

- `reports/` root — polished top-100 distal CSVs per report + combined
  xlsx/csv/parquet. Filter: distance_category ≠ promoter, target_mean_cpm ≥ 1.0,
  n_target_patients ≥ 2.
- `reports/detail/` — HTML reports with sortable tables (batch/leak highlights),
  annotated full CSVs, top50 distal CSVs, top30 distal BEDs, raw-score parquets,
  index.html gallery.
- `reports/detail/per_celltype/` — same for each cell type.

### Pre-release working versions (not shipped)

- v0.5 — CATLAS-only peak reference, cohort-level pseudobulks only, no
  off-target metrics. Replaced because the atlas can't claim to cover GBM
  chromatin if 74% of accessible peaks aren't in the reference set.
- v0.6 — pool rule fix (same peak reference, cohort-level still). Verified
  the pool rule logic but didn't ship due to pivot to patient-level.
- v0.7 — patient-level pseudobulks, CPM floor, per-non-target-cell-type
  contamination metrics. Became v2026-10-11 after reorganization + rename.
