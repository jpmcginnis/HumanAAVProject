# CHANGELOG — Human Glioblastoma Enhancer Atlas (snATAC-seq)

Semantic versioning: major = schema change, minor = new data / new reports, patch = fixes with same schema.

## v1.0.0 — 2026-10-08 (deposit-ready)

**Headline:** First Zenodo-eligible release. Schema frozen (see `SCHEMA.md`). Fixes
the CNV-missing bug that excluded 12 patients in fragments-mode cohorts.

### Data

- Pan-malignant pool rebuilt from 45/51 patients × 8 cohorts (575,971 cells; CNV | label union).
  Pre-fix pool was 33/51 (438,994 cells from gbm_space only).
- Matrix: `pan_malignant_matrix_v4.parquet` (10.4M rows), `pan_malignant_scores_v4.parquet` (544,735 peaks).
- Candidate scores + enhancer candidate matrix unchanged in content; schema is frozen.

### Code (`src/flatten_cohort.py`)

- **FIX:** `_parallel_chr7_10_ratio` now runs in-process during flattening for fragments-mode cohorts
  (tcga_scatac, gse276177). Previously `malignant_cnv` silently defaulted to 0, dropping these
  cohorts from any CNV-based malignant analysis.

### Reports

- Three HTML reports: `pan_malignant_report_v3.html` (primary), `per_celltype_report.html`
  (10 cell types including pan_malignant), `pan_myeloid_report.html` (TAM ∪ microglia).
- **Dropped the 100 kb upper distance cap** on all shortlists. For AAV targeting, the enhancer
  is extracted from genomic context; native-genome distance is irrelevant. 2 kb floor retained
  (excludes promoter-proximal non-enhancers). Added `distance_category` as context annotation
  (near / distal / far-distal / gene-desert). Impact: pan-malignant AAV-clean pool 461 → 617
  (+34%), pan-myeloid SHARED pool 367 → 380.
- **Added Daigle Z-score column** (`z_daigle = (strength_target − mean_other) / sd_other`,
  `passes_daigle_z2`) to all three reports. Z ≥ 2 is the Allen Institute Armamentarium
  criterion for cell-type specificity. 11,641 peaks atlas-wide pass Z ≥ 2 for pan-malignant.

### QC

- 11 QC checks compiled in `qc/FINDINGS.md`. All triaged or resolved. No blockers.
- Per-patient audit regenerated post-fix as `per_patient_audit_v2.{csv,html}`.
  Pre-fix audits removed to avoid confusion.

### Documentation

- Deposit-ready docs added: `README.md` (project orientation), `SCHEMA.md` (frozen schema),
  `cohort_attribution.csv` (reviewer-ready attribution), `LICENSE_DATA` (CC-BY-4.0),
  `LICENSE_CODE` (MIT), `zenodo.json` (metadata).

### Tests

- `tests/` added with **27 unit tests** covering schema contract (parquet columns + dtypes),
  computations (distance_category, Daigle Z, selectivity, aav_score), and filter logic
  (pool size, bounds, cell-type coverage). Run with `pytest tests/`.

### Known slug misnomer

- Cohort slug `gse276177_khan_astro` is actually **Sojka et al. 2025 Nature Cell Biology**,
  not Khan. The slug is kept for internal consistency with historic column values but
  renaming would require regenerating the matrix. See `cohort_attribution.csv` for the
  real citation.

## v0.3 — 2026-10-07 (internal — pre-CNV-fix pan_malignant v3 report)

- First pan-malignant report with distal + selectivity + non-chr7 + non-housekeeping filters.
- 33/51 patients (pre-fix). 461 clean AAV candidates.
- Found and documented the fragments-mode CNV bug; fix implemented and v4 pool rebuilt same day.

## v0.2 — 2026-10-06

- 8-cohort matrix build complete (`enhancer_candidate_matrix.parquet`, 209M rows).
- Hierarchical Bayesian pass-2 scoring on top 2000 per cell type.
- First pan-malignant report (v1).

## v0.1 — 2026-10-05 and earlier

- Pipeline bootstrap on EC2 (r6i.4xlarge + 500 GB EBS at `/data`).
- CATLAS 544K-peak reference + GBmap integration.
- Per-cohort QC, CNV, label-transfer scaffolding.
- Pipeline source in `src/*.py` on GitHub (`jpmcginnis/HumanAAVProject`).
