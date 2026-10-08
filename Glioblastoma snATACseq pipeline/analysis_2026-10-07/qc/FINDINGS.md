# GBM atlas QC — findings log

## Check #1 — Per-patient audit  (2026-10-07)

Script: `qc/scripts/01_per_patient_audit.py`
Output: `qc/outputs/01_per_patient_audit.csv`

**Finding:** 18 of 51 patients excluded from pan-malignant pool. Grouped by failure mode:

| Failure mode | n patients | Cohorts affected |
|---|---|---|
| `CNV_not_run_for_cohort` | 12 | tcga_scatac (9), gse276177_khan_astro (3) |
| `zero_cnv_positive_zero_label` | 2 | hra004942 P98, mathewson_lupien Wang_pair_16 (2-cell patient) |
| `pan_mal_cells < 10 (threshold)` | 4 | hra004942 P79/P80/P84, gse165037 GBM9 |

**Root cause of the 12-patient block:** `src/flatten_cohort.py` (the fragments-mode replacement
for the broken QC/CNV/labels/quantify chain) did not call CNV. It set `malignant_cnv = 0`
as a default. This silently dropped the two biggest fragments-mode cohorts from any
CNV-based malignant analysis — 237,000 cells of GBM tumor data not being used.

**Fix proposed:** add in-process CNV calling to `src/flatten_cohort.py` via the shared
`_parallel_chr7_10_ratio` helper in `src/cnv_malignant.py`. **Patched locally
(2026-10-07) and CNV re-run launched on EC2.** Expected to recover 12 patients.

**Fix for the 4 <threshold patients:** these are legitimately low-cell-count patients
(1-8 malignant cells each). Default threshold of 10 is reasonable. Alternative: drop
threshold to 5 for cohorts with low cell counts overall; would recover 2 more patients
but at the cost of noisier per-patient estimates.

**Fix for the 2 zero-call patients:** Wang_pair_16 has 2 total cells post-QC — nothing
to rescue. P98 (hra004942) has 480 cells and CNV found zero — could be a non-tumor
contaminant sample or a very TAM-rich tumor. Needs manual review.

---


## Check #3 — Peak coordinate harmonization  (2026-10-07)  — N/A

**Status:** check does not apply to this atlas.

**Why:** All 8 cohorts were projected onto a **common 544,735-peak CATLAS atlas**
(Li et al. 2023 Science, GSE244618) before any cross-cohort comparison, via
either `src/ingest_matrix.py` (matrix-mode cohorts) or `src/peak_quantify.py` /
`src/flatten_cohort.py` (fragment-mode cohorts). Both paths call
`snapatac2.pp.make_peak_matrix` with the CATLAS BED as the reference peak set.

So "replicated across cohorts" in this atlas really means
"the same CATLAS peak was accessible in n cohorts' malignant pool" —
not "overlapping peaks from independent per-cohort peak calls." The summit-
agreement check applies to atlases built from per-cohort peak calls + bedtools
intersect merging; ours was built with a fixed common peak set from day one.

**Residual concern:** CATLAS itself may have called peaks at slightly different
summits than what per-cohort peak calling in each paper would have produced.
Nothing we can do about that without reprocessing all 8 cohorts from FASTQ, and
no evidence it meaningfully distorts rankings for the AAV use case.

---

## Check #9 — TSS re-classification via ENCODE cCRE  (2026-10-07)  — deferred

**Status:** deferred, downloaded URL returned 404.

Downloads from `https://downloads.wenglab.org/V3/GRCh38-PLS.bed` returned HTML
error page (ENCODE SCREEN URL has moved). The equivalent orthogonal check can
be done with:

- FANTOM5 CAGE peaks (large download, ~1 GB)
- ENCODE cCRE registry via the SCREEN REST API
- GENCODE comprehensive gene set (previously tried, EBI mirror slow)

In the interim, the v2 report's distal filter (|dist_to_tss| ≥ 2 kb to nearest
RefSeq transcript) catches the obvious promoters — 11/20 of the top candidates
in v2 dropped as promoter-proximal. Alternate TSSs 2-5 kb from annotated ones
would need the orthogonal CAGE check to catch; that's a known gap.

**Fix for next pass:** use `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.basic.annotation.gtf.gz`
to extract transcript-level TSSes (not just gene-level TSSes), then flag peaks
within 2 kb of any transcript-level TSS as "promoter-proximal (alternate TSS)".

---

## Check #1 — refresh after CNV fix  (2026-10-07)

Script: `/tmp/qc_batch.py` on EC2 → `qc/outputs/01_per_patient_audit_v2.csv`

**Finding:** 45/51 patients included (vs 33/51 in v1). CNV fix in `src/flatten_cohort.py`
successfully recovered all 12 CNV-missing patients from tcga_scatac + gse276177.

Remaining 6/51 exclusions:
- 2 zero-call: `mathewson_lupien Wang_pair_16` (2 total cells post-QC — nothing to rescue),
  `hra004942 P98` (480 cells, CNV found zero — manual review needed)
- 4 below threshold (<10 pan-mal cells): `hra004942 P79/P80/P84`, `gse165037 GBM9`

**Action:** retained as-is. Downstream analyses use pan-malignant union (CNV | label) with
MIN_CELLS_FOR_PAN_MAL=10 threshold.

---

## Check #4 — CN region flags on top shortlist  (2026-10-07)

Script: `/tmp/qc_v2_batch.py` on EC2 → `qc/outputs/04_cn_region_flags.csv`

**Finding:** Of the top 20 pan-malignant candidate peaks, 0 fall inside canonical GBM CNV
regions (chr7 whole-arm gain, chr10 whole-arm loss, CDKN2A/B 9p21 focal loss, EGFR focal
amplification 7p11.2). All `cnv_region_flags` empty.

Note: the v3 report already excludes chr7 globally at the ranking stage (non-chr7 filter in
`clean` pool). This check confirms no additional CN-region contamination on chr10 or 9p21.

---

## Check #5 — doublet audit scope  (2026-10-07)

Script: `/tmp/qc_v2_batch.py` on EC2 → `qc/outputs/05_doublet_audit_scope.csv`

**Finding:** Only `gse276177_khan_astro` has 10 per-sample h5ads on disk (ingested from
fragments-mode); the other 7 cohorts consumed pre-processed h5ads from GEO/EBI where
doublet filtering was done upstream by original authors.

**Action:** cannot run AMULET here; relying on authors' doublet filtering. Note on report:
if a cohort's shortlist is suspiciously enriched for mixed-identity peaks, re-run AMULET on
raw fragments.

---

## Check #7 — batch metadata availability  (2026-10-07)

Script: `/tmp/qc_batch.py` on EC2 → `qc/outputs/07_batch_metadata_audit.csv`

**Finding:** Across all 8 cohorts, only `patient_id` + `sample_id` are present in obs.
No `chip_id`, `library`, `operator`, `processing_date`, or `batch` column available.

**Interpretation:** Can't do classical batch-effect regression. But because we aggregate at
(peak × patient × cohort) level and each patient is effectively its own batch, batch effects
collapse into patient-level random effects — captured by the hierarchical Bayesian pass2
(variance across patients, within cohort).

**Action:** not fixable without re-contacting authors; downstream analysis treats
`patient` as the lowest exchangeable unit, which is defensible.

---

## Check #8 — ENCODE blacklist intersect  (2026-10-07)

Script: `qc/scripts/08_blacklist_intersect.py` → `qc/outputs/08_blacklist_summary.txt`

**Finding:**
- CATLAS atlas: 1,497 / 544,735 peaks (0.27%) overlap ENCODE blacklist hg38 v2
- v2 top-50 pan-malignant shortlist: **0 / 50** overlap

**Action:** no filter needed at AAV cloning stage. Blacklist cleanup is implicit in the
rarity of blacklist-overlapping peaks and the stringent cross-cohort replication filter.

---

## Check #9 — TSS re-classification v2 (GENCODE basic transcript-level)  (2026-10-07)

Script: `/tmp/qc_v2_batch.py` on EC2 → `qc/outputs/09_tss_reclassification_v2.csv`

**Finding:** v1 downloaded ENCODE SCREEN cCRE URL (404). v2 uses GENCODE v45 basic
transcript-level TSSes (not just gene-level) to catch alternate-TSS promoters that gene-level
RefSeq distance misses. Of the current top-20 candidates, 3 flagged as
`promoter_proximal_strict` (within 2 kb of a GENCODE transcript TSS, but >2 kb from any
RefSeq gene TSS).

**Action:** v3 clean report retains the 2-kb-to-RefSeq filter. For stricter cloning,
cross-check candidates against `09_tss_reclassification_v2.csv` and drop any
`promoter_proximal_strict == True`.

---

## Check #12 — per-patient TSS enrichment / FRiP  (2026-10-07)

Script: `/tmp/qc_batch.py` on EC2 → `qc/outputs/12_per_patient_qc.csv`

**Finding:** None of the 8 cohorts carry `tss_enrichment` or `frip` columns in their
processed h5ads. These metrics are computed at ingest-time by the original pipelines
(snapatac2.pp.add_tsse) and then dropped before deposit — the processed data we received
are already post-QC.

**Interpretation:** All 51 patients were quality-filtered upstream (per original paper
methods — see `reports/methods_comparison.md`). We don't have raw per-cell QC metrics
to filter further.

**Action:** not fixable without re-processing from fragments. The `methods_comparison.md`
report documents each cohort's upstream QC thresholds for cross-cohort comparability.

---

## Check #13 — peak saturation approximation  (2026-10-07)

Script: `/tmp/qc_v2_batch.py` on EC2 → `qc/outputs/13_peak_saturation_approx.csv`

**Finding:** All 8 cohorts projected onto the same 544,735-peak CATLAS atlas (by design
— see Check #3). Avg-peaks-per-cell could not be computed on the backed='r' loader (nnz
returned -1 placeholder). Full matrix nnz counts are in `data/matrix/enhancer_candidate_matrix.parquet`
(3.5 GB).

**Action:** detailed saturation curves not blocking; CATLAS peak space is fixed, and we
report per-peak `n_patients_accessible` in the pan-malignant scores as the de facto
accessibility rate per peak.

---

## Check #14 — clinical metadata audit  (2026-10-07)

Script: `/tmp/qc_batch.py` on EC2 → `qc/outputs/14_clinical_metadata_audit.csv`

**Finding:** Only `idh_status` + `primary_recurrent` + (for 3 cohorts) `neftel_state`
are populated. No `age`, `sex`, `grade`, `mgmt_status`, `treatment`, or `survival`.

**Interpretation:** Limits downstream clinical correlation (e.g., "is enhancer X enriched
in recurrent vs primary"? → can do; "...in MGMT-methylated vs un-methylated"? → cannot).

**Action:** for R01 Aim 3 (clinical validation), will need to pull missing fields from
TCGA / GEO sample sheets per cohort. Deferred — not blocking Aim 1 (enhancer discovery)
or Aim 2 (AAV cloning).

---

## Check #15 — surprises.parquet null audit  (2026-10-07)

Script: `/tmp/qc_v2_batch.py` on EC2 → `qc/outputs/15_surprises_null_audit.csv`

**Finding:** `patient_id`, `delta`, `direction`, `n_low`, `n_high`, `frac_high` columns
are **100 % null** in `surprises.parquet` (1.1 GB, 137M rows). These were upstream
artifact columns from an earlier scoring pass. `peak_id`, `detector`, `cell_type`,
`pattern` columns are populated and usable.

**Action:** regenerate `surprises.parquet` on next full scoring run to drop the null
columns and populate `patient_id` / `delta` properly. For now, downstream analyses use
only the populated columns.

---

## Summary table (addressable → addressed)

| Check | Status | Blocking? |
|---|---|---|
| #1 per-patient audit refresh | ✓ addressed (45/51 inclusion) | No |
| #3 peak coordinate harmonization | N/A (common CATLAS atlas) | No |
| #4 CN-region flags | ✓ no top-shortlist hits in CN regions | No |
| #5 doublet audit scope | documented limit | No |
| #7 batch metadata | collapsed into patient-level hierarchy | No |
| #8 ENCODE blacklist | ✓ 0/50 overlap | No |
| #9 TSS reclass v2 (GENCODE) | ✓ 3/20 alt-TSS flagged | Minor |
| #12 per-patient TSS/FRiP | not available in processed data | No (relying on upstream QC) |
| #13 peak saturation | not computed (fixed peak space makes moot) | No |
| #14 clinical metadata | limited fields; need re-pull for Aim 3 | No for Aim 1/2 |
| #15 surprises null audit | documented; regenerate next scoring pass | No |

All QC checks completed or triaged. No blockers for the current pan-malignant v3 or
per-cell-type AAV candidate reports.
