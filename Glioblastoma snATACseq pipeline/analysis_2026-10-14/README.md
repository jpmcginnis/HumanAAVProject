# Harmonized human glioblastoma chromatin accessibility atlas (snATAC-seq)

**v2026-10-14 — 2026-10-14**
**Scope:** 48 patients (unique (cohort, patient) pseudobulks with ≥50 cells; counted from the matrix, not the raw submission counts), 8 published GBM scATAC/snMultiome cohorts, 1.52M nuclei, 1,958,753 union peaks.
**Deliverable focus:** AAV cell-type-targeting enhancer peak candidates with patient-level replication, cohort batch sanity, per-non-target-cell-type contamination, and arm-level CNV correction.
**License:** CC-BY-4.0 for data (`LICENSE_DATA`), MIT for code (`LICENSE_CODE`)
**Contact:** JP McGinnis (jpmcginnis1@gmail.com), Baylor College of Medicine, Department of Neurosurgery.
**Pipeline source:** <https://github.com/jpmcginnis/HumanAAVProject> → `Glioblastoma snATACseq pipeline/`.

---

## Quick start (one screen)

**What this is.** A per-cell-type scored snATAC-seq enhancer-candidate atlas across 8 published GBM cohorts, ranked for AAV capsid cell-type restriction.

### Where to find each piece

| Component | Location | How to access |
|---|---|---|
| Pipeline source (scripts, tests, docs) | **This GitHub repo** | `git clone` |
| Union peak reference (`union_peak_reference.bed`, 84 MB) | **Zenodo DOI [10.5281/zenodo.23252324](https://doi.org/10.5281/zenodo.23252324)** | `scripts/zenodo_upload.py` or web |
| Per-peak provenance (`peak_provenance.parquet`, 55 MB) | ↑ same Zenodo DOI | ↑ |
| Patient × cell-type pseudobulk matrix (`enhancer_candidate_matrix.parquet`, 366 MB) | ↑ same Zenodo DOI | ↑ |
| Greenwald fine-label backup (`greenwald_fine_pseudobulk_backup.parquet`, 495 MB) | ↑ same Zenodo DOI | ↑ |
| Per-report scored parquets (16 files, 19–35 MB each) | ↑ same Zenodo DOI | ↑ |
| Polished top-100 CSVs + combined HTML/xlsx/parquet | ↑ same Zenodo DOI | ↑ |
| Original per-cohort snATAC h5ads | Original studies' deposits | see `cohort_attribution.csv` for GEO/ArrayExpress/HRA/dbGaP accessions |

### Regenerate the deliverables from the Zenodo data

```bash
# 1. Download the Zenodo deposit; place its contents under data/matrix/ and reports/detail/
#    (zenodo_upload.py also pulls deposits by DOI)
python3 scripts/zenodo_upload.py --download 10.5281/zenodo.23252324 --out .

# 2. Rebuild the polished deliverables from the per-report parquets
python3 scripts/build_reports.py          # annotated CSVs + per-report HTMLs + top30 distal BEDs
python3 scripts/build_polished_reports.py # 12 top-100 polished CSVs + top100_distal_combined.{csv,xlsx,parquet} + top10_per_arm_combined
python3 scripts/build_combined_html.py    # sortable + filterable + hover-definition HTML
python3 -m pytest tests/                  # 17/17 regression tests — see REGEN_LOG.md for the latest receipt
```

The scoring step (`enhancer_candidate_matrix.parquet` → 16 per-report scored parquets) runs on the host where the full matrix lives; see `CHANGELOG.md` for the scoring pipeline description.

---

## Design decisions (judgment calls behind this release)

Grouped by where in the pipeline the call was made. Full rationale and the versions where each decision was made or revised are in `CHANGELOG.md`.

### Scope and inclusion

- **Public data only, no deferrals.** Every cohort that was publicly available at build time was included (open GEO/ArrayExpress/HRA and dbGaP-authorized access). No cohort was skipped on grounds that "another cohort already covers the same region or cell type" — study diversity is what drives cross-cohort replicability.
- **GBM snATAC / snMultiome only.** No bulk ATAC, no scATAC from non-GBM brain tumors, no methylation (ATAC and methylation specificity criteria don't share a null and would need separate ranking machinery).
- **Mixed fresh + postmortem.** Both are included; the fresh-vs-postmortem axis is not modeled explicitly and is called out as a limitation. The alternative (dropping postmortem) would collapse the atlas to ~3 cohorts.
- **No distance cap on candidate reporting.** The 100 kb promoter-distance cap used for many enhancer atlases is deliberately omitted: in AAV cassettes the enhancer is extracted from genomic context, so native-genome distance to the nominally regulated gene is informational, not a filter. See `nearest_gene` / `dist_to_tss_signed` which are annotations only.

### Peak reference

- **hg38 throughout.** All cohorts coordinate-converted to GRCh38 (hg38); cross-checked by CATLAS peak overlap per cohort.
- **Union of per-cohort MACS3 calls ∪ CATLAS, NOT CATLAS-only.** 74% of the union set (1,453,136 of 1,958,753) is de-novo peaks not in CATLAS. A CATLAS-only reference would discard most of the biology this atlas is designed to surface. (Rejected as v0.5; see CHANGELOG.)
- **No transitive chaining when merging.** Union is a single-pass bedtools intersect across cohorts' called peaks, so merged intervals stay at single-cohort peak-call precision (mean width 594 bp, max 6.7 kb). `peak_provenance.parquet` traces every merged peak back to its contributing source coordinates.

### Pseudobulks

- **Patient-level, not cohort-level.** Columns are `<cohort>__<patient_id>__<label>` so patient is the unit of replication (`n_target_patients`) and cohort is recoverable for batch sanity (`n_target_cohorts`, `cohort_balance`). (Cohort-level pseudobulks rejected as v0.6.)
- **Malignant pool rule (CNV-aware).** A cell with `malignant_cnv == 1` OR `cell_type == "malignant_unresolved"` is bucketed as `malignant`, not as its gbmap cell-type label. Fixes the "hidden malignant" contamination of astrocyte / OPC / neuron pseudobulks that affected earlier versions.
- **MIN_CELLS = 50 per (patient × label) group.** Smaller groups skipped as too noisy to pseudobulk.
- **Greenwald label transfer via Multiome barcode sharing.** Greenwald is 10x snMultiome, so the snRNA barcodes are the same as snATAC barcodes — labels transfer directly from the paired snRNA h5ad into the ATAC h5ad with no reference-anchor method required.
- **Coarse-label remap for Greenwald; Neftel subtypes retained as detail; fine labels preserved for future.** Greenwald's 99 fine labels (AC-progenitor-like-1-4, OPC-NPC-like-1-2, Pro-inflammatory TAMs, Astrocyte-like TAMs, etc.) are remapped to the coarse schema used by the other 7 cohorts so cross-cohort reports work; the 4 Neftel-state buckets (`malignant_AC` / `MES` / `OPC` / `NPC`) are preserved as separate drill-down reports in `detail/per_celltype/`; and the full 99-label pseudobulk is preserved at `data/matrix/greenwald_fine_pseudobulk_backup.parquet` so future releases can re-remap and generate super-fine reports without re-reading the 10.5B-nnz ATAC h5ad.

### Scoring

- **Daigle Z ≥ 2 (Allen Armamentarium specificity criterion).** `(target_mean − μ_bg) / σ_bg` across non-target columns. Z ≥ 2 is the Allen Institute's published cutoff for cell-type-specific regulatory elements.
- **Arm-adjusted Z is the primary ranking column.** Null distribution computed per hg38 chromosome arm, so chr7-gain / chr10-loss in malignant cells no longer inflates every peak on those arms against a genome-wide null. Raw Z is retained alongside for transparency.
- **Per-patient CNV flag for malignant targets.** `patient_fraction_chr7_elevated` surfaces peaks whose apparent malignant-specificity is driven by patient-level CNV rather than cell-type biology.
- **CPM floors: 0.5 raw, 1.0 polished.** Floor applied before ranking kills mathematical outliers where near-zero signal in both target and background produces artifactually huge Z scores. The stricter 1.0 floor is used in the headline polished top-100 view.
- **MIN_PATIENTS ≥ 2 for scoring.** Primary replication filter; applied at the per-report parquet level, not just in the polished view.

### Reports

- **12 headline reports in the top-level polished view.** 10 non-malignant cell types + `pan_malignant` + `pan_myeloid`. `pan_malignant` is the single primary malignant deliverable.
- **4 Neftel-state malignant reports as detail-only.** `malignant_AC` / `MES` / `OPC` / `NPC` live in `detail/per_celltype/` as drill-down resources. They are intentionally excluded from the headline combined view because they are Greenwald-dominated (only Greenwald has Neftel resolution across the atlas).
- **Nearest-TSS annotation only, not a regulatory claim.** `nearest_gene`, `dist_to_tss_signed`, `distance_category` are annotations only. No HiC / enhancer-promoter linking.

## For a reviewer: four files that tell the whole story

1. **`reports/top100_distal_combined.html`** — scrollable, sortable, filterable table of 1200 polished top-100 distal AAV enhancer candidates across the 12 headline reports (10 non-malignant cell types + `pan_malignant` + `pan_myeloid`). Hover any column header or filter label for its equation + plain-language definition. Also written as `.xlsx`/`.csv`/`.parquet`.
2. **`reports/detail/index.html`** — browsable gallery; each tile links to a per-report HTML with sortable table and automatic batch/leak highlighting. Includes drill-down detail reports for the 4 Neftel-state malignant subtypes (`malignant_AC`/`MES`/`OPC`/`NPC`) in `detail/per_celltype/`.
3. **`cohort_attribution.csv`** — one row per cohort: citation, accession, data-use conditions.
4. **`CHANGELOG.md`** — rationale for the design decisions (peak reference, pool rule, patient-level replication, CPM floor, off-target contamination, arm-adjusted Z).

## Scoring metrics (every peak in every report)

| metric | meaning |
|---|---|
| `target_mean_cpm` | strength in target cell type (patient-averaged) |
| `bg_mean_cpm` | mean across non-target columns |
| `daigle_z_raw` | `(target_mean − μ_bg) / σ_bg` — Allen Armamentarium specificity criterion against a genome-wide null (Z ≥ 2) |
| `daigle_z_arm_adjusted` | **primary ranking column.** Same formula but μ/σ computed per chromosome arm. **Standardized per chromosome arm and not directly comparable across arms or to `daigle_z_raw`** — a z_arm_adj of 97 does not mean "97× more specific than background", it means that peak sits 97 arm-conditional SDs above the arm's non-target mean. |
| `n_target_patients` | **primary replication** — unique (cohort, patient_id) tuples with signal in any target label |
| `n_target_pseudobulks` | distinct patient × label pseudobulk columns with signal. Only differs from `n_target_patients` for `pan_malignant` and `pan_myeloid`, where a patient contributes multiple target columns (e.g. a Greenwald patient has 5 malignant-label columns) |
| `n_target_cohorts` | **batch sanity** — distinct cohorts those patients span |
| `cohort_balance` | fraction of target signal from top-contributing cohort (close to 1 = cohort-dominated) |
| `greenwald_fraction` | fraction of target signal (CPM-normalized) contributed by Greenwald. **Important for `pan_malignant`**: Greenwald is the only cohort with Neftel-split malignant labels (5 label types × 12 patients = 60 target columns vs. the other cohorts' single `malignant` label × ~28 patients = ~28 columns), so the pan_malignant denominator is Greenwald-weighted (median greenwald_fraction = 0.78). Pan_myeloid is cross-cohort balanced (median 0.25). |
| `offtarget_max_cpm` | **AAV safety** — worst non-target cell type's mean CPM |
| `offtarget_top_cell_type` | which non-target lineage leaks most |
| `offtarget_pooled_mean` | mean across non-target cell-type means |
| `n_offtarget_ct_above_1` | number of non-target cell types with mean CPM > 1 |
| `aav_safety_ratio` | `target_mean_cpm / offtarget_max_cpm` — directly actionable |
| `chromosome_arm` | hg38 centromere-split arm (e.g. `7p`, `7q`, `10p`, `10q`) |
| `chromosome_cnv_class` | `chr7_gain` / `chr10_loss` / `none` — known-malignant arm-level CNVs |
| `patient_fraction_chr7_elevated` | (malignant targets only) fraction of target patients where this peak falls on an arm with per-patient elevated signal — flags CNV-driven apparent specificity |
| `nearest_gene`, `dist_to_tss_signed`, `distance_category` | gene annotation (promoter / proximal / distal / far-distal) |

A good AAV capsid candidate = high `target_mean_cpm`, high `daigle_z_arm_adjusted`, `n_target_patients` ≥ 10 across ≥ 3 cohorts, low `offtarget_max_cpm`, low `n_offtarget_ct_above_1`, `cohort_balance` ≤ 0.7, and (for malignant targets) low `patient_fraction_chr7_elevated`.

## Reports structure

**12 headline reports** go into the top-level polished view (`reports/top100_distal_combined.*`):
10 non-malignant cell types (astrocyte, endothelial, GABA_neuron, microglia, neuron, oligodendrocyte, OPC, T_cell, TAM, unassigned) + **`pan_malignant`** (all malignant pooled across Neftel states + unresolved) + `pan_myeloid` (TAM ∪ microglia).

**4 Neftel-state detail reports** (`malignant_AC`, `malignant_MES`, `malignant_OPC`, `malignant_NPC`) are available in `reports/detail/per_celltype/` for drill-down. They are intentionally excluded from the headline combined view — `pan_malignant` is the primary malignant deliverable, with the Neftel split as a secondary resource.

The fine-grained Greenwald 99-label pseudobulk (full annotation granularity: AC-progenitor-like-1-4, OPC-NPC-like-1-2, pro-inflammatory TAMs, astrocyte-like TAMs, CD8+ cytotoxic T cells, Exc L4-5 IT, etc.) is preserved at `data/matrix/greenwald_fine_pseudobulk_backup.parquet` so future releases can re-remap and generate those super-fine reports without re-reading the 10.5B-nnz Greenwald h5ad.

## Directory contents

```
.
├── README.md                               ← you are here
├── CHANGELOG.md                            ← version history + rationale
├── SCHEMA.md                               ← frozen deposit schema
├── HOW_TO_WORK.md                          ← playbook for future analysis projects
├── cohort_attribution.csv                  ← per-cohort citation + accession
├── zenodo.json                             ← deposit metadata
├── LICENSE_DATA  LICENSE_CODE
│
├── data/
│   └── matrix/
│       ├── union_peak_reference.bed              1.96M peaks (bed4: chrom, start, end, source list)
│       ├── peak_provenance.parquet               per-peak contributing cohort + source coords
│       ├── enhancer_candidate_matrix.parquet     1.96M peaks × 337 (patient × cell_type) pseudobulks across 8 cohorts
│       ├── greenwald_fine_pseudobulk_backup.parquet  99-label fine-grained Greenwald pseudobulk (for future re-remapping)
│       └── slug_canonical.tsv                    legacy → canonical slug map (incl. Sojka ≠ Khan fix)
│
├── reports/                                ← polished deliverables
│   ├── top100_distal_combined.html         scrollable sortable all-reports view with hover tooltips
│   ├── top100_distal_combined.{csv,xlsx,parquet}
│   ├── top10_per_arm_combined.{csv,xlsx}   top-10 per chromosome arm per report (CNV-diversity view)
│   ├── <celltype>_top100.csv               × 10 non-malignant cell types + pan_malignant + pan_myeloid
│   └── detail/                             ← full artifact set
│       ├── index.html                      gallery
│       ├── pan_malignant.{html,parquet}    + _annotated.csv + _top100.csv + _top50_distal.csv + _top30_distal.bed
│       ├── pan_myeloid.* (same)
│       ├── top50_combined.{csv,xlsx,parquet}        unfiltered top-50 any-distance across all 16 reports
│       ├── top50_distal_combined.{csv,xlsx,parquet} top-50 distal only
│       └── per_celltype/
│           ├── per_celltype_<ct>.html         per-cell-type HTML with sortable table
│           ├── per_celltype_<ct>.parquet      raw score
│           ├── per_celltype_<ct>_annotated.csv full peak list with annotations
│           ├── per_celltype_<ct>_top100.csv   top-100 any-distance
│           ├── per_celltype_<ct>_top50_distal.csv
│           └── per_celltype_<ct>_top30_distal.bed
│           (14 cell types: 10 non-malignant + 4 Neftel-state malignant drill-downs)
│
├── scripts/
│   ├── build_reports.py                    regenerates reports/detail/ from parquets (annotation + per-report HTML + top50/top30 distal)
│   ├── build_polished_reports.py           regenerates reports/*_top100.csv + top100_distal_combined.{csv,xlsx} + top10_per_arm_combined
│   ├── build_combined_html.py              regenerates reports/top100_distal_combined.html (sortable + filterable + hover tooltips)
│   ├── pileup_region.py                    per-cohort fragment pileup for any region
│   └── zenodo_upload.py                    resumable Zenodo CLI uploader
│
├── reference/                              hg38 refGene BED, ENCODE blacklist (unchanged)
├── workflow/
│   └── REPLICATE_FOR_NORMAL_CORTEX.md      protocol for the normal-cortex equivalent atlas
└── tests/                                  pytest regression tests (15/15 green: schema + filters + computations)
```

## How to regenerate the deliverables

```bash
python3 scripts/build_reports.py            # annotated CSVs + per-report HTMLs + top50/top30 BED in reports/detail/
python3 scripts/build_polished_reports.py   # top-100 polished CSVs + top100_distal_combined.{csv,xlsx} + top10_per_arm_combined
python3 scripts/build_combined_html.py      # top100_distal_combined.html with sortable table + hover tooltips
python3 -m pytest tests/                    # 15/15 regression tests
```

## How to interrogate a specific peak

```bash
python3 scripts/pileup_region.py chr5:44411500-44412700 --out pileup.png
```

Pulls each cohort's raw fragments overlapping the region, builds per-cohort coverage tracks, writes a stacked PNG. Useful when a merged union interval needs sanity-checking against the underlying signal in each contributing cohort.

## Known limitations

- **pan_malignant is Greenwald-weighted** (median `greenwald_fraction` = 0.78). Greenwald is the only cohort with Neftel-resolved malignant labels, which gives it 60 of ~88 total pan_malignant target columns. The scoring is still cross-cohort-validated (seven other cohorts contribute signal), but reviewers should treat the pan_malignant top-100 as "Greenwald-primary, cross-cohort-corroborated" rather than fully cohort-balanced. Use the `greenwald_fraction` column to filter if a strictly cross-cohort view is needed.
- **TAM top-5 includes chrX pile-up driven by Greenwald.** Three of the top-5 TAM peaks (ZCCHC12, SASH3, ARMCX2) are on chrXq with Greenwald contributing 61–69% of signal; only 3 cohorts have myeloid-label pseudobulks on chrX. Greenwald's snMultiome gives cleaner TAM labels than other cohorts' label-transfer labels, which inflates its apparent contribution at low-coverage loci. Pattern does not reflect sex confounding; Greenwald's TAM contribution is uniform across its 12 patients.
- **sundaram_2024 sample X008** skipped in peak-matrix build due to a known snapatac2 Rust panic on its obs structure. 8/9 sundaram samples contribute. Explicitly: of 50 patients across the 8 cohorts' deposit records (12 greenwald + 3 sojka + 9 sundaram + 11 zhao + 5 mathewson + 3 guo + 4 guilhamon + 3 raviram), **48 made it into scored pseudobulks** — sundaram X008 excluded (snapatac2 panic) and raviram GBM9 excluded (MIN_CELLS=50 floor). See CHANGELOG §"Patient-count provenance" for the audit trail (the zhao row was transcribed as 7 patients in a prior cohort_attribution.csv — actual count from h5ad is 11, 4-patient attribution error corrected in v2026-10-11).
- **Arm-adjusted Z nulls are stable at small-arm sizes.** We checked whether the per-arm null estimation suffered on tiny arms (chr21p=47 peaks, chr22p=9, chrYp=34, chrYq=83, chr15p=1 in the matrix); the top-100 of each report has at most 3 peaks on chr21q and ≤2 on chr22q, with no peaks from the smaller arms, so the σ estimates are not driving ranking.
- **Fresh-tissue vs postmortem.** Cohorts mix fresh surgical specimens with postmortem tissue; batch effects attributable to this split are not explicitly modeled. `cohort_balance` partially controls, but residual bias is possible.
- **Replication count metric: detection, not strength.** `n_target_patients` counts patients with any nonzero signal. A peak weakly detected in 15 patients is treated the same as one strongly detected in 15. The CPM axis carries the strength component.
- **Nearest-TSS annotation does not imply regulatory contact.** `nearest_gene` / `dist_to_tss_signed` are informational labels only — no HiC / enhancer-promoter linking.
- **Neftel-state assignment is cohort-dependent.** Only Greenwald's snRNA labels resolved to Neftel states; the other 7 cohorts contribute to `pan_malignant` but not to the `malignant_AC/MES/OPC/NPC` detail reports. Treat the Neftel detail reports as Greenwald-dominated (secondary resource).
