# A harmonized human glioblastoma enhancer atlas (snATAC-seq, 51 patients, 8 cohorts, 1.52M nuclei)

**Deposit title:** A harmonized human glioblastoma enhancer atlas (snATAC-seq, 51 patients, 8 cohorts, 1.52M nuclei)
**Short name:** GBM Enhancer Atlas v1.0.0
**Deliverable focus:** AAV cell-type-targeting peak candidates (pan-malignant, per-cell-type, pan-myeloid) with Daigle Z-score cell-type-specificity annotation.
**Version:** 1.0.0 (2026-10-08)
**License:** CC-BY-4.0 for data (`LICENSE_DATA`), MIT for code (`LICENSE_CODE`)
**DOI:** pending Zenodo deposit
**Pipeline source:** <https://github.com/jpmcginnis/HumanAAVProject> → `Glioblastoma snATACseq pipeline/`
**Contact:** JP McGinnis (jpmcginnis1@gmail.com), Baylor College of Medicine, Department of Neurosurgery.

---

## What this is

A harmonized AAV-targeting enhancer atlas built from **8 published scATAC-seq GBM cohorts**
(51 patients, 1.52M nuclei, 544,735 CATLAS peaks), with:

- A **pan-malignant** peak shortlist (575,971 cells across 45/51 patients; CNV | marker-peak-label union).
- **Per-cell-type** peak shortlists for all 10 non-noise cell types (microglia, TAM, astrocyte, oligodendrocyte, OPC, neuron, GABA_neuron, endothelial, T cell, pan-malignant).
- A **pan-myeloid** peak shortlist (TAM ∪ microglia) — because for AAV targeting of the myeloid compartment there's no therapeutic reason to split those apart.
- All peaks carry the **Allen Institute Armamentarium Daigle Z-score** (`z_daigle`, `passes_daigle_z2`) so results are reportable in the field-standard vocabulary.

The deposit is **ready to query** without any pipeline rerun: all shortlists are in CSV + BED, all underlying scores are in `data/matrix/*.parquet`, HTML reports render in any browser.

## For a reviewer: four files that tell the whole story

1. **`reports/pan_malignant_report_v3.html`** — primary deliverable. Headline numbers, methods, caveats, top-50 cloning shortlist.
2. **`cohort_attribution.csv`** — one row per cohort: citation, accession, data-use conditions, usage justification.
3. **`qc/FINDINGS.md`** — 11 QC checks (per-patient audit, ENCODE blacklist, CN-region flags, TSS reclass, batch metadata, etc.); all triaged or resolved.
4. **`reports/methods_comparison.md`** — per-cohort snATAC methods extracted verbatim from source papers (tissue harvest → nuclei isolation → library prep → sequencing). Flags for low-cell cohorts and orthogonal technology.

## Directory contents

```
.
├── README.md                                 ← you are here
├── SCHEMA.md                                 ← frozen column definitions for every parquet/CSV
├── CHANGELOG.md                              ← version history (v0.1 → v1.0.0)
├── cohort_attribution.csv                    ← reviewer-ready cohort metadata table
├── zenodo.json                               ← Zenodo deposit metadata (upload with the record)
├── LICENSE_DATA                              ← CC-BY-4.0
├── LICENSE_CODE                              ← MIT
│
├── data/
│   ├── matrix/
│   │   ├── MATRIX_SCHEMA.md                  ← detailed schema for enhancer_candidate_matrix
│   │   ├── pan_malignant_matrix_v4.parquet   (188 MB, 10.4M rows, peak × patient × cohort — pan-mal pool)
│   │   ├── pan_malignant_scores_v4.parquet   (14 MB, 544,735 peaks — per-peak summary)
│   │   ├── panmal_v4_run.log                 (log of pan-malignant v4 build)
│   │   ├── candidate_scores.parquet          (134 MB, 5.6M peaks × cell_type scores with pass-2 posterior for top 2000 each)
│   │   ├── enhancer_candidate_matrix.parquet (3.5 GB, 209M rows — atlas workhorse, peak × patient × cohort × cell_type)
│   │   ├── surprises.parquet                 (1.1 GB, pattern detector outputs — some columns null, see FINDINGS #15)
│   │   ├── top_candidates.parquet, provenance.parquet, dataset_issues.parquet
│   └── processed_h5ads/                      (69 GB, per-cohort CATLAS-quantified h5ads with CNV + labels)
│       ├── gbm_space/                        (62 GB; 12 patients, 1.04M cells)
│       ├── gse276177_khan_astro/             (2.9 GB; see slug note below — actually Sojka et al. 2025)
│       ├── sundaram_gbm/                     (2.1 GB; = tcga_scatac in the matrix cohort column)
│       ├── gbm_tme_atlas_hra004942/          (1.1 GB; low cells-per-patient)
│       ├── mathewson_lupien/, gse138794_guo/, guilhamon/, gse165037/
├── reports/
│   ├── pan_malignant_report_v3.html          ← PRIMARY (45 patients, methods caveats, Daigle Z)
│   ├── per_celltype_report.html              (10 cell types incl. pan_malignant top section)
│   ├── pan_myeloid_report.html               (TAM ∪ microglia, Table B = SHARED is primary)
│   ├── methods_comparison.md                 (per-cohort methods appendix)
│   ├── pan_malignant_top50_distal_selective_v3.csv  ← AAV cloning shortlist
│   ├── pan_malignant_top100_all_v3.csv              ← fuller top-100 (includes promoters/chr7 for reference)
│   ├── pan_malignant_top30_distal_selective_v3.bed  ← GRCh38 BED, Benchling-ready
│   ├── pan_myeloid_top50_shared.csv + _top30_shared.bed   ← PRIMARY pan-myeloid
│   ├── pan_myeloid_top50_distal_selective.csv + _top30.bed  ← strict pan-myeloid
│   ├── per_celltype/*_top50.csv + *_top30.bed  (one pair per cell type)
│   ├── per_patient_audit_v2.html + .csv       ← post-fix 45/51 audit
├── qc/
│   ├── FINDINGS.md                           ← 11 QC checks, synopsis + triage
│   ├── scripts/                              (QC audit scripts)
│   └── outputs/                              (12 CSVs: per-patient, blacklist, TSS, saturation, clinical, batch, …)
├── scripts/
│   ├── build_pan_malignant_report_v3.py
│   ├── build_per_celltype_report.py
│   ├── build_pan_myeloid_report.py
│   ├── per_patient_audit_v2.py               (post-CNV-fix)
│   └── per_patient_audit.py                  (original)
├── tests/                                     ← pytest harness; run `pytest tests/` from this dir
│   ├── test_schema.py                        (schema contract tests on parquet files)
│   ├── test_computations.py                  (unit tests for distance_category / Daigle Z / selectivity / aav_score)
│   └── test_filters.py                       (filter logic + pool size assertions)
└── reference/
    ├── refgene_hg38.bed                      (nearest-gene annotation source — GENCODE basic v45 TSSes)
    ├── encode_blacklist_hg38.v2.bed          (ENCODE blacklist used in QC #8)
    └── refGene.txt.gz                        (UCSC refGene source for refgene_hg38.bed)
```

## Reproducing from source

1. Clone pipeline: `git clone https://github.com/jpmcginnis/HumanAAVProject.git`
2. Install env: `conda env create -f atacseq_env.yml` (in `Glioblastoma snATACseq pipeline/`)
3. Point `config/paths.yml` at your local raw-fragments storage (originals from the 8 cohorts' accessions; see `cohort_attribution.csv`).
4. Run `snakemake --use-conda --cores 16` on an r6i.4xlarge-equivalent box (we used EC2 with a 500 GB EBS volume; takes ~24 h end-to-end).
5. Re-run the three report builders: `python scripts/build_pan_malignant_report_v3.py && python scripts/build_per_celltype_report.py && python scripts/build_pan_myeloid_report.py` (<10 min from parquets).

For a one-shot restore of the full compute environment, the EBS snapshot ID is `snap-0e28e2cef18747e24` (us-east-1). See `RESTORE.md` in the GitHub pipeline repo.

## Known cohort slug misnomer

The folder/column slug **`gse276177_khan_astro`** should be read as **Sojka et al. 2025 Nature Cell Biology**. The slug came from an early typo and was kept rather than regenerating the matrix. Use the real citation in any publication derived from this deposit. See `cohort_attribution.csv` row for `gse276177_khan_astro`.

## Known limits (summary — see FINDINGS.md for detail)

- **chr7 excluded** from pan-malignant clean pool: universal chr7+ CNV in GBM inflates strength independent of cell-type specificity.
- **hra004942** has 20-250 malignant cells per patient (vs 10K-80K in gbm_space). Low-cell replication contributes less information.
- **gse165037** is sci-ATAC (combinatorial indexing), not 10x. Peaks per cell are shallower; co-opening patterns can diverge. Replication here is orthogonal-technology evidence.
- **guilhamon** confirmed hg38 (pan-malignant pool matches CATLAS peaks across 4 patients × 500-700 malignant cells each).
- **surprises.parquet has 100% null in some columns** (`patient_id`, `delta`, `direction`, `n_low`, `n_high`, `frac_high`) — upstream artifact. Use populated columns only.
- **Batch metadata absent** — only patient_id/sample_id across all 8 cohorts. Batch effects collapse into patient-level random effects (handled by the Bayesian pass-2 model).
- **Per-patient TSS enrichment / FRiP not in processed data.** Rely on upstream paper QC thresholds; cross-cohort QC compared in methods_comparison.md.
- **CNV calling** uses the chr7+/chr10- ratio method (not full inferCNV); mirrored into `malignant_cnv` for all cohorts (fragments-mode fix in v1.0.0).

## Citation

```
McGinnis JP (2026). GBM Enhancer Atlas v1.0.0: cross-cohort scATAC-seq
derived AAV enhancer candidates. Zenodo. https://doi.org/<pending>
```

Please also cite the eight upstream cohorts (see `cohort_attribution.csv`).

## For the next Claude / future reader

- **Start here:** this file + `reports/pan_malignant_report_v3.html`.
- **Primary parquets for ad-hoc queries:** `data/matrix/pan_malignant_scores_v4.parquet` and `data/matrix/candidate_scores.parquet`. Both polars-queryable at ~10M rows/sec locally.
- **Full atlas (bigger):** `data/matrix/enhancer_candidate_matrix.parquet` (209M rows).
- **For a new cohort to add:** pipeline is in GitHub (`jpmcginnis/HumanAAVProject`); spin up EC2 from snapshot `snap-0e28e2cef18747e24` (see repo's `RESTORE.md`).
- **For a bug or refinement in scoring:** report builders are self-contained in `scripts/build_*.py`; re-run locally in the atacseq conda env.
- **Avoid pushing data to Google Drive from a script.** Drive File Stream in Mirror mode can trigger adjacent-folder re-materialization that fills the disk. Direct-upload via the Drive app or S3 instead.
- **This folder is the working copy.** The authoritative version at upload is on Zenodo; GitHub holds pipeline source. Avoid auto-sync between the two without explicit intent.
