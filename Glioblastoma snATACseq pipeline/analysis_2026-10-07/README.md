# GBM Enhancer Atlas — 2026-10-07 (v4)

**Project:** AAV-delivered cell-type-specific enhancer discovery for GBM. Backs an R01 due June 2026.
**Owner:** JP McGinnis (jpmcginnis1@gmail.com), BCM Neurosurgery.
**Pipeline repo:** <https://github.com/jpmcginnis/HumanAAVProject> → `Glioblastoma snATACseq pipeline/`
**Atlas scope:** 8 cohorts, 51 patients (45 with ≥10 pan-malignant cells), 1,523,668 nuclei, 544,735 CATLAS peaks.

---

## Where to start if you are a new Claude session

1. **Open `reports/pan_malignant_report_v3.html`** — current primary deliverable. Top-50 AAV cloning shortlist for pan-malignant targeting, with full methods caveats.
2. **Open `reports/per_celltype_report.html`** — same workflow but across 10 cell types (pan-malignant section at top, then microglia, TAM, astrocyte, oligodendrocyte, OPC, neuron, GABA_neuron, endothelial, T_cell).
3. **Open `reports/pan_myeloid_report.html`** — TAM + microglia merged for pan-myeloid AAV targeting (JP's call: no therapeutic reason to hit one without the other).
4. **Read `qc/FINDINGS.md`** — synopsis of 11 QC checks, all triaged or resolved. No blocking issues.
5. **Read `reports/methods_comparison.md`** — per-cohort snATAC methods extracted from GEO/EBI/papers. Flags for low-cell cohorts (hra004942), orthogonal tech (gse165037 sci-ATAC), and cohort slug misnomers (gse276177 is Sojka et al., not Khan).

## Current state (what's deliverable now — 2026-10-08 refresh)

- **Pan-malignant pool:** 575,971 cells across 45 patients × 8 cohorts (CNV | label union). Pre-CNV-fix (v3 and earlier) was 33 patients due to a bug in `src/flatten_cohort.py` that silently defaulted `malignant_cnv=0` for fragments-mode cohorts (tcga_scatac + gse276177). Fixed 2026-10-07, patch is in the GitHub repo.
- **Distance filter:** |dist to TSS| ≥ 2 kb floor (excludes promoter-proximal non-enhancers). **No upper cap** — for AAV the enhancer is extracted from genomic context and placed next to a minimal promoter, so native-genome distance is irrelevant. `distance_category` column (near / distal / far-distal / gene-desert) is informational, not a filter.
- **Daigle Z-score:** every peak carries `z_daigle = (strength_target - mean_strength_other) / sd_strength_other` and a `passes_daigle_z2` boolean. Z ≥ 2 is the Allen Institute Armamentarium criterion for cell-type specificity. **11,641 peaks** pass Z ≥ 2 atlas-wide for pan-malignant; subset within the AAV-clean pool is highlighted in each report.
- **Candidate AAV peaks (pan-malignant, |TSS|≥2kb + sel≥2× + non-chr7 + non-housekeeping):** **617 clean candidates** (up from 461 after dropping the 100 kb cap). Top-ranked near CPNE4, ADAMTSL1, ZFP36L1, GPNMB (4.1× selectivity), LINC01235.
- **Per-cell-type candidates:** 690 for pan-malignant, 684 for neuron, 406 for endothelial, 373 for TAM, 358 for GABA_neuron, 346 for astrocyte, 208 for OPC, 161 for microglia, 108 for T_cell, 57 for oligodendrocyte. Strict mode (sel ≥ 1.0×) holds for OPC + neuron + pan-malignant; fallback mode for the glia-confounded types ranks by selectivity × composite.
- **Pan-myeloid:** 380 peaks accessible in both TAM AND microglia (≥5% strength each; Table B is primary). Strict selectivity table (Table A, sel ≥ 1.0×) has 10 peaks because reactive glia confound the non-myeloid baseline.

## Directory contents

```
.
├── README.md                                 ← you are here
├── data/
│   ├── matrix/
│   │   ├── enhancer_candidate_matrix.parquet (3.5 GB, 209M rows, peak × patient × cohort × cell_type)
│   │   ├── candidate_scores.parquet          (134 MB, per-peak × per-cell-type scores)
│   │   ├── surprises.parquet                 (1.1 GB, pattern-detector outputs — some columns 100% null, see FINDINGS #15)
│   │   ├── pan_malignant_matrix_v4.parquet   (188 MB, pan-malignant long-form, 45 patients)
│   │   ├── pan_malignant_scores_v4.parquet   (14 MB, per-peak pan-mal scores — primary for v3 report)
│   │   ├── panmal_v4_run.log
│   │   ├── top_candidates.parquet
│   │   ├── provenance.parquet
│   │   ├── dataset_issues.parquet
│   │   └── MATRIX_SCHEMA.md                  ← start here for schema questions
│   └── processed_h5ads/                      (69 GB, per-cohort CATLAS-quantified h5ads with CNV + labels)
│       ├── gbm_space/                        (62 GB — the whale; 1.04M cells, 12 patients)
│       ├── gse276177_khan_astro/             (2.9 GB; Sojka et al. 2025, 3 patients)
│       ├── sundaram_gbm/                     (2.1 GB; TCGA scATAC, 9 patients)
│       ├── gbm_tme_atlas_hra004942/          (1.1 GB; 7 patients, low cells per patient — see QC caveats)
│       ├── mathewson_lupien/                 (260 MB; 5 patients)
│       ├── gse138794_guo/                    (270 MB; 3 patients)
│       ├── guilhamon/                        (140 MB; 4 patients, hg38 verified)
│       └── gse165037/                        (80 MB; 2 patients, sci-ATAC not 10x)
├── reports/
│   ├── pan_malignant_report_v3.html          ← PRIMARY report (45 patients, methods caveats)
│   ├── per_celltype_report.html              (10 cell types incl. pan-malignant)
│   ├── pan_myeloid_report.html               (TAM+microglia merged)
│   ├── methods_comparison.md                 (per-cohort snATAC methods appendix)
│   ├── pan_malignant_top50_distal_selective_v3.csv
│   ├── pan_malignant_top100_all_v3.csv
│   ├── pan_malignant_top30_distal_selective_v3.bed
│   ├── pan_myeloid_top50_shared.csv          (primary pan-myeloid output)
│   ├── pan_myeloid_top30_shared.bed
│   ├── pan_myeloid_top50_distal_selective.csv (strict/selective)
│   ├── pan_myeloid_top30.bed
│   ├── per_celltype/                         (per-cell-type top-50 CSVs + top-30 BEDs)
│   ├── per_patient_audit.csv / .html
├── qc/
│   ├── FINDINGS.md                           ← synopsis of 11 QC checks
│   ├── scripts/                              (QC scripts: per-patient audit, blacklist intersect, TSS reclass)
│   └── outputs/                              (12 CSVs — per-patient, CN flags, batch metadata, blacklist, TSS, saturation, clinical, surprises-null)
└── scripts/
    ├── build_pan_malignant_report_v3.py
    ├── build_per_celltype_report.py
    ├── build_pan_myeloid_report.py
    └── per_patient_audit.py
```

## Where everything else lives (not in this folder)

- **S3 archive (authoritative, encrypted at rest):** `s3://jpm-atacseq-archive-2026/final_atlas_2026_10_07/`
- **GitHub (pipeline source):** <https://github.com/jpmcginnis/HumanAAVProject> → `Glioblastoma snATACseq pipeline/`
- **EBS snapshot (full compute env for restore):** `snap-0e28e2cef18747e24` (us-east-1, 2026-10-08). Restore procedure in the GitHub repo's `RESTORE.md`. Scheduled check 2026-11-09 to move to Archive tier (~75% cheaper storage).
- **EC2:** currently spun down. No running instance.
- **Google Drive mirror** (previously used as working copy; JP is moving away from it due to Mirror-mode sync surprises that wipe the Desktop folder).

## Known caveats (short list — see methods_comparison.md and FINDINGS.md for detail)

- **chr7 excluded from pan-malignant clean pool.** Universal chr7+ CNV in GBM inflates strength independent of cell-type specificity.
- **hra004942** has 20-250 malignant cells per patient (vs 10K-80K in gbm_space). Low-cell replication contributes less information. Suspect low-throughput platform.
- **gse165037** is sci-ATAC (combinatorial indexing), not 10x. Peaks per cell are shallower; co-opening patterns can diverge from droplet. Replication here is orthogonal technology — extra-strong evidence when it agrees.
- **gse276177 slug is misnomer.** Actually Sojka et al. 2025 Nature Cell Biol, not Khan. Data content and analysis are correct; just the folder name.
- **guilhamon confirmed hg38** (pan-malignant pool matches CATLAS peaks across 4 patients × 500-700 malignant cells each).
- **surprises.parquet has 100% null in some columns** (`patient_id`, `delta`, `direction`, `n_low`, `n_high`, `frac_high`) — upstream artifact. Use the populated columns (peak_id, detector, cell_type, pattern) only. Regenerate on next scoring pass.
- **batch metadata absent.** Only patient_id/sample_id across all 8 cohorts. Batch effects collapse into patient-level random effects (handled by the Bayesian pass2).
- **Per-patient TSS enrichment / FRiP not in processed data.** Relying on upstream paper QC thresholds; cross-cohort QC compared in methods_comparison.md.

## Open to-do

- (Non-blocking) Pilot TCGA scATAC end-to-end before scaling further.
- (Non-blocking) Resolve TCGA scATAC — pull from GDC manifest.
- (Non-blocking) Wang Sci Adv DNBelab adapter if we want to add that cohort.
- JP to email Marco Gallo for spatial_epigenomic data (external dependency).
- R-side RDS extraction (ArchR/Signac/SnapATAC → h5ad) — in progress.

---

## For the next Claude session

- **Start by reading this file.** Then the three HTML reports in order: pan_malignant_report_v3, per_celltype_report, pan_myeloid_report.
- **Primary matrix file for ad-hoc queries:** `data/matrix/pan_malignant_scores_v4.parquet` + `data/matrix/candidate_scores.parquet`. Both polars-queryable at ~10M rows/sec locally.
- **For a new cohort to add:** the pipeline is in GitHub HumanAAVProject. Spin up EC2 from snapshot `snap-0e28e2cef18747e24` (see repo RESTORE.md).
- **For a bug or refinement in the scoring:** report builders are self-contained in `scripts/build_*.py` and read parquets from `data/matrix/`. Re-run locally in `atacseq` conda env (`source /usr/local/Caskroom/miniconda/base/etc/profile.d/conda.sh && conda activate atacseq`).
- **Avoid pushing data to Google Drive from a script** — JP's experience is that Drive File Stream in Mirror mode can trigger adjacent-folder re-materialization that fills the disk. Direct-upload via the Drive app or S3 instead.
- **This folder is the working copy for JP.** The S3 archive is the authoritative version; GitHub holds pipeline source. Avoid auto-sync'ing between the two without JP's explicit say-so.
