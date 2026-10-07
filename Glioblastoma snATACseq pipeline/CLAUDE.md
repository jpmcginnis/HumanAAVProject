# GBM Enhancer Atlas — Claude Code project brief

**Owner:** JP McGinnis, MD, PhD (jpmcginnis1@gmail.com) — BCM Neurosurgery
**Purpose:** Harmonize every publicly accessible primary-tissue IDH-WT GBM snATAC-seq / 10x Multiome dataset into one re-queryable candidate enhancer matrix, so JP can mine it for AAV-delivery enhancer questions time and again.
**Deadline context:** June NIH R01 — the output of this pipeline decides whether the grant is written for *validation* (if existing data is sufficient) or *discovery* (if we need our own n=24 Multiome cohort).

---

## What success looks like

1. **A persistent candidate enhancer matrix** — a parquet tensor keyed on `(peak_id, cell_type, cohort, patient)` with accessibility, specificity, consistency, and (where available) coordinated-RNA-expression scores. This is the artifact JP returns to for every future GBM enhancer question.
2. **A per-cell-type report** naming the enhancers best positioned to drive cell-type-specific expression — strong, consistent across patients, and selective against other cell types.
3. **A "clone these first" list** — 5–10 enhancer candidates ranked by composite score, with coordinates ready to order from Genewiz.
4. **A "surprises" section** — findings that weren't the question being asked but are worth knowing (unexpected cell-type co-accessibility, cross-dataset discrepancies, novel candidate regions).
5. **A "do we need our own data?" verdict** — explicit statement on whether existing public data suffices for Aim 3, or whether JP should write the R01 for prospective n=24 Multiome generation.

---

## Known data discrepancies (2026-10-04) — resolve before first full run

JP flagged four issues on initial review of `config/datasets.yaml`; a second-pass
literature sweep then surfaced three more datasets that had been missed. Treat
the config as provisional until all seven items are closed out.

**Datasets added in the second-pass sweep (not yet verified end-to-end):**

- **`gbm_tme_atlas_hra004942`** — Nature Neuroscience 2026 paper. 100 primary
  GBM patients, 121 samples across spatial + scRNA + scATAC + patch-seq. This
  is the largest primary GBM cohort found. Deposited at NGDC GSA (Chinese
  archive) under HRA004942; processed data on Zenodo. **IDH status, ATAC-sample
  fraction, and NGDC-access-from-US-institutions are all unknowns to resolve
  before committing to this dataset.**
- **`spatial_epigenomic_niches`** — bioRxiv 2025 (DOI 10.1101/2025.05.09.653178).
  24 patients, 28 primary IDH-WT GBM, mostly spatial ATAC (28) with small
  snATAC (3) and Multiome (2) subsets. No GEO/dbGaP accession found yet; check
  Data Availability section directly.
- **`wang_sciadv_tc_ptb`** — Wang 2024 Science Advances. DNBelab C4 platform
  (not 10x), paired tumor-core vs peritumoral-brain samples. Needs a DNBelab
  ingestion handler. Accession likely CNGB or CNSA; verify.

**Original four discrepancies:**

1. **GSE165037 citation is wrong.** Was labeled "Raviram et al. 2023 PNAS" but
   that paper is scHi-C/Paired-seq, not snATAC. The actual GSE165037 title is
   "Single Cell Analysis of Chromatin Accessibility Reveals Genetic and
   Regulatory Heterogeneity in Glioblastomas" — a Guilhamon-lab-adjacent
   deposit. Cohort slug renamed from `raviram` → `gse165037`. **Dedup against
   Guilhamon GSE139136 (4 patients) is required** — sample IDs look distinct
   (`GBM4` vs `GBM_4349`) but computational CNV-fingerprint dedup is the only
   way to prove non-overlap. Previous "RL3 multi-region case (4 sections)"
   attribution was incorrect and has been removed; no such samples exist in
   GSE165037.

2. **TCGA scATAC S3 contents unverified.** The `integrated.h5ad` at
   `s3://jpm-atacseq-archive-2026/gbm_pipeline_complete/tcga_scatac/` is
   assumed from a 2026-04-23 memory note that confirmed the folder exists
   with 6 samples / 22,004 cells, but not the specific AnnData file. The
   fetch step now `aws s3 sync`s and then checks for `integrated.h5ad`. If
   missing, it falls back to looking for `.arrow` files; if those are also
   absent, it fails fast. If only `.arrow` files exist, we need an
   ArchR → AnnData conversion step (not yet implemented).

3. **GSE174554 (Mathewson/Lupien) IDH-WT subset unverified.** The deposit
   has 113 samples total, 10 of which are snATAC, but IDH status and
   primary-vs-recurrent are not in sample titles. Table S1 of Mathewson 2022
   Nature Cancer must be downloaded and parsed to filter the 10 snATAC
   samples. The claim that "4 primary IDH-WT from Mathewson 2021 Nature are
   included as reanalysis" is also unverified — those 4 may live under a
   separate GEO accession. Ingest handler reads
   `mathewson_tableS1_IDH_status.tsv` from the raw directory; if absent,
   `idh_status='unknown'` propagates and matrix_build will exclude them.

4. **Expected cohort counts are soft.** Prior claim of ~23 IDH-WT primary
   snATAC samples across 4 existing cohorts assumed Mathewson 2022 contributed
   all 10 scATAC samples as IDH-WT primary. If Table S1 shows fewer qualifying
   samples, this number drops and the "do we need our own data?" verdict in
   `src/dataset_issues.py` becomes more likely to recommend prospective
   generation. Reflect any revised counts in the memory file
   `/projects/.../overview.md` after the first real ingest run.

---

## Hard rules

- **Primary human adult GBM tissue only.** No GSCs, organoids, PDX, cell lines, mouse, or FACS-sorted fractions. The AAV sees everything in the head; the reference must too.
- **IDH-wildtype only** (GBM 2021 WHO definition).
- **Common peak coordinate system = CATLAS 544K peaks** (resolves the peak-coordinate-mismatch problem that killed 7/8 TCGA-only candidates in prior work).
- **CNV-based malignant calling** (chr7+/chr10- ratio ≥ 2.0) before any cell-type label transfer; malignant cells go through a separate Neftel-state classifier, not the TME label transfer.
- **Label transfer:** TME labels from GBmap (disease-matched scRNA, 109 patients, 330K cells). Neuronal subtype labels from CATLAS (snATAC-native, 107 cell types). Report both with explicit provenance per cell.
- **For Multiome cohorts**, require both accessibility AND coordinated target-gene expression in the same cell as a stricter replication filter.

---

## Code conventions (JP's standing preferences)

- **Full replacement files only.** No sed/grep/one-liners; when fixing a script, print the complete runnable file.
- **Verbose progress output.** Every analysis script prints cell counts, filter-pass counts, per-cell-type stats, and top hits as computed. Never suppress.
- **Tufte visualization defaults.** Minimal chartjunk, data-ink maximized, no 3D, no unnecessary gridlines. HTML/Playwright-rendered widgets preferred over matplotlib for presentation figures.
- **Scripts run via** `python -u script.py 2>&1 | tee logfile.txt`.
- **Results CSVs** go to `/data/projects/atacseq/results/`. Downloaded outputs to local Mac at `/Users/jpmcginnis1/atacseq_project/results/`.
- **Long-running jobs run in tmux session named `work`.**

---

## Repo layout

```
gbm-enhancer-atlas/
├── CLAUDE.md                    # this file — read first every session
├── README.md                    # for humans
├── infra/
│   ├── aws_spinup.sh            # start EC2, mount EBS, activate conda, attach tmux
│   ├── aws_shutdown.sh          # stop EC2 cleanly
│   └── s3_layout.md             # bucket organization
├── config/
│   ├── datasets.yaml            # filtered primary-tissue cohort catalog
│   ├── pipeline.yaml            # thresholds, references, scoring formula
│   └── celltypes.yaml           # canonical cell-type ontology + GBmap/CATLAS mapping
├── Snakefile                    # orchestrator
├── src/
│   ├── __init__.py
│   ├── fetch.py                 # per-cohort download logic
│   ├── ingest.py                # cohort → standardized AnnData
│   ├── qc.py                    # TSS enrichment, fragment size, nuclei filtering
│   ├── cnv_malignant.py         # chr7+/chr10- malignant calling
│   ├── label_transfer.py        # GBmap + CATLAS dual label transfer
│   ├── peak_quantify.py         # quantify all cohorts against CATLAS peak atlas
│   ├── matrix_build.py          # assemble the persistent (peak × celltype × cohort × patient) tensor
│   ├── scoring.py               # consistency, selectivity, specificity scoring
│   ├── ranking.py               # composite ranking → top-N for cloning
│   ├── surprises.py             # anomaly / unexpected-pattern detection
│   ├── dataset_issues.py        # per-dataset QC issue log → "do we need our own data?" verdict
│   └── report.py                # HTML report rendering
├── reports/
│   └── templates/
│       ├── base.html
│       ├── celltype_report.html
│       ├── top_candidates.html
│       ├── surprises.html
│       └── dataset_issues.html
├── matrix/                      # persistent artifacts (gitignored, synced to S3)
│   ├── enhancer_candidate_matrix.parquet
│   ├── candidate_scores.parquet
│   ├── provenance.parquet
│   └── MATRIX_SCHEMA.md         # how to query the matrix
├── notebooks/                   # exploratory, versioned
└── tests/                       # schema and smoke tests
```

---

## AWS environment

- EC2 `i-000dddf603f4bde12`, us-east-1a, security group `sg-0ab95934ea158c30a`
- SSH key: `~/atacseq_project/bioinfo-key.pem`, user `ubuntu`
- Conda env `atacseq` on `/data/miniconda3`
- S3 archive bucket: `jpm-atacseq-archive-2026`
- `infra/aws_spinup.sh` handles the full restart. Use that; don't run the steps by hand.

---

## Running the pipeline

```bash
# One-time setup
./infra/aws_spinup.sh

# Full pipeline (all cohorts, all stages)
snakemake --cores all

# Just rebuild the matrix after adding a cohort
snakemake matrix_build --cores 32

# Rank candidates with a different scoring formula (edit config/pipeline.yaml first)
snakemake rank_top_n score_candidates --forcerun

# Generate reports
snakemake report
```

---

## Querying the persistent matrix later

```python
import polars as pl
m = pl.read_parquet("matrix/enhancer_candidate_matrix.parquet")

# "Show me the top TAM-specific enhancers that replicate in >=3 cohorts"
tam = (m.filter(pl.col("cell_type") == "TAM")
         .filter(pl.col("n_cohorts_reproducing") >= 3)
         .sort("composite_score", descending=True)
         .head(20))
```

See `matrix/MATRIX_SCHEMA.md` for the full column list and example queries.
