# Glioblastoma snATAC-seq Enhancer Atlas Pipeline

End-to-end pipeline that harmonizes every publicly accessible primary-tissue
IDH-wildtype glioblastoma snATAC-seq and 10x Multiome dataset into one
re-queryable **candidate enhancer matrix** — a parquet tensor keyed on
``(peak_id, cell_type, cohort, patient_id, region)``.

Built to drive AAV-delivered cell-type-specific enhancer selection for a BCM
Neurosurgery R01, but the matrix is general-purpose for any question that
asks "which GRCh38 regions are open, in which cells, in which patients, how
consistently?"

## Latest atlas release

> **Current deposit: [`analysis_2026-10-14/`](analysis_2026-10-14/) — v2026-10-11.**
> Zenodo DOI [10.5281/zenodo.23252324](https://doi.org/10.5281/zenodo.23252324).
> Previous Oct-7 draft is in git history (prior to this release); the Oct-7 draft had three issues reviewer-caught against the Oct-11 state — see "Changes vs. the Oct-7 draft" below and `analysis_2026-10-14/CHANGELOG.md` for the full audit trail.

## Final atlas (as of 2026-10-11, canonical first-author slugs)

| Cohort (slug) | First author (last author) | Journal (year) | Accession | Modality | Patients | Nuclei |
|---|---|---|---|---|---|---|
| `greenwald_2024` | Greenwald CA (Suva) | Cell (2024) | E-MTAB-17183 (ArrayExpress) | 10x snMultiome | 12 | 1,041,869 |
| `sojka_2025` | Sojka I (Goldman) | Nat Cell Biol (2025) | GSE276177 | 10x snATAC | 3 | 289,242 |
| `sundaram_2024` | Sundaram S (Andersen) | Cell Rep (2024) | dbGaP phs001287 (TCGA PanCanAtlas subset) | 10x snATAC | 8 (X008 excluded) | 157,685 |
| `zhao_2024` | Zhao et al. | Cell Rep Med (2024) | HRA004942 (CNCB-GSA) | snATAC | 11 | 6,735 |
| `mathewson_2021` | Mathewson KL (Lupien) + Wang 2024 | Cell (2021) / Sci Adv (2024) | GSE163654 + DNBelab | snATAC / DNBelab | 5 | 14,378 |
| `guo_2020` | Guo X | Nat Commun (2020) | GSE138794 | 10x snATAC | 3 | 6,614 |
| `guilhamon_2021` | Guilhamon P (Lupien) | eLife (2021) | GSE139860 | 10x snATAC | 4 | 3,818 |
| `nagaraja_2021` | Nagaraja S (Monje) | 2021 | GSE165037 | sci-ATAC | 2 | 3,327 |
| **TOTAL** | | | | | **48 patients** | **1,523,668 nuclei** |

### Changes vs. the Oct-7 draft

- **Greenwald now scored.** Oct-7 draft flagged Greenwald as "pending dbGaP access"; the Multiome snRNA label-transfer path removes that blocker, and Greenwald contributes 12 patients (Multiome donors) to the Oct-11 atlas.
- **`zhao_2024` patient count correction.** Oct-7 `cohort_attribution.csv` said 7 patients; the h5ad actually contains 11 unique donor IDs. 4-patient attribution error, corrected.
- **X008 explicit exclusion.** `sundaram_2024` lists 9 patients in deposit, 8 scored (X008 excluded due to a snapatac2 Rust panic on its obs structure). Noted in `cohort_attribution.csv` and `analysis_2026-10-14/README.md` Known Limitations.
- **`mathewson_2021`** and **`nagaraja_2021`** each lose 1 patient to the MIN_CELLS = 50 floor in the new patient-level pseudobulk pool rule (patient × label groups below 50 cells are dropped as too noisy to pseudobulk).
- **New scoring path.** `pilot_v5/score.py` replaces `src/scoring.py` as the active scorer: 1,958,753 union peaks (74% de-novo vs CATLAS) replaces the CATLAS-only 544k reference; arm-adjusted Daigle Z replaces genome-wide Z (removes chr7-gain/chr10-loss inflation of malignant peaks); patient-level replication replaces cohort-level; off-target per-non-target-cell-type contamination metrics added for AAV safety.
- **n_target_patients semantics fix.** Previously counted pseudobulk COLUMNS (patient × label); now counts unique (cohort, patient_id) tuples, bounded by atlas size. The old quantity preserved as `n_target_pseudobulks`.

Two cohorts remain flagged but not included:
- **Wang 2024 Sci Adv (DNBelab C4)** — pending DNBelab-platform ingestion adapter.
- **spatial_epigenomic_niches (Kint/Gallo 2025 bioRxiv)** — no public deposit yet; email corresponding author.

## What's in this repo

```
├── Snakefile                 # orchestrator — fetch → ingest → QC/CNV/labels → quantify → matrix → scoring → rank → reports
├── config/
│   ├── datasets.yaml         # cohort catalog (slug, accession, URL, patient count, platform, access tier)
│   ├── pipeline.yaml         # thresholds + scoring formula
│   ├── celltypes.yaml        # canonical cell type ontology + GBmap / CATLAS mappings
│   └── marker_tsses.tsv      # curated marker-gene TSSes for cell type scoring (CATLAS-space)
├── src/
│   ├── fetch.py              # per-cohort downloaders (GEO, ArrayExpress, Zenodo, GDC)
│   ├── ingest.py             # fragments → per-sample snapatac2 h5ads (fragment-mode cohorts)
│   ├── ingest_matrix.py      # already-quantified MTX / per-sample peaks → CATLAS-space h5ad (matrix-mode)
│   ├── flatten_cohort.py     # per-sample snapatac2 h5ads → cohort-level CATLAS-space h5ad (fragment-mode terminal)
│   ├── qc.py                 # TSS enrichment, fragment-size, nuclei filters (matrix-mode only; fragment-mode uses flatten_cohort)
│   ├── cnv_malignant.py      # chr7+/chr10- ratio → malignant_cnv obs column
│   ├── label_transfer.py     # GBmap scRNA + CATLAS snATAC dual label transfer (matrix-mode)
│   ├── marker_peak_scoring.py # fallback: cell type labels from marker-gene TSS ± 2kb accessibility
│   ├── peak_quantify.py      # CATLAS peak lift (matrix-mode; fragment-mode handled by flatten_cohort)
│   ├── matrix_build.py       # cohort h5ads → long-format (peak × cell_type × cohort × patient) parquet tensor
│   ├── pan_malignant_matrix.py  # same, but pools ALL tumor cells (CNV flag OR malignant_unresolved label) into one bucket — bypasses the Neftel-subtype gap and the marker-peak mis-labeling of CNV-malignant cells in gbm_space
│   ├── scoring.py            # Pass 1 polars point-estimates + Pass 2 PyMC hierarchical ADVI posteriors
│   ├── ranking.py            # composite ranking + 200-2000 bp + CNV-safe constraints → top candidates per cell type
│   ├── surprises.py          # patient-subpopulation, multiome discordance, cross-dataset outlier detectors
│   ├── dataset_issues.py     # cohort coverage audit → verdict on whether existing data suffices
│   ├── report.py             # Jinja2-based HTML renderer
│   └── extract_rds.R         # R-side ArchR/Signac/SnapATAC .rds → h5ad converter (controlled-access cohorts)
├── infra/
│   ├── atacseq_env.yml       # conda env spec
│   ├── aws_spinup.sh         # EC2 r6i.4xlarge + 500GB gp3 attach + conda env setup
│   ├── aws_shutdown.sh       # clean stop, keeps EBS
│   ├── aws_restore.sh        # restore from EBS snapshot
│   └── watchdog.sh           # 20-min cron: restarts snakemake up to 3x on crash
├── pilot_v5/                 # Oct-11 scoring revisions (supersedes src/scoring.py for the current deposit — see pilot_v5/README.md)
│   ├── requant{,_greenwald,_greenwald_parallel}.py  # pseudobulk builders against the union-peak reference
│   ├── greenwald_finish.py   # 8-worker finisher for the 10.5B-nnz Greenwald h5ad
│   ├── remap_greenwald_coarse.py  # Greenwald's 99 fine labels → 14 coarse schema (keeps 4 Neftel malignant states)
│   ├── cnv_patch.py          # hg38 arm coords + arm-conditional Z helpers + chr7-elevated flag
│   └── score.py              # per-report scorer: arm-adjusted Daigle Z + CNV + off-target metrics
├── analysis_2026-10-14/      # the Oct-11 atlas deposit (full README, SCHEMA, CHANGELOG, polished scripts, tests, hg38 reference, cohort attribution, Zenodo deposit metadata)
├── CLAUDE.md                 # project brief (useful if resuming with an LLM collaborator)
├── DEPOSIT_LAYOUTS.md        # per-cohort GEO/dbGaP/Zenodo deposit layouts
├── RESTORE.md                # EBS-restore quickstart
└── requirements.txt          # pip deps (used alongside the conda env for a few extras)
```

## Running the pipeline

### One-time setup

```bash
# AWS (or any Linux box with ~128 GB RAM, 500 GB SSD):
./infra/aws_spinup.sh
conda env create -f infra/atacseq_env.yml && conda activate atacseq
```

### Full pipeline

```bash
# Fetch references + all cohorts, then build atlas.
snakemake --cores all

# Just rebuild the matrix after adding a cohort (fragment-mode cohorts chain
# through flatten_cohort.py, not through qc/cnv/labels/peak_quantify).
snakemake matrix_build --cores 32

# Rank candidates with a different scoring formula (edit config/pipeline.yaml first).
snakemake rank_top_n --forcerun

# Generate reports.
snakemake report
```

### Querying the matrix (the main product)

```python
import polars as pl
m = pl.read_parquet("matrix/enhancer_candidate_matrix.parquet")

# Top TAM-specific enhancers that replicate in ≥3 cohorts
top_tam = (
    m.filter(pl.col("cell_type") == "TAM")
     .group_by("peak_id")
     .agg([
         pl.col("frac_accessible").mean().alias("mean_frac"),
         pl.col("cohort").n_unique().alias("n_cohorts"),
         pl.col("patient_id").n_unique().alias("n_patients"),
     ])
     .filter(pl.col("n_cohorts") >= 3)
     .sort("mean_frac", descending=True)
     .head(20)
)

# TAM-selective (open in TAM, closed in homeostatic microglia)
tam_pivot = m.pivot(
    index="peak_id", on="cell_type",
    values="frac_accessible", aggregate_function="mean",
)
selective = tam_pivot.filter(
    (pl.col("TAM") > 0.3) & (pl.col("microglia") < 0.05)
)
```

Full schema in `matrix/MATRIX_SCHEMA.md`.

## Design decisions worth knowing

**Common peak coordinate system = CATLAS 544,735 peaks.** Every cohort is
projected onto the same peak atlas (Li et al. 2023 Science, GSE244618) before
any cross-cohort comparison. This resolves the peak-coordinate-mismatch problem
that killed 7/8 TCGA-only candidates in prior work with ad-hoc coordinate
systems.

**CNV-based malignant calling** (chr7+/chr10- ratio ≥ 2.0) runs before any
cell type label transfer. Malignant cells go through a separate Neftel-state
classifier (NPC-like / OPC-like / AC-like / MES-like); non-malignant cells
through TME label transfer. Both labels are kept with explicit provenance per
cell.

**Label transfer:** TME labels from GBmap (Ruiz-Moreno 2025, 330K scRNA cells
over 109 patients) via gene-activity bridge; neuronal subtype labels from
CATLAS (107 cell types, snATAC-native). Fragment-mode cohorts without
gene-activity data fall back to marker-peak scoring at ±2 kb of a curated
marker-gene TSS list (`config/marker_tsses.tsv`).

**Scoring is two-pass.** Pass 1 is point-estimate (strength / selectivity /
consistency / n_cohorts) across all 5.6 M (peak × cell_type) pairs via polars.
Pass 2 fits a single vectorized PyMC hierarchical Bayesian model over the top
2000-per-cell-type shortlist (22K pairs) with ADVI — per-cell-type logit
baselines, per-pair deviations with partial pooling, per-(pair × cohort)
offsets, Binomial likelihood against aggregated per-cohort counts. 30k ADVI
iterations ≈ 15-20 min on a 16-core box.

**Ranking applies real-world constraints:** 200-2000 bp length (AAV ITR
budget), CNV-safe (not in chronic chr7+/chr10- regions), PmlI counter-selection
site-free. The `top_candidates.parquet` has 10 ranked candidates per cell type,
ready to order from Genewiz.

## Reproducibility caveats

Three cohorts' labels were produced by **marker-peak scoring** on
`config/marker_tsses.tsv` because they lack coordinated gene-activity (snATAC
only, no paired RNA). Changing the marker list changes labels. Keep the TSV
versioned.

**Pan-malignant mis-labeling.** Marker-peak scoring routinely assigns
CNV-called malignant cells to TME types (TAM, microglia, neuron) when the
tumor cell's chromatin incidentally opens at a marker-gene promoter. In
GBM-Space this hit 438,994 cells (42% of the cohort) — ALL CNV-malignant,
NONE labeled as `malignant_unresolved`. For questions about pan-tumor
accessibility, `src/pan_malignant_matrix.py` rebuilds the aggregate using
``(malignant_cnv == 1) OR (cell_type == "malignant_unresolved")`` as the
inclusion rule. This recovers ~2.7× more pan-malignant cells and takes
cross-cohort replication from max 3 cohorts to 6 cohorts for the top hits.

The 4 TCGA GBM samples behind dbGaP phs000178 are **not** in the public
pipeline run. If you have dbGaP access, drop their BAMs → fragments → under
`/data/raw/tcga_scatac/Cancer_scATACseq_data/` and re-run `snakemake
flatten_cohort --forcerun`.

## Data deposition

The pipeline produces artifacts suitable for:
- **Zenodo**: `enhancer_candidate_matrix.parquet` + `candidate_scores.parquet`
  + `top_candidates.parquet` + `MATRIX_SCHEMA.md` + the four HTML reports.
  One zip → one DOI → citable in the R01.
- **CellxGene**: per-cohort `catlas_quantified.h5ad` (scverse-standard,
  cell-level annotations). Curation takes weeks; longer-term track.

## License

MIT — see LICENSE file.

## Citation

If you use this pipeline or the candidate matrix, please cite the eight
underlying datasets (table above) *and* this repo. A proper DOI will follow
the Zenodo deposit of the final matrix.
