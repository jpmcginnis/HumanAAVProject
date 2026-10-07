# Glioblastoma snATAC-seq Enhancer Atlas Pipeline

End-to-end pipeline that harmonizes every publicly accessible primary-tissue
IDH-wildtype glioblastoma snATAC-seq and 10x Multiome dataset into one
re-queryable **candidate enhancer matrix** — a parquet tensor keyed on
``(peak_id, cell_type, cohort, patient_id, region)``.

Built to drive AAV-delivered cell-type-specific enhancer selection for a BCM
Neurosurgery R01, but the matrix is general-purpose for any question that
asks "which GRCh38 regions are open, in which cells, in which patients, how
consistently?"

## Final atlas (as of 2026-10-07)

| Cohort | First author … last author | Journal (year) | Accession | Modality | Patients | Nuclei |
|---|---|---|---|---|---|---|
| GBM-Space | De Jong & Saraswat et al. | bioRxiv (2025) | E-MTAB-17183 | multiome | 12 | 1,041,869 |
| Khan — astrocyte trajectory | — (GEO-only) | GEO (2024) | GSE276177 | snATAC | 3 | 289,242 |
| TCGA scATAC (GBM subset) | Terekhanova … Corces/Chang | Science (2024), PMID 39236169 | TCGA-ATAC-Seq-2024 | snATAC | 9 | 157,685 |
| Wang — "under therapy" | Wang L … Diaz AA | Nature Cancer (2022), PMID 36539501 | GSE174554 | snATAC | 6 | 14,378 |
| Lin — TME atlas | Lin, Chen, Li, Chen, Fang … Qu | Nature Neuroscience (2026) | HRA004942 (NGDC) | snATAC | 11 | 6,735 |
| Wang — "phenotypes" | Wang L … Diaz AA | Cancer Discovery (2019), PMID 31554641 | GSE138794 | snATAC | 3 | 6,614 |
| Guilhamon | Guilhamon P … Lupien M | eLife (2021), PMID 33427645 | GSE139136 | snATAC | 4 | 3,818 |
| GSE165037 | — (GEO-only, no PMID linked) | GEO | GSE165037 | snATAC | 3 | 3,327 |
| **TOTAL** | | | | | **51 patients** | **1,523,668 nuclei** |

Three cohorts flagged but not included:
- **Wang 2024 Sci Adv (DNBelab C4)** — pending DNBelab-platform ingestion adapter.
- **spatial_epigenomic_niches (Kint/Gallo 2025 bioRxiv)** — no public deposit yet; email corresponding author.
- **Greenwald 2024 Cell** — pending dbGaP access.

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
├── reports/
│   └── templates/            # Jinja2 HTML report templates
├── notebooks/                # exploratory
├── tests/                    # smoke tests
├── CLAUDE.md                 # project brief (useful if resuming with an LLM collaborator)
└── DEPOSIT_LAYOUTS.md        # per-cohort GEO/dbGaP/Zenodo deposit layouts
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
