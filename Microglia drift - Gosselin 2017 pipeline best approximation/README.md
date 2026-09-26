# gosselin-drift-pipeline

A reproducible pipeline for measuring how much a human brain microglia
population has drifted from its ex-vivo transcriptome, benchmarked
against Gosselin et al. 2017 (*Science* eaal3222).

Given paired snRNA-seq samples from a "fresh" (D0) and a "cultured"
(later timepoint) condition, this pipeline reports the fraction of
Gosselin's 477-gene conserved microglia identity signature that dropped
> 2-fold at FDR < 0.05, using limma-voom exactly as described in
Gosselin's supplementary methods. The reference number is Gosselin's
**33 %** at day 7 in dissociated culture.

The pipeline is **validated against Gosselin's own data** — see
`validation/`. On Gosselin's paired ex-vivo → 7-day in-vitro human
microglia samples, this pipeline reproduces his primary published
headline within 1 percentage point (32.1 % on his data vs his
published 33 %) and reproduces his whole-transcriptome DOWN count to
within 7 genes out of 3,702.

---

## What the pipeline does, step by step

The pipeline is limma-voom applied to per-donor pseudobulk counts.
Full walkthrough in [`docs/PIPELINE_STEPS.md`](docs/PIPELINE_STEPS.md).

**Inputs**
1. A **gene × sample integer count matrix** — pseudobulk of microglia
   from each (donor × timepoint), or bulk RNA-seq counts of sorted
   microglia. Rows = HGNC gene symbols, columns = samples.
2. A **sample metadata table** — `donor` and `condition` for each
   column. `condition` should be `ex` for the fresh/day-0 timepoint
   and `iv` for the cultured timepoint.

**Steps** (all live in `pipeline/02_limma_voom.R`)
1. `filterByExpr(y, design = ~donor + condition)` — drop
   near-undetected genes.
2. `calcNormFactors(y, method = "TMM")` — TMM library-size norm.
3. `voom(y, design)` — log2-CPM + per-observation precision weights.
4. `lmFit(v, design = ~donor + condition)` — paired linear model
   (donor as blocking factor, condition as effect).
5. `eBayes(fit)` — empirical-Bayes variance moderation.
6. `topTable(fit, coef = "conditioniv", …)` — per-gene log2FC and
   BH-FDR.
7. Threshold: **FDR < 0.05 AND |log2FC| > 1** (Gosselin's stated
   cutoff).
8. Intersect the significantly-DOWN gene set with Gosselin's
   477-conserved signature to get the drift headline percentage.

**Outputs**
- Per-gene DE table (log2FC, p, adjusted p, expression) as CSV.
- Whole-transcriptome DE counts at > 2-fold and > 10-fold, both
  directions.
- The drift headline: fraction of the 477-conserved signature
  dropping > 2-fold at FDR < 0.05.
- Comparison bar chart vs Gosselin's 33 %.

---

## Repo layout

```
gosselin-drift-pipeline/
├── README.md                            # this file
├── requirements-python.txt              # Python deps
├── requirements-R.txt                   # R deps and how to install
├── pipeline/
│   ├── 01_h5_to_pseudobulk.py           # snRNA-seq h5 → pseudobulk counts
│   ├── 02_limma_voom.R                  # pseudobulk → DE stats
│   ├── 03_drift_summary.py              # DE stats → headline + figures
│   └── config_example.yaml
├── references/
│   ├── README.md                        # provenance of gene lists
│   ├── gosselin_conserved_signature_477.csv
│   └── gosselin_881_signature.csv
├── validation/
│   ├── README.md                        # validation writeup
│   ├── dump_gosselin_counts.py          # Table S2 → CSV
│   ├── run_validation.R                 # pipeline on Gosselin's own data
│   └── expected_numbers.md              # published targets to reproduce
└── docs/
    └── PIPELINE_STEPS.md                # detailed step-by-step walkthrough
```

---

## Quick start

```bash
# One-time setup (Python 3.11 recommended)
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-python.txt

# R dependencies (R >= 4.0 recommended; needs limma, edgeR from Bioconductor)
Rscript -e 'install.packages("BiocManager"); \
            BiocManager::install(c("limma","edgeR"), update=FALSE, ask=FALSE)'
```

```bash
# Step 1: snRNA-seq h5 files → pseudobulk counts
python pipeline/01_h5_to_pseudobulk.py \
    --config pipeline/config_example.yaml \
    --outdir output/

# Step 2: pseudobulk counts → per-gene DE stats
Rscript pipeline/02_limma_voom.R \
    --counts output/pseudobulk_counts.csv \
    --metadata output/sample_metadata.csv \
    --outdir output/

# Step 3: DE stats → drift headline + figures
python pipeline/03_drift_summary.py \
    --de output/limma_voom_de.csv \
    --signature references/gosselin_conserved_signature_477.csv \
    --outdir output/
```

---

## Validation

The validation folder reruns this pipeline on Gosselin's published
human RNA-seq data (Table S2, paired ex-vivo → 7-day in-vitro
microglia). It reports the numbers our pipeline gets vs his published
numbers, so any user can confirm the pipeline is faithful before
applying it to their own samples.

Head-to-head on Gosselin's own N=5 paired data (see
`validation/expected_numbers.md`):

| metric | Gosselin published | this pipeline | delta |
|---|---:|---:|---|
| WT DOWN >2-fold, FDR<0.05 | 3,702 | 3,709 | +7 (99.8 %) |
| WT UP >2-fold, FDR<0.05 | 2,263 | 2,913 | +29 % (see caveats) |
| WT DOWN >10-fold, FDR<0.05 | "more than 700" | 650 | slight undershoot |
| WT UP >10-fold, FDR<0.05 | "more than 300" | 561 | satisfies |
| **477-sig DOWN >2-fold, FDR<0.05** | **33.0 %** | **32.1 %** | **0.9 pp** |

The primary drift metric (DOWN direction on the 477 signature) matches
within 1 percentage point. Whole-transcriptome UP counts run 25-30 %
higher than Gosselin published, an unresolved methodological gap
attributable to some detail of his 2016-vintage count-generation
pipeline that we cannot fully reconstruct. All discussed in the
validation writeup.

---

## Requirements

- Python 3.11+
- R 4.0+ with limma and edgeR from Bioconductor
- (For Step 1 only) scanpy, anndata, pandas, numpy — see
  `requirements-python.txt`

## Citation

If you use this pipeline, please cite:

- Gosselin D et al. **An environment-dependent transcriptional network
  specifies human microglia identity.** *Science* 356, eaal3222 (2017).
  DOI: 10.1126/science.aal3222
- [Your paper, when published, goes here.]

## License

MIT.
