# Pipeline steps, in detail

A step-by-step walkthrough of the analysis, assuming no prior
bioinformatics background. The pipeline is identical for the
validation on Gosselin's data and for user data. Only the **inputs**
(which count matrix, which sample metadata) differ.

## The overall goal

For each donor, given microglia at two timepoints — a "fresh" one and
a "cultured" one — get two numbers out:

1. **Whole-transcriptome differential expression counts:** how many
   individual genes moved > 2-fold (or > 10-fold), in either direction,
   at FDR < 0.05.
2. **The drift headline:** what fraction of Gosselin's 477-gene
   conserved microglia identity signature (his Table S5) dropped
   > 2-fold in the cultured timepoint at FDR < 0.05. This is the same
   metric Gosselin's 33 % headline uses.

The two datasets this repo handles:

- **Gosselin 2017 (validation target).** Ex-vivo = day 0, freshly
  FACS-sorted microglia from resected human brain. In-vitro = day 7,
  same cells cultured 7 days in DMEM/F12 + IL-34. Bulk RNA-seq of the
  sorted microglia.
- **User data (e.g. OSCM).** Any paired D0-vs-later-timepoint
  microglia experiment. In our own use, D0 = fresh flash-frozen
  brain tissue and D14 = 14 days in organotypic slice culture,
  snRNA-seq of the whole slice with microglia identified by leiden
  clustering + marker score, per-donor-per-timepoint pseudobulk.

Both go through **the same limma-voom pipeline** independently. The
comparison at the end is the two resulting drift percentages next to
each other.

---

## Input the pipeline needs

**A gene × sample integer count matrix.**

- Rows: HGNC gene symbols.
- Columns: samples. One column per (donor, timepoint) — or per
  bulk RNA-seq library.
- Values: raw integer read counts (or for snRNA-seq, pseudobulk sums
  of raw UMI counts across all microglial nuclei from that donor ×
  timepoint).

**A sample metadata table** — one row per sample column, with:
- `donor`: an identifier for the patient (e.g. `HMG008`, `5-6-24`).
- `condition`: `ex` for the fresh/day-0 timepoint, `iv` for the
  cultured timepoint.

That's everything the DE pipeline needs.

For snRNA-seq input, `pipeline/01_h5_to_pseudobulk.py` handles the
h5 → pseudobulk conversion first (see that script's docstring).

---

## Step 1 — Filter to sensibly expressed genes

**What.** Drop genes with essentially no reads across all samples.

**Why.** Roughly a third of annotated human genes are undetected in
any given tissue or cell type. Testing them adds noise to the
multiple-testing correction (BH-FDR penalizes based on the total
number of tests) without contributing signal.

**How.** `edgeR::filterByExpr(y, design = ~donor + condition)`. Looks
at the design and keeps genes with at least ~10 counts in enough
samples for the smallest group in the design.

**Typical effect.** ~24 K genes → ~15 K genes.

---

## Step 2 — Compute library-size normalization factors (TMM)

**What.** Each sample gets one number — a "norm factor" — that
accounts for total sequencing depth and compositional imbalance.

**Why.** Two samples with identical biology can have very different
total read counts because of variable sequencing depth. Naive
per-total-reads scaling ("CPM") breaks when a few very high-expression
genes dominate one sample — every *other* gene in that sample then
looks artifactually lower. TMM (Trimmed Mean of M-values) accounts for
this by ignoring extreme fold changes when computing the scale
factor.

Relevant here because cultured microglia strongly induce a handful of
cytokines/stress genes (IL1B ~10-fold up, TNF, HSPA family). Without
TMM those induced genes would eat library-share, making every other
gene look artifactually lower at the in-vitro timepoint.

**How.** `edgeR::calcNormFactors(y, method = "TMM")`.

**Result.** Each sample gets a norm-factor number stored alongside
its library size. No count data has been modified yet.

---

## Step 3 — voom: raw counts → log2-CPM with per-observation precision weights

**What.** Transform integer counts into log2-scale expression values,
AND compute a "reliability weight" for each count in each sample.

**Why.** Count data has an intrinsic problem: low-count genes are
inherently noisier than high-count genes. A gene with 3 counts vs 6
counts across two samples is 2-fold in raw magnitude but very
unreliable; the same 2-fold change from 500 to 1000 counts is highly
reliable. Statistical tests that treat all measurements as equally
trustworthy over-count noisy low-expression genes as significant.

voom fits a lowess curve to the mean-variance relationship across all
genes and assigns each individual measurement a precision weight
proportional to how reliable that measurement is at its expression
level. Downstream regression uses those weights.

**How.** `limma::voom(y, design = ~donor + condition)`.

**Result.** For each (gene × sample) cell — a log2-CPM value plus a
weight. Both feed into Step 4.

---

## Step 4 — Fit a linear model per gene (this is where "paired" happens)

**What.** For each gene, fit the model

```
log2CPM  =  intercept  +  donor_effect  +  condition_effect
```

The `condition_effect` is the coefficient we care about — the average
log2 fold change between the cultured (`iv`) and fresh (`ex`) timepoint.

**Why the two-term design.** Including the `donor` term absorbs each
patient's baseline expression differences (patient A might just
express gene X higher than patient B at baseline). With `donor` in
the model, the `condition` term captures *within-donor* change from
ex-vivo to in-vitro — this is a **paired analysis** in the classical
sense. Without the donor term, patient-to-patient variance would
swamp the timepoint effect.

**How.** `limma::lmFit(v, design = ~donor + condition)`.

**Result.** For each gene, a fitted coefficient for the
`conditioniv` contrast (average log2 fold change) plus a raw per-gene
variance estimate.

---

## Step 5 — Empirical-Bayes variance moderation (the step that makes this work at small n)

**What.** Shrink each gene's per-gene variance estimate toward a
global prior fit across all genes.

**Why this is critical at small sample sizes.** At n = 3–5, per-gene
variance estimates are terrible. Random chance produces a gene with
near-zero variance across 5 numbers, making a small effect look
highly significant (false positive), or a gene with random huge
variance that hides a real effect (false negative).

Empirical Bayes borrows information across the whole transcriptome:
with ~15 K genes we can fit a global mean-variance trend and assume
each gene's true variance is drawn from that population. Each
individual gene's raw variance estimate is pulled toward the global
prior. This gives each gene *effective degrees of freedom in the tens
or hundreds* instead of the naive n−1.

This is why a naive paired t-test at n = 3 across 15 K tests returns
zero significant genes after BH-FDR (n−1 = 2 df; not enough power),
while limma-voom on the exact same n = 3 data returns thousands.

**How.** `limma::eBayes(fit)`.

**Result.** Moderated t-statistics, moderated p-values, and effective
degrees of freedom per gene.

---

## Step 6 — Extract results and correct for multiple testing

**What.** Per-gene results table: log2 fold change, raw p-value,
BH-adjusted p-value (FDR), average expression, moderated
t-statistic.

**Why FDR.** With ~15 K genes and raw p < 0.05 we'd expect ~750 false
positives by chance. BH-FDR controls the expected fraction of false
positives among calls (not the per-gene rate). Convention: call a
gene "significant" if BH-adjusted p < 0.05.

**How.** `limma::topTable(fit, coef = "conditioniv", number = Inf,
sort.by = "none", adjust.method = "BH")`.

---

## Step 7 — Apply Gosselin's exact thresholds

Call a gene **"significantly changed"** if BOTH:
- **FDR < 0.05** (multiple-testing corrected), AND
- **|log2 fold change| > 1** (i.e., > 2-fold in either direction).

Split into **DOWN** (logFC < −1 AND FDR < 0.05) and **UP**
(logFC > 1 AND FDR < 0.05). Also count > 10-fold (log2FC > log2(10) ≈
3.32) in both directions.

Verbatim from Gosselin's supplementary methods: *"Significance was
assessed at a false discovery rate (FDR) of 0.05 using the
Benjamini-Hochberg method and an effect-size cutoff of 2-fold change
in expression."* Applied identically to both datasets.

---

## Step 8 — Compute the drift headline number

**What.** Intersect the significantly-DOWN gene set with Gosselin's
477-gene conserved microglia signature (his Table S5). Divide by the
number of the 477 that were testable in the data (i.e., passed the
Step 1 filter).

**Result — the drift headline.** Fraction of the conserved microglia
identity signature that dropped > 2-fold at FDR < 0.05.

For Gosselin (his 5 paired donors, day 0 → day 7 dissociated
culture): **33.0 % published; 32.1 % reproduced by this pipeline.**

For user data: whatever number the same pipeline produces on the
user's D0 → later-timepoint samples. Put it next to 33 % — a lower
number means the user's culture preserves microglia identity better
than Gosselin's dissociated culture; a higher (or similar) number
means it doesn't.

---

## Interpreting the whole-transcriptome DE counts

Alongside the 477-signature headline, the pipeline reports:

- Number of genes changed > 2-fold at FDR < 0.05, DOWN and UP.
- Number of genes changed > 10-fold at FDR < 0.05, DOWN and UP.

Gosselin's published numbers for the ex-vivo → 7 d comparison:
3,702 DOWN and 2,263 UP at > 2-fold; "more than 700" DOWN and "more
than 300" UP at > 10-fold. Ours reproduce the DOWN counts to within
7 genes (99.8 %) and match the > 10-fold bounds. UP counts run
25-30 % higher; discussed in `validation/README.md`.

---

## Choices we deliberately did NOT vary

The following are held fixed at Gosselin's stated methodology and
should not be tweaked to match different results:

- Test: **limma-voom + eBayes** (not paired-t, not Wilcoxon, not
  DESeq2 — those all give different-scale answers).
- Design: **~donor + condition** (paired, blocking on donor).
- Filter: **edgeR::filterByExpr** with defaults.
- Normalization: **TMM**.
- Thresholds: **FDR < 0.05 AND |log2FC| > 1**.
- Signature: **Gosselin's 477-gene conserved signature from Table S5**.

Users applying this pipeline to their own data should not modify
these choices to force a particular comparison outcome. The pipeline's
value comes from matching Gosselin's methodology as closely as we can
verify.
