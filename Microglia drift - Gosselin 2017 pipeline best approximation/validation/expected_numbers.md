# Expected reproduction of Gosselin 2017's published numbers

Running `validation/run_validation.sh` on Gosselin 2017 Table S2 read
counts, restricted to the N=5 paired ex-vivo → 7-day in-vitro human
microglia design (see `dump_gosselin_counts.py` for the exact sample
list), should produce the following numbers.

If your local run lands within a small tolerance of these, the
pipeline is faithful on your machine.

## Primary drift headline (the "33 %" figure Gosselin published)

| metric                                                | Gosselin published | This pipeline | tolerance |
|-------------------------------------------------------|-------------------:|--------------:|----------:|
| fraction of 477-conserved signature DOWN > 2-fold FDR<.05 | **33.0 %**        | **32.1 %**   | ± 2 pp    |

Reproduces to within 1 percentage point.

## Whole-transcriptome DE counts (paired ex-vivo → 7 d, FDR<0.05)

| metric                           | Gosselin published | This pipeline | notes |
|----------------------------------|-------------------:|--------------:|-------|
| genes DOWN > 2-fold, FDR<0.05    | **3,702**          | **3,709**    | +7 (99.8 %) |
| genes UP   > 2-fold, FDR<0.05    | 2,263              | 2,913        | +29 % (see caveats) |
| genes DOWN > 10-fold, FDR<0.05   | "more than 700"    | 650          | slight undershoot |
| genes UP   > 10-fold, FDR<0.05   | "more than 300"    | 561          | satisfies bound |

## Caveats

The whole-transcriptome UP count runs 25-30 % higher than Gosselin
published. Extensive testing (see the repo history) confirmed this
is robust to every reasonable methodological variant of the stated
pipeline (voom vs voomWithQualityWeights, standard vs stricter
`filterByExpr`, TMM vs TMMwsp vs upperquartile vs no normalization,
with and without a > 200 nt gene-length filter). The residual gap is
attributable to some detail of Gosselin's 2016-vintage
count-generation pipeline that we cannot fully reconstruct from the
published Table S2 counts.

This does not compromise the drift comparison. The DOWN direction
(which is the direction the drift headline measures) reproduces
essentially exactly.

## Source of the "33 %"

Verbatim from Gosselin et al. 2017, *Science* eaal3222, p. 7:
> "Thirty-three percent of the conserved signature gene set and 31 %
> of the genes associated with AD risk variants exhibited > 2-fold
> reduction in expression in vitro (Figs. 4, B to D)."

And from the supplementary methods:
> "Within-species differential gene expression analysis was performed
> using the R package limma using the voom normalization (69).
> Significance was assessed at a false discovery rate (FDR) of 0.05
> using the Benjamini-Hochberg method and an effect-size cutoff of
> 2-fold change in expression unless otherwise noted."
