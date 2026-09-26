#!/usr/bin/env Rscript
# Step 2 of the pipeline: pseudobulk counts + sample metadata → per-gene
# DE stats using limma-voom, exactly as described in Gosselin 2017's
# supplementary methods.
#
# Method:
#   * edgeR::filterByExpr with the paired design
#   * calcNormFactors(method = "TMM")
#   * voom(y, design)
#   * lmFit(v, design = ~donor + condition)
#   * eBayes(fit)
#   * topTable with BH-FDR
#
# Metadata `condition` values must be exactly "ex" and "iv" (ex-vivo /
# fresh vs in-vitro / cultured); the tested coefficient is "conditioniv".
#
# Outputs:
#   limma_voom_de.csv        per-gene log2FC, p, adj.P.Val, AveExpr, t
#   limma_voom_summary.csv   whole-transcriptome DE counts summary
#
# Usage:
#   Rscript 02_limma_voom.R --counts pseudobulk_counts.csv \
#                            --metadata sample_metadata.csv \
#                            --outdir output/

suppressPackageStartupMessages({
  library(optparse)
  library(limma)
  library(edgeR)
})

opts <- OptionParser() |>
  add_option(c("--counts"),   type="character", help="pseudobulk_counts.csv") |>
  add_option(c("--metadata"), type="character", help="sample_metadata.csv") |>
  add_option(c("--outdir"),   type="character", help="output directory") |>
  add_option(c("--min-count"), type="integer", default=NULL,
             help="override filterByExpr min.count (default: edgeR default)") |>
  parse_args()

stopifnot(!is.null(opts$counts), !is.null(opts$metadata), !is.null(opts$outdir))
dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

# Load
counts <- read.csv(opts$counts, row.names=1, check.names=FALSE)
meta   <- read.csv(opts$metadata, stringsAsFactors=FALSE)
stopifnot(all(c("sample","donor","condition") %in% colnames(meta)))
stopifnot(all(meta$condition %in% c("ex","iv")))
stopifnot(all(meta$sample %in% colnames(counts)))
counts <- as.matrix(counts[, meta$sample, drop=FALSE])
storage.mode(counts) <- "integer"
cat(sprintf("input: %d genes x %d samples (%d donors)\n",
             nrow(counts), ncol(counts), length(unique(meta$donor))))

# Enforce factor levels: ex is reference, iv is the test coefficient
meta$condition <- factor(meta$condition, levels=c("ex","iv"))
meta$donor     <- factor(meta$donor)

# Step 1: build DGEList and filter
y <- DGEList(counts=counts, samples=meta, group=meta$condition)
design <- model.matrix(~ donor + condition, data=meta)
if (is.null(opts$`min-count`)) {
  keep <- filterByExpr(y, design=design)
} else {
  keep <- filterByExpr(y, design=design, min.count=opts$`min-count`)
}
cat(sprintf("filterByExpr kept %d / %d genes\n", sum(keep), length(keep)))
y <- y[keep, , keep.lib.sizes=FALSE]

# Step 2: TMM norm factors
y <- calcNormFactors(y, method="TMM")

# Step 3: voom
v <- voom(y, design, plot=FALSE)

# Steps 4-5: fit + eBayes
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Step 6: extract, BH-FDR
coef_name <- grep("conditioniv", colnames(design), value=TRUE)
stopifnot(length(coef_name) == 1)
tt <- topTable(fit, coef=coef_name, number=Inf, sort.by="none",
                adjust.method="BH")
tt$gene <- rownames(tt)
tt <- tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")]
write.csv(tt, file.path(opts$outdir, "limma_voom_de.csv"), row.names=FALSE)

# Step 7: whole-transcriptome DE counts summary
sig <- tt$adj.P.Val < 0.05
summary_row <- data.frame(
  n_donors=length(unique(meta$donor)),
  n_samples=nrow(meta),
  n_genes_tested=nrow(tt),
  n_fdr_lt_05=sum(sig),
  wt_down_2fold  = sum(sig & tt$logFC < -1),
  wt_up_2fold    = sum(sig & tt$logFC >  1),
  wt_down_10fold = sum(sig & tt$logFC < -log2(10)),
  wt_up_10fold   = sum(sig & tt$logFC >  log2(10))
)
write.csv(summary_row, file.path(opts$outdir, "limma_voom_summary.csv"),
           row.names=FALSE)

cat("\n=== Whole-transcriptome DE summary (FDR<0.05) ===\n")
cat(sprintf("  genes tested:       %d\n", nrow(tt)))
cat(sprintf("  genes at FDR<.05:   %d\n", sum(sig)))
cat(sprintf("  >2-fold  DOWN=%d   UP=%d\n",
             summary_row$wt_down_2fold, summary_row$wt_up_2fold))
cat(sprintf("  >10-fold DOWN=%d   UP=%d\n",
             summary_row$wt_down_10fold, summary_row$wt_up_10fold))
cat(sprintf("\nwrote: %s/limma_voom_de.csv\n", opts$outdir))
cat(sprintf("       %s/limma_voom_summary.csv\n", opts$outdir))
