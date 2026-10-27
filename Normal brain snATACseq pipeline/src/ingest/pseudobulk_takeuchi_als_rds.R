#!/usr/bin/env Rscript
# Takeuchi/Nagano ALS motor cortex + spinal cord — pseudobulk from Seurat RDS.
# Source: figshare 10.6084/m9.figshare.29628653
#   brain_so_integrated_LinkePeaks.rds        (5.6 GB; multiome w/ ATAC assay 159,872 × 35,739)
#   spinalcord_so_integrated_LinkePeaks.rds   (3.2 GB; 136,971 × 26,061)
#
# Author-peak coords are the row names of the ATAC assay. Author cell barcodes don't match CxG obs
# directly:
#   Brain:  "MU01__AAACAGCCAAATTCGT-1" (double underscore) → CxG "MU01_AAACAGCCAAATTCGT-1"
#   Spinal: "AAACAGCCACAGGAAT-1_1"    (no prefix, _N suffix) → CxG "MU04_AAACAGCCACAGGAAT-1"
# We normalize to CxG format using md$Project (the Seurat object's own sample tag) and join to the
# CxG obs TSV to get the cell_type → major_class map.
#
# donor_id is tagged with tissue so same subject in brain + spinal cord contributes separately.
#
# Reconstructed from conversation transcript 2026-10-27.
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(GenomicRanges)
})

setwd("/private/tmp/takeuchi_als")  # adjust to local workspace with the two RDS files
UNION_BED <- "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v1.0.0/union_peaks.bed"
CXG_OBS   <- "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/audit_2026-10-25/acquisition/cxg_new/cxg-takeuchi-als_obs.tsv"
OUT_TSV   <- "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/h5ads_v2.2_nnz_rebuild/cxg-takeuchi-als_direct_nnz1k.tsv"
STUDY_ID  <- "cxg-takeuchi-als"
NNZ_FLOOR <- 1000
MIN_CELLS_PER_PEAK <- 3

# CxG cell_type → atlas major_class
CT_MAP <- c(
  "oligodendrocyte" = "Oligodendrocyte",
  "glutamatergic neuron" = "Excitatory",
  "astrocyte" = "Astrocyte",
  "GABAergic neuron" = "Inhibitory",
  "microglial cell" = "Microglia",
  "oligodendrocyte precursor cell" = "OPC",
  "endothelial cell" = "Endothelial",
  "central nervous system macrophage" = "Microglia",
  "fibroblast" = "Mural",
  "pericyte" = "Mural",
  "somatomotor neuron" = "Motor_neuron",
  "ependymal cell" = "Ependymal",
  "perivascular cell" = "Mural"
)
canonical <- c(paste0("chr", 1:22), "chrX", "chrY")

cat("[1/6] load union peaks\n")
peaks <- fread(UNION_BED, col.names = c("chrom", "start", "end", "name"))
peaks[, peak_num := as.integer(sub("peak_", "", name))]
peaks <- peaks[chrom %in% canonical]
setorder(peaks, chrom, start)
peaks[, ord := .I]
union_gr <- GRanges(peaks$chrom, IRanges(peaks$start + 1L, peaks$end))

cat("[2/6] load CxG obs + build cell_type map\n")
cxg_obs <- fread(CXG_OBS)
cxg_obs[, major_class := unname(CT_MAP[as.character(cell_type)])]

process_one <- function(rds_path, tissue_tag) {
  cat(sprintf("\n[3/6] load %s (tissue=%s)\n", rds_path, tissue_tag))
  suppressPackageStartupMessages({ library(Seurat); library(Signac) })
  s <- readRDS(rds_path)
  counts <- GetAssayData(s, assay = "ATAC", slot = "counts")
  md <- s@meta.data
  common <- intersect(colnames(counts), rownames(md))
  stopifnot(length(common) >= 0.9 * ncol(counts))
  counts <- counts[, common, drop = FALSE]
  md <- md[common, , drop = FALSE]

  if (tissue_tag == "brain") {
    cxg_bc <- gsub("__", "_", colnames(counts), fixed = TRUE)
  } else {
    trimmed <- sub("_[0-9]+$", "", colnames(counts))
    cxg_bc <- paste0(md$Project, "_", trimmed)
  }
  cxg_idx <- match(cxg_bc, cxg_obs$barcode)
  cell_donor <- ifelse(!is.na(cxg_idx), cxg_obs$donor_id[cxg_idx], md$SubjectID)
  cell_class <- ifelse(!is.na(cxg_idx), cxg_obs$major_class[cxg_idx], NA_character_)

  # Parse peak coords from rownames
  rn <- rownames(counts)
  m <- regmatches(rn, regexec("^(chr[0-9XYM]+)[-:](\\d+)[-](\\d+)$", rn))
  has <- vapply(m, length, 1L) == 4L
  author_df <- data.table(
    idx = seq_along(rn),
    chrom = vapply(m, function(x) if (length(x) == 4L) x[2] else NA_character_, ""),
    start = as.integer(vapply(m, function(x) if (length(x) == 4L) x[3] else NA_character_, "")),
    end = as.integer(vapply(m, function(x) if (length(x) == 4L) x[4] else NA_character_, ""))
  )
  author_df <- author_df[!is.na(start) & chrom %in% canonical]
  author_gr <- GRanges(author_df$chrom, IRanges(author_df$start + 1L, author_df$end))

  cat("[4/6] overlap author peaks → union peaks\n")
  hits <- findOverlaps(author_gr, union_gr, type = "any")
  M <- sparseMatrix(i = author_df$idx[queryHits(hits)],
                    j = subjectHits(hits), x = 1,
                    dims = c(nrow(counts), length(union_gr)))

  cat("[5/6] build cell × union_peak binary + nnz filter\n")
  counts_bin <- as(counts > 0, "CsparseMatrix")
  cell_union <- t(counts_bin) %*% M
  cell_union_bin <- as(cell_union > 0, "CsparseMatrix")
  rownames(cell_union_bin) <- colnames(counts)  # explicit, matmul doesn't preserve
  rm(counts, counts_bin, cell_union, M, s); gc(verbose = FALSE)
  cell_nnz <- rowSums(cell_union_bin)
  keep <- cell_nnz >= NNZ_FLOOR & !is.na(cell_class)

  list(bin = cell_union_bin, donor = cell_donor, class = cell_class, keep = keep)
}

accumulate <- function(bins, donors, classes, keeps) {
  cat("\n[6/6] aggregate (donor, class) across tissues\n")
  rows_list <- list()
  for (tissue in names(bins)) {
    bin <- bins[[tissue]]
    donors_kept <- paste0(donors[[tissue]][keeps[[tissue]]], "_", tissue)
    classes_kept <- classes[[tissue]][keeps[[tissue]]]
    bin_kept <- bin[keeps[[tissue]], , drop = FALSE]
    grp_key <- paste(donors_kept, classes_kept, sep = "|")
    for (g in unique(grp_key)) {
      idx <- which(grp_key == g)
      sub <- bin_kept[idx, , drop = FALSE]
      n_access <- colSums(sub)
      active <- which(n_access >= MIN_CELLS_PER_PEAK)
      if (length(active) == 0) next
      parts <- strsplit(g, "|", fixed = TRUE)[[1]]
      rows_list[[g]] <- data.table(
        peak_num = peaks$peak_num[active],
        major_class = parts[2],
        donor_id = parts[1],
        n_cells_access_direct = as.integer(n_access[active]),
        n_cells_total_direct = length(idx),
        fraction_open_direct = n_access[active] / length(idx),
        study_id = STUDY_ID
      )
    }
  }
  out <- rbindlist(rows_list)
  fwrite(out, OUT_TSV, sep = "\t")
  cat(sprintf("wrote %s (%d rows, %d donors)\n", OUT_TSV, nrow(out),
              length(unique(out$donor_id))))
}

bins <- list(); donors <- list(); classes <- list(); keeps <- list()
for (tt in list(c("brain", "brain_so_integrated_LinkePeaks.rds"),
                c("spinalcord", "spinalcord_so_integrated_LinkePeaks.rds"))) {
  r <- process_one(tt[2], tt[1])
  bins[[tt[1]]] <- r$bin; donors[[tt[1]]] <- r$donor
  classes[[tt[1]]] <- r$class; keeps[[tt[1]]] <- r$keep
  gc(verbose = FALSE)
}
accumulate(bins, donors, classes, keeps)
