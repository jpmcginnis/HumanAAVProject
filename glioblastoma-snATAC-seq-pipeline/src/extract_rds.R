#!/usr/bin/env Rscript
# Extract an R-serialized single-cell ATAC object → h5ad for the Python pipeline.
#
# Supported object classes:
#   - Seurat with ChromatinAssay (Signac) — e.g. HRA004942's gbm_atac_signac_ArchR.rds
#   - SnapATAC v1 "snap" S4 object — e.g. GSE165037's GSMxxxx_*_IDHWT_snATAC.rds.gz
#   - ArchRProject — ArchR's own project class
#
# Output: an h5ad whose X is the peak × cell counts matrix (CSC), with .var
# carrying chrom/start/end and .obs carrying all cell-level metadata including
# cell_type labels when present. Downstream ingest handles the normalization
# step that lifts per-sample peaks onto CATLAS space (via src/ingest_matrix.py).
#
# Usage:
#   Rscript src/extract_rds.R --rds <input.rds[.gz]> --out <output.h5ad> \
#       [--sample-id GBM4] [--idh-status wildtype] [--primary-recurrent primary]
#
# Writes a plain anndata-format h5ad using the anndata R package (installed
# alongside Signac). Falls back to writing mtx+barcodes+peaks files if anndata
# isn't available, which src/ingest_matrix.py already knows how to read.

suppressWarnings({
  .libPaths(c("/data/R_libs", .libPaths()))
})

# Define a stub `snap` S4 class so readRDS of GSE165037's SnapATAC v1 objects
# works without the SnapATAC package installed (its install depends on raster
# → terra → system GDAL/PROJ which are not easy to install cleanly). We only
# need slot access, not any SnapATAC methods, so defining the class with the
# slots we touch is sufficient.
suppressWarnings(
  if (!isClass("snap")) {
    methods::setClass("snap", representation = list(
      pmat     = "ANY",
      bmat     = "ANY",
      peak     = "ANY",
      feature  = "ANY",
      barcode  = "ANY",
      metaData = "ANY",
      cluster  = "ANY",
      des      = "ANY",
      file     = "ANY",
      sample   = "ANY"
    ))
  }
)

parse_args <- function(argv) {
  out <- list(rds=NULL, out=NULL, sample_id=NA, idh_status="wildtype",
              primary_recurrent="primary", cohort=NA)
  i <- 1
  while (i <= length(argv)) {
    k <- argv[i]
    v <- argv[i+1]
    if (k == "--rds") out$rds <- v
    else if (k == "--out") out$out <- v
    else if (k == "--sample-id") out$sample_id <- v
    else if (k == "--idh-status") out$idh_status <- v
    else if (k == "--primary-recurrent") out$primary_recurrent <- v
    else if (k == "--cohort") out$cohort <- v
    i <- i + 2
  }
  if (is.null(out$rds) || is.null(out$out))
    stop("--rds and --out are required")
  out
}

log_msg <- function(msg) cat(sprintf("[extract-rds] %s\n", msg))

load_rds_any <- function(path) {
  log_msg(sprintf("Loading %s (%.1f MB)", path, file.info(path)$size / 1e6))
  # GSE165037's `.rds.gz` files are DOUBLE-gzipped: a gzip wrapper around an
  # already gzip-serialized RDS. R's readRDS handles single-nested gzip, so we
  # un-wrap the outer gzip here first (if a nested-gzip is detected) and pass
  # a plain file handle to readRDS.
  if (grepl("\\.gz$", path)) {
    # Peek at the inner stream: if the first two bytes after the outer gunzip
    # are the gzip magic 1f 8b, we have double-gzip.
    con_peek <- gzfile(path, "rb")
    hdr <- readBin(con_peek, "raw", n = 2L)
    close(con_peek)
    double_gzip <- length(hdr) == 2L && hdr[1] == as.raw(0x1f) && hdr[2] == as.raw(0x8b)

    if (double_gzip) {
      log_msg("  detected double-gzip; unwrapping outer gzip to a plain .rds first")
      plain_path <- sub("\\.gz$", "", path)
      if (!file.exists(plain_path)) {
        system2("gunzip", c("-c", path), stdout = plain_path)
      }
      obj <- readRDS(plain_path)
    } else {
      con <- gzfile(path, "rb")
      obj <- readRDS(con)
      close(con)
    }
  } else {
    obj <- readRDS(path)
  }
  log_msg(sprintf("Loaded object of class: %s", paste(class(obj), collapse=", ")))
  obj
}

# ---------- extractors per object class ----------

extract_from_seurat <- function(obj) {
  # Signac Seurat: pull the ChromatinAssay's counts matrix + peak coords
  if (!requireNamespace("Signac", quietly=TRUE))
    stop("Signac R package not installed")
  log_msg("Extracting Seurat ChromatinAssay")
  # Prefer the 'peaks' assay; fall back to 'ATAC' or whatever is default
  assay_names <- Seurat::Assays(obj)
  log_msg(sprintf("  assays: %s", paste(assay_names, collapse=", ")))
  atac_assay <- NULL
  for (a in c("peaks", "ATAC", "atac", "chromatin")) {
    if (a %in% assay_names) { atac_assay <- a; break }
  }
  if (is.null(atac_assay)) atac_assay <- assay_names[1]
  log_msg(sprintf("  using assay: %s", atac_assay))

  # Seurat 5 deprecated `slot=` in GetAssayData; try new API first, fall back
  counts <- tryCatch(
    SeuratObject::LayerData(obj[[atac_assay]], layer="counts"),
    error = function(e) tryCatch(
      Seurat::GetAssayData(obj, assay=atac_assay, layer="counts"),
      error = function(e2) Seurat::GetAssayData(obj, assay=atac_assay, slot="counts")
    )
  )
  log_msg(sprintf("  counts: %d peaks × %d cells  nnz=%d", nrow(counts), ncol(counts), length(counts@x)))

  # Peak coordinates — Signac stores these in the ChromatinAssay's ranges
  granges <- Signac::granges(obj[[atac_assay]])
  var_df <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(granges)),
    start = as.integer(GenomicRanges::start(granges)),
    end   = as.integer(GenomicRanges::end(granges)),
    stringsAsFactors = FALSE
  )
  obs_df <- as.data.frame(obj@meta.data)
  list(counts=counts, var=var_df, obs=obs_df)
}

extract_from_snap <- function(obj) {
  # SnapATAC v1 S4 object. Use methods::slot() which works even when the
  # `snap` class is defined only via our stub above — we never need any
  # SnapATAC method, just raw slot access.
  log_msg("Extracting SnapATAC v1 object via raw slot access")

  pmat <- tryCatch(methods::slot(obj, "pmat"), error = function(e) NULL)
  bmat <- tryCatch(methods::slot(obj, "bmat"), error = function(e) NULL)
  peak_gr <- tryCatch(methods::slot(obj, "peak"), error = function(e) NULL)
  feature_gr <- tryCatch(methods::slot(obj, "feature"), error = function(e) NULL)
  barcode <- tryCatch(methods::slot(obj, "barcode"), error = function(e) NULL)
  metaData <- tryCatch(methods::slot(obj, "metaData"), error = function(e) NULL)
  cluster <- tryCatch(methods::slot(obj, "cluster"), error = function(e) NULL)

  # Choose matrix: prefer @pmat (peaks), fall back to @bmat (5kb bins) if peaks
  # are empty — some SnapATAC objects publish only the bin matrix.
  mat <- NULL; gr <- NULL; src <- NULL
  if (!is.null(pmat) && is(pmat, "Matrix") && length(pmat@x) > 0) {
    mat <- pmat; gr <- peak_gr; src <- "pmat"
  } else if (!is.null(bmat) && is(bmat, "Matrix") && length(bmat@x) > 0) {
    mat <- bmat; gr <- feature_gr; src <- "bmat"
    log_msg("  @pmat empty, falling back to @bmat (5kb bins)")
  } else {
    stop("Neither @pmat nor @bmat has data in this snap object")
  }

  # SnapATAC convention: matrix is cells × features. Transpose to features × cells.
  mat <- Matrix::t(mat)
  log_msg(sprintf("  using @%s transposed: %d features × %d cells  nnz=%d",
                  src, nrow(mat), ncol(mat), length(mat@x)))

  if (is.null(gr) || length(gr) == 0)
    stop(sprintf("@%s present but matching coord slot (@peak/@feature) is empty", src))

  var_df <- data.frame(
    chrom = as.character(GenomicRanges::seqnames(gr)),
    start = as.integer(GenomicRanges::start(gr)),
    end   = as.integer(GenomicRanges::end(gr)),
    stringsAsFactors = FALSE
  )

  obs_df <- data.frame(
    barcode = as.character(if (is.null(barcode)) seq_len(ncol(mat)) else barcode),
    stringsAsFactors = FALSE
  )
  if (!is.null(metaData) && is.data.frame(metaData) && nrow(metaData) == nrow(obs_df)) {
    obs_df <- cbind(obs_df, as.data.frame(metaData))
  }
  if (!is.null(cluster) && length(cluster) == nrow(obs_df)) {
    obs_df$cluster <- as.character(cluster)
  }
  # SnapATAC v1 does not enforce unique barcodes within an object (multiplexed
  # libraries may collide). Make rownames unique via make.unique() so the
  # downstream obs_df can be rowname-indexed.
  rownames(obs_df) <- make.unique(obs_df$barcode, sep = "-dup")
  list(counts=mat, var=var_df, obs=obs_df)
}

extract_from_archr <- function(obj) {
  stop("ArchRProject extraction not yet implemented — use ArchR::getMatrixFromProject")
}

# ---------- writer ----------

write_h5ad <- function(data, out_path, extra_obs=list()) {
  # Augment obs with cohort-level provenance fields the Python pipeline expects
  obs <- data$obs
  for (k in names(extra_obs)) {
    val <- extra_obs[[k]]
    if (!is.na(val)) obs[[k]] <- val
  }
  # Ensure required columns exist
  if (!"sample_id" %in% colnames(obs) && !is.na(extra_obs$sample_id)) {
    obs$sample_id <- extra_obs$sample_id
  }
  if (!"region" %in% colnames(obs)) obs$region <- "unknown"
  if (!"modality" %in% colnames(obs)) obs$modality <- "snATAC_rds_extracted"

  # Var: name as chr:start-end so downstream ingest_matrix can parse
  var <- data$var
  var$peak_id <- paste0(var$chrom, ":", var$start, "-", var$end)
  rownames(var) <- var$peak_id

  # Use anndata R package if present; otherwise write mtx + bed + tsv
  if (requireNamespace("anndata", quietly=TRUE)) {
    log_msg("Writing h5ad via anndata R")
    a <- anndata::AnnData(
      X = Matrix::t(data$counts),       # anndata expects cells × features
      obs = obs,
      var = var
    )
    anndata::write_h5ad(a, out_path, compression="gzip")
    log_msg(sprintf("Wrote %s (%.1f MB)", out_path, file.info(out_path)$size / 1e6))
  } else {
    # Fallback: write as mtx + bed + tsv so src/ingest_matrix.py can read it
    log_msg("anndata R not available; writing mtx + bed + tsv (src/ingest_matrix.py reads these)")
    out_dir <- gsub("\\.h5ad$", "", out_path)
    dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
    Matrix::writeMM(data$counts, file.path(out_dir, "matrix.mtx"))
    system2("gzip", c("-f", file.path(out_dir, "matrix.mtx")))
    write.table(var[, c("chrom","start","end")], file=gzfile(file.path(out_dir, "peaks.bed.gz")),
                sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(data.frame(barcode=rownames(obs)), file=gzfile(file.path(out_dir, "barcodes.tsv.gz")),
                sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
    write.table(obs, file=gzfile(file.path(out_dir, "obs.tsv.gz")),
                sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
    log_msg(sprintf("Wrote mtx/bed/tsv under %s", out_dir))
  }
}

# ---------- main ----------

args <- parse_args(commandArgs(trailingOnly=TRUE))
log_msg(sprintf("=== extract_rds: %s ===", basename(args$rds)))
obj <- load_rds_any(args$rds)
klass <- class(obj)[1]

data <- NULL
# Check snap FIRST — `inherits(obj, "Seurat")` on a snap object triggers S4
# class resolution, which attempts to load the SnapATAC package (which we
# deliberately avoid installing; see the stub class at the top of this file).
if (klass == "snap") {
  data <- extract_from_snap(obj)
} else if (inherits(obj, "Seurat")) {
  data <- extract_from_seurat(obj)
} else if (inherits(obj, "ArchRProject")) {
  data <- extract_from_archr(obj)
} else {
  stop(sprintf("Unsupported object class: %s", paste(class(obj), collapse=", ")))
}

extra <- list(
  sample_id = args$sample_id,
  idh_status = args$idh_status,
  primary_recurrent = args$primary_recurrent,
  cohort = args$cohort
)
write_h5ad(data, args$out, extra_obs=extra)
log_msg("=== DONE ===")
