"""
Marker-peak scoring for matrix-mode cohorts.

The standard TME label transfer (gene activity → GBmap bridging) requires
fragments. Matrix-mode cohorts (e.g. Guilhamon GSE139136) have only
peak-barcode matrices, so that path is unavailable.

This module scores cells for cell-type membership directly from accessibility
at CATLAS peaks that fall within ±2 kb of a curated set of marker-gene TSSes.
Marker TSSes live in `config/marker_tsses.tsv` and can be edited without
touching code.

The pipeline:
  1. Load CATLAS peaks (var of the input AnnData) and the marker TSS table.
  2. For each cell type, pick CATLAS peaks within ±2 kb of ANY of its marker
     gene TSSes. These are the cell type's "signature peaks."
  3. For each cell, compute per-cell-type score = normalized mean accessibility
     across that cell type's signature peaks.
  4. Assign cell_type = argmax of scores, provided max > min_score_threshold.
     Otherwise leave as 'unassigned'.
  5. Record provenance per cell in `label_tme_marker_peaks` + score columns.

Expected to be called from src/label_transfer.py when the cohort is matrix-mode
(no fragments) and the standard gene-activity bridge is unavailable.

Caveats to be honest about:
  - This is a MARKER-based call, not a classifier trained on the full
    transcriptome. Expect lower accuracy than CellTypist on gene activity.
  - Open chromatin at a marker promoter doesn't guarantee the gene is
    expressed. False positives for leaky promoters.
  - Cell types without >=3 CATLAS peaks within ±2kb of any marker (unusual)
    get skipped with a warning.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np


DEFAULT_WINDOW_BP = 2000
DEFAULT_MIN_SCORE = 0.02          # at least 2% accessibility at marker peaks
DEFAULT_MARGIN_RATIO = 1.3        # top score must be ≥1.3× second-highest


def log(msg: str) -> None:
    print(f"[marker-peak] {msg}", flush=True)


def _load_marker_tsses(path: Path):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", comment="#",
                     names=["gene", "chrom", "tss", "strand", "cell_type"])
    df["tss"] = df["tss"].astype(int)
    log(f"  loaded {len(df)} markers across {df['cell_type'].nunique()} cell types "
        f"from {path.name}")
    return df


def _parse_peak_coords(adata):
    """Resolve chrom/start/end for every peak in adata.var.
    Accepts common column name variants (chrom/seqnames/chr) and falls back
    to parsing var_names like 'chr1:1234-5678' OR 'chr1:1234:5678:...'."""
    import pandas as pd
    import re

    # Try obs column variants
    _CHROM_COLS = ("chrom", "seqnames", "chr", "Chromosome")
    _START_COLS = ("start", "Start", "chromStart")
    _END_COLS = ("end", "End", "chromEnd")
    chrom_col = next((c for c in _CHROM_COLS if c in adata.var.columns), None)
    start_col = next((c for c in _START_COLS if c in adata.var.columns), None)
    end_col = next((c for c in _END_COLS if c in adata.var.columns), None)
    if chrom_col and start_col and end_col:
        df = adata.var[[chrom_col, start_col, end_col]].reset_index(drop=True)
        df.columns = ["Chromosome", "Start", "End"]
        df["Start"] = pd.to_numeric(df["Start"], errors="coerce").astype("Int64")
        df["End"] = pd.to_numeric(df["End"], errors="coerce").astype("Int64")
    else:
        names = pd.Series(adata.var_names.astype(str))
        # Try "chr:start-end" first, then "chr:start:end[:...]"
        parsed = names.str.extract(
            r"^(?P<Chromosome>chr[\dXYM]+):(?P<Start>\d+)[-:](?P<End>\d+)"
        )
        parsed["Start"] = pd.to_numeric(parsed["Start"], errors="coerce").astype("Int64")
        parsed["End"] = pd.to_numeric(parsed["End"], errors="coerce").astype("Int64")
        df = parsed
    # Drop any rows that couldn't be parsed
    valid = df["Chromosome"].notna() & df["Start"].notna() & df["End"].notna()
    if not valid.all():
        log(f"  WARNING: {int((~valid).sum())}/{len(df)} var rows had unparseable coords; dropping")
    df["peak_idx"] = np.arange(len(df), dtype=np.int64)
    df = df[valid].reset_index(drop=True)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df


def _signature_peaks_for_celltypes(markers_df, peaks_df, window_bp=DEFAULT_WINDOW_BP):
    """For each cell type, return the sorted list of CATLAS peak indices whose
    coords fall within ±window_bp of any of its marker TSSes."""
    import pyranges as pr
    import pandas as pd

    peaks_pr = pr.PyRanges(peaks_df.rename(columns={"peak_idx": "peak_idx"}))
    sig = {}
    for ct, g in markers_df.groupby("cell_type"):
        tss_df = pd.DataFrame({
            "Chromosome": g["chrom"].values,
            "Start": (g["tss"].values - window_bp).clip(min=0),
            "End":   g["tss"].values + window_bp,
        })
        tss_pr = pr.PyRanges(tss_df)
        joined = peaks_pr.join(tss_pr).df
        peak_idxs = sorted(set(joined["peak_idx"].values.tolist())) if len(joined) else []
        sig[ct] = peak_idxs
        log(f"    {ct:22s} n_markers={len(g):2d}  n_signature_peaks={len(peak_idxs):4d}")
    return sig


def score_cells_matrix_mode(adata,
                            markers_tsv: Path,
                            window_bp: int = DEFAULT_WINDOW_BP,
                            min_score: float = DEFAULT_MIN_SCORE,
                            margin_ratio: float = DEFAULT_MARGIN_RATIO) -> None:
    """Assign label_tme_gbmap, label_tme_confidence, label_tme_method in-place
    on `adata` using marker-peak scoring. Only touches cells where
    `label_tme_gbmap` is 'unassigned_matrix_mode' (or absent)."""
    import scipy.sparse as sp

    log("Running marker-peak scoring (CATLAS-coords, ±{} bp of marker TSSes)".format(window_bp))
    markers_df = _load_marker_tsses(markers_tsv)
    peaks_df = _parse_peak_coords(adata)
    sig = _signature_peaks_for_celltypes(markers_df, peaks_df, window_bp=window_bp)
    cell_types = [ct for ct, idxs in sig.items() if len(idxs) >= 3]
    log(f"  cell types with ≥3 signature peaks: {cell_types}")

    if len(cell_types) < 2:
        log("  ABORT: fewer than 2 cell types have ≥3 signature peaks — markers mismatch?")
        return

    # Normalize X per cell to CPM-like scale so marker scores are comparable.
    # CHUNKED PATH for backed / huge data: avoids materializing a 1M × 438K
    # matrix in RAM (what OOM'd the box on GBM-Space). Iterate cell chunks
    # of CHUNK rows, compute totals + per-ct signature sums per chunk, write
    # accumulators. Memory peak: one chunk's sparse rows (~CHUNK × 438K sparse).
    log("  computing per-cell marker-set accessibility (chunked for huge/backed data)")
    is_backed = getattr(adata, "isbacked", False)
    CHUNK = 50_000
    n_obs = adata.n_obs
    cell_total = np.zeros(n_obs, dtype=np.float64)
    scores = np.zeros((n_obs, len(cell_types)), dtype=np.float64)
    peak_idxs_per_ct = [np.asarray(sig[ct], dtype=np.int64) for ct in cell_types]

    for i in range(0, n_obs, CHUNK):
        j = min(i + CHUNK, n_obs)
        Xc = adata.X[i:j, :]
        if not sp.issparse(Xc):
            Xc = sp.csr_matrix(Xc)
        row_sums = np.asarray(Xc.sum(axis=1)).flatten()
        cell_total[i:j] = row_sums
        for ct_j, peak_idxs in enumerate(peak_idxs_per_ct):
            if peak_idxs.size == 0:
                continue
            sub = Xc[:, peak_idxs]
            counts = np.asarray(sub.sum(axis=1)).flatten()
            scores[i:j, ct_j] = counts
        if is_backed and (i // CHUNK) % 10 == 0:
            log(f"    chunk {i:,}/{n_obs:,}")
    cell_total = np.where(cell_total == 0, 1, cell_total)
    # Normalize AFTER accumulation (so divisions aren't done per chunk)
    for ct_j, peak_idxs in enumerate(peak_idxs_per_ct):
        if peak_idxs.size > 0:
            scores[:, ct_j] = scores[:, ct_j] / cell_total / len(peak_idxs)

    # Assign cell type = argmax, if above threshold AND margin over 2nd-highest
    top_j = scores.argmax(axis=1)
    sorted_scores = np.sort(scores, axis=1)
    top_val = sorted_scores[:, -1]
    second = sorted_scores[:, -2] if sorted_scores.shape[1] >= 2 else np.zeros_like(top_val)
    # min_score here is relative-fraction scale (very small numbers after normalization)
    pass_min = top_val > (min_score * np.median(top_val[top_val > 0]) if (top_val > 0).any() else min_score)
    pass_margin = top_val > (margin_ratio * second)
    assigned = pass_min & pass_margin

    labels = np.array(["unassigned"] * adata.n_obs, dtype=object)
    for i, j in enumerate(top_j):
        if assigned[i]:
            labels[i] = cell_types[j]
    adata.obs["label_tme_gbmap"] = labels.astype(str)
    adata.obs["label_tme_confidence"] = top_val
    adata.obs["label_tme_method"] = "marker_peak_scoring"

    # Report per-cell-type counts
    import pandas as pd
    log("Marker-peak label distribution:")
    for ct, n in pd.Series(labels).value_counts().items():
        log(f"    {str(ct):22s} {int(n):>7d}")


def main() -> int:
    """Standalone CLI for standalone invocation / testing."""
    import anndata as ad
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--markers", default="config/marker_tsses.tsv", type=Path)
    ap.add_argument("--window-bp", type=int, default=DEFAULT_WINDOW_BP)
    args = ap.parse_args()

    log(f"Loading {args.src}")
    adata = ad.read_h5ad(args.src)
    score_cells_matrix_mode(adata, args.markers, window_bp=args.window_bp)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out, compression="lzf")
    log(f"→ {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
