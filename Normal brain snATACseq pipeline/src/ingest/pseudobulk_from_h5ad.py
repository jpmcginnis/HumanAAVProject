#!/usr/bin/env python3
"""Build nnz-filtered per-donor-per-class per-peak pseudobulk from an h5ad.

For each cohort's annotated per-cell h5ad:
  1. Load cells × peaks sparse matrix (CSR, int8/int)
  2. Compute per-cell nnz = (X > 0).sum(axis=1)
  3. Keep cells with nnz >= 1000 (universal snATAC QC floor)
  4. For each (donor_id, major_class):
       n_cells_total = cells in that donor-class after nnz filter
       for each peak:
          n_cells_access = cells with X[cell, peak] > 0
          fraction_open  = n_cells_access / n_cells_total
  5. Write a parquet in the same column schema as the existing direct_parquets:
       peak_num, major_class, donor_id, n_cells_access_direct, n_cells_total_direct,
       fraction_open_direct, study_id
     (So nb_rescore_gbm_style.py works without changes.)

Processes cells in chunks to stay memory-friendly on a laptop.

Usage:
    python3 nb_pseudobulk_from_h5ad.py <input.h5ad> <output.parquet> [--nnz 1000] [--chunk 50000]
"""
import argparse
import gc
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import polars as pl
import scipy.sparse as sp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Input h5ad path")
    ap.add_argument("output", help="Output parquet path")
    ap.add_argument("--nnz", type=int, default=1000, help="Per-cell nnz floor (default 1000)")
    ap.add_argument("--chunk", type=int, default=50000, help="Cells per chunk (default 50000)")
    args = ap.parse_args()

    t0 = time.time()
    print(f"[1/4] Loading h5ad (backed) {args.input}", flush=True)
    a = ad.read_h5ad(args.input, backed="r")
    print(f"  shape: {a.shape} (cells × peaks)", flush=True)
    print(f"  study_id values: {a.obs['study_id'].unique().tolist() if 'study_id' in a.obs.columns else '(missing)'}", flush=True)
    print(f"  n unique donors: {a.obs['donor_id'].nunique()}", flush=True)
    print(f"  n unique major_class: {a.obs['major_class'].nunique()}", flush=True)

    n_cells, n_peaks = a.shape

    # peak_num convention: scored parquet uses "peak_N" strings; direct_parquets use integer peak_num
    # var rows are peak_N in order; peak_num = integer part
    if "peak_id" in a.var.columns:
        peak_ids = a.var["peak_id"].to_numpy()
    else:
        peak_ids = np.array(a.var_names)
    # extract integer
    peak_nums = np.array([int(p.replace("peak_", "")) for p in peak_ids], dtype=np.int64)

    # Pull donor_id, major_class, study_id as arrays for all cells
    donor_arr = a.obs["donor_id"].to_numpy()
    class_arr = a.obs["major_class"].to_numpy()
    study_arr = a.obs["study_id"].to_numpy() if "study_id" in a.obs.columns else np.array(["unknown"] * n_cells)
    # Single study per cohort usually; take mode
    study_id = pl.Series(study_arr).mode()[0]
    print(f"  study_id: {study_id}", flush=True)

    print(f"\n[2/4] Streaming nnz filter (keep nnz >= {args.nnz}), {args.chunk} cells/chunk", flush=True)
    pass_mask = np.zeros(n_cells, dtype=bool)
    cell_nnz = np.zeros(n_cells, dtype=np.int32)
    n_chunks = (n_cells + args.chunk - 1) // args.chunk
    for ci, start in enumerate(range(0, n_cells, args.chunk)):
        end = min(start + args.chunk, n_cells)
        X_chunk = a.X[start:end]
        if not sp.issparse(X_chunk):
            X_chunk = sp.csr_matrix(X_chunk)
        chunk_nnz = X_chunk.getnnz(axis=1)
        cell_nnz[start:end] = chunk_nnz
        pass_mask[start:end] = chunk_nnz >= args.nnz
        print(f"  chunk {ci+1}/{n_chunks}: {end-start:,} cells, nnz>={args.nnz}: {pass_mask[start:end].sum():,} ({100*pass_mask[start:end].sum()/(end-start):.1f}%)", flush=True)

    print(f"  TOTAL kept: {pass_mask.sum():,} / {n_cells:,} ({100*pass_mask.sum()/n_cells:.1f}%)", flush=True)

    print(f"\n[3/4] Streaming per-(donor,class) × peak aggregation", flush=True)
    # Group cells by (donor, class). Build index arrays.
    # Donor-class key = donor_id | major_class (string)
    dc_keys = np.array([f"{d}|{c}" for d, c in zip(donor_arr, class_arr)])
    unique_dc = np.unique(dc_keys[pass_mask])
    print(f"  unique (donor, class) groups passing filter: {len(unique_dc):,}", flush=True)

    # Build dict: dc_key -> list of cell indices (post-filter)
    dc_cells = {}
    for i in range(n_cells):
        if not pass_mask[i]:
            continue
        k = dc_keys[i]
        dc_cells.setdefault(k, []).append(i)

    # For each (donor, class), load those cells, aggregate per peak
    rows_out = []  # list of dicts for parquet output
    for dc_idx, dc_key in enumerate(unique_dc):
        donor_id, major_class = dc_key.split("|", 1)
        cell_idx = np.array(dc_cells[dc_key])
        n_cells_total = len(cell_idx)
        if n_cells_total == 0:
            continue
        # Fetch submatrix for these cells
        X_sub = a.X[cell_idx]
        if not sp.issparse(X_sub):
            X_sub = sp.csr_matrix(X_sub)
        # n_cells_access per peak = count of cells with X[cell, peak] > 0
        n_cells_access = np.asarray((X_sub > 0).sum(axis=0)).ravel()
        # Only emit peaks with n_cells_access >= 3 (matches add_surgical_agreement MIN_PER_PEAK)
        active = np.where(n_cells_access >= 3)[0]
        if len(active) == 0:
            continue
        for pi in active:
            rows_out.append({
                "peak_num": int(peak_nums[pi]),
                "major_class": major_class,
                "donor_id": donor_id,
                "n_cells_access_direct": int(n_cells_access[pi]),
                "n_cells_total_direct": int(n_cells_total),
                "fraction_open_direct": float(n_cells_access[pi]) / n_cells_total,
                "study_id": study_id,
            })
        if (dc_idx + 1) % 20 == 0 or dc_idx == len(unique_dc) - 1:
            print(f"  processed {dc_idx+1}/{len(unique_dc)} groups, {len(rows_out):,} output rows so far ({time.time()-t0:.0f}s)", flush=True)

    print(f"\n[4/4] Writing {args.output}", flush=True)
    out = pl.DataFrame(rows_out)
    out.write_parquet(args.output)
    print(f"  wrote {out.shape[0]:,} rows, {Path(args.output).stat().st_size/1_048_576:.1f} MB", flush=True)
    print(f"  Total elapsed: {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
