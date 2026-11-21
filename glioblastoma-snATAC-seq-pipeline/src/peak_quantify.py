"""
Quantify every cohort against the CATLAS 544K peak atlas (common coordinate system).

Two paths:

1. Fragment-mode (small cohorts): snap.pp.make_peak_matrix on the fragments.
   Loads adata unbacked. Only used for cohorts whose h5ad carries
   fragments_bed_path in .uns AND is small enough to fit in RAM.

2. Peak-matrix lift (huge pre-matrixed cohorts, e.g. GBM-Space 1M cells × 438K
   author peaks): lift X from the author's peak coordinates onto CATLAS peak
   coordinates via pyranges reciprocal overlap. Streamed chunk-by-chunk so
   RAM stays bounded. Uses the SAME 25% reciprocal overlap rule as
   src/ingest_matrix.py so the resulting matrix is directly compatible with
   other matrix-mode cohorts in the atlas.
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np


def log(msg): print(f"[quantify] {msg}", flush=True)


def _load_catlas_peaks(bed: Path):
    """Return (DataFrame with chr/start/end/catlas_idx, pyranges)."""
    import pandas as pd
    import pyranges as pr
    df = pd.read_csv(bed, sep="\t", header=None, names=["Chromosome", "Start", "End"])
    df["catlas_idx"] = np.arange(len(df), dtype=np.int64)
    log(f"  CATLAS peaks: n={len(df)}")
    return df, pr.PyRanges(df)


def _var_to_pyranges(adata):
    """Parse chrom/start/end from adata.var into a pyranges. Handles common
    column name variants (chrom/seqnames/chr/Chromosome) and falls back to
    parsing var_names like 'chr1:1234-5678' OR 'chr1:1234:5678:...'."""
    import pandas as pd
    import pyranges as pr

    _CHROM_COLS = ("chrom", "seqnames", "chr", "Chromosome")
    _START_COLS = ("start", "Start", "chromStart")
    _END_COLS = ("end", "End", "chromEnd")
    v = adata.var
    chrom_col = next((c for c in _CHROM_COLS if c in v.columns), None)
    start_col = next((c for c in _START_COLS if c in v.columns), None)
    end_col = next((c for c in _END_COLS if c in v.columns), None)
    if chrom_col and start_col and end_col:
        df = pd.DataFrame({
            "Chromosome": v[chrom_col].astype(str).values,
            "Start": pd.to_numeric(v[start_col], errors="coerce").values,
            "End": pd.to_numeric(v[end_col], errors="coerce").values,
        })
    else:
        parsed = v.index.to_series().str.extract(
            r"^(?P<Chromosome>chr[\dXYMT]+):(?P<Start>\d+)[-:](?P<End>\d+)"
        )
        df = pd.DataFrame({
            "Chromosome": parsed["Chromosome"].astype(str).values,
            "Start": pd.to_numeric(parsed["Start"], errors="coerce").values,
            "End": pd.to_numeric(parsed["End"], errors="coerce").values,
        })
    # Drop any unparseable rows
    valid = df["Chromosome"].notna() & df["Start"].notna() & df["End"].notna()
    if not valid.all():
        log(f"  WARNING: {int((~valid).sum())}/{len(df)} var rows had unparseable coords; dropping")
    df["src_idx"] = np.arange(len(df), dtype=np.int64)
    df = df[valid].reset_index(drop=True)
    df["Start"] = df["Start"].astype(int)
    df["End"] = df["End"].astype(int)
    return df, pr.PyRanges(df)


def _build_lift_matrix(src_df, catlas_df, catlas_pr, overlap_frac=0.25):
    """Return a sparse (n_src, n_catlas) matrix M such that
    X_catlas = X_src @ M. Each src peak with reciprocal 25% overlap into a
    CATLAS peak contributes 1.0 (summed across many src peaks per CATLAS
    peak). Uses pyranges' join for the overlap compute."""
    import pyranges as pr
    import scipy.sparse as sp
    src_pr = pr.PyRanges(src_df.rename(columns={"src_idx": "SrcIdx"}))
    cat_pr = pr.PyRanges(catlas_df.rename(columns={"catlas_idx": "CatIdx"}))
    joined = src_pr.join(cat_pr).df
    if joined.empty:
        return sp.csr_matrix((len(src_df), len(catlas_df)), dtype=np.float32)
    inter_lo = np.maximum(joined["Start"].values, joined["Start_b"].values)
    inter_hi = np.minimum(joined["End"].values, joined["End_b"].values)
    inter = np.maximum(inter_hi - inter_lo, 0)
    src_len = joined["End"].values - joined["Start"].values
    cat_len = joined["End_b"].values - joined["Start_b"].values
    recip_mask = (inter >= overlap_frac * src_len) & (inter >= overlap_frac * cat_len)
    joined = joined[recip_mask]
    log(f"  overlaps: {len(joined):,} reciprocal-{int(overlap_frac*100)}% src→CATLAS pairs")
    rows = joined["SrcIdx"].values.astype(np.int64)
    cols = joined["CatIdx"].values.astype(np.int64)
    data = np.ones(len(joined), dtype=np.float32)
    M = sp.csr_matrix((data, (rows, cols)),
                      shape=(len(src_df), len(catlas_df)))
    recall = int((M.sum(axis=0) > 0).sum()) / len(catlas_df)
    log(f"  CATLAS peak recall after lift: {recall:.2%}")
    return M


def _lift_backed(adata, M, out_path, chunk=50_000):
    """Chunked lift that streams BOTH compute and output to disk — avoids
    the 100+ GB in-memory accumulation that OOM'd the box.

    Strategy: for each cell chunk, compute X_catlas_chunk = X_src_chunk @ M
    and write it to its own temp .h5ad. At end, build an on-disk CSR by
    opening each temp with h5py and concatenating indptr/indices/data arrays
    into pre-sized datasets in the output. Only one chunk at a time lives
    in RAM."""
    import scipy.sparse as sp
    import h5py as _h5
    try:
        from anndata.io import write_elem as _write_elem
    except ImportError:
        from anndata.experimental import write_elem as _write_elem

    n_obs = adata.n_obs
    n_catlas = M.shape[1]
    log(f"  streaming lift: {n_obs:,} cells × {M.shape[0]:,} src → {n_catlas:,} CATLAS")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_path.parent / f".{out_path.stem}_chunks"
    tmp_dir.mkdir(exist_ok=True)

    # Pass 1: compute each chunk's lifted X, write to its own .npz on disk.
    # Memory peak: one chunk's sparse matrix (~1-2 GB).
    chunk_meta = []  # list of (path, nnz, n_rows_in_chunk)
    for i in range(0, n_obs, chunk):
        j = min(i + chunk, n_obs)
        Xc = adata.X[i:j, :]
        if not sp.issparse(Xc):
            Xc = sp.csr_matrix(Xc)
        X_cat = (Xc @ M).tocsr().astype(np.float32)
        tmp_npz = tmp_dir / f"chunk_{i:08d}.npz"
        sp.save_npz(tmp_npz, X_cat, compressed=False)
        chunk_meta.append((tmp_npz, int(X_cat.nnz), int(X_cat.shape[0])))
        if (i // chunk) % 5 == 0:
            log(f"    chunk {i:,}/{n_obs:,}  (nnz={X_cat.nnz:,})")
        del Xc, X_cat
    total_nnz = sum(m[1] for m in chunk_meta)
    log(f"  all chunks on disk. total nnz = {total_nnz:,}  (expect ~{total_nnz*8/1e9:.1f} GB final)")

    # Pass 2: build output h5ad with sparse CSR groups. Create preallocated
    # indptr/indices/data arrays and fill from the per-chunk npz files.
    log(f"  writing output h5ad to {out_path} (one-pass stream)")
    import pandas as pd
    import anndata as ad
    # Parse CATLAS peak coordinates from the BED file so var carries real
    # chrom/start/end (needed by matrix_build and scoring). h5py can't write
    # None, so these must be concrete strings/ints.
    import pandas as pd
    catlas_bed_df = pd.read_csv(
        "/data/projects/atacseq/ref/catlas_peaks.bed",
        sep="\t", header=None, names=["chrom", "start", "end"]
    )
    assert len(catlas_bed_df) == n_catlas, f"CATLAS BED has {len(catlas_bed_df)} but n_catlas={n_catlas}"
    out_var = catlas_bed_df.copy()
    out_var["chrom"] = out_var["chrom"].astype(str)
    out_var["start"] = out_var["start"].astype(int)
    out_var["end"] = out_var["end"].astype(int)
    out_var.index = [f"{c}:{s}-{e}" for c, s, e in zip(out_var["chrom"], out_var["start"], out_var["end"])]

    with _h5.File(out_path, "w") as hf:
        # Write obs/var via anndata's elem writer so an.read_h5ad roundtrips cleanly
        obs_copy = adata.obs.copy()
        _write_elem(hf, "obs", obs_copy)
        _write_elem(hf, "var", out_var)
        # Write X as a CSR group manually, filling datasets chunk-by-chunk
        xg = hf.create_group("X")
        xg.attrs["encoding-type"] = "csr_matrix"
        xg.attrs["encoding-version"] = "0.1.0"
        xg.attrs["shape"] = np.array([n_obs, n_catlas], dtype=np.int64)
        d_data = xg.create_dataset("data", shape=(total_nnz,), dtype=np.float32,
                                   chunks=True)
        d_idx = xg.create_dataset("indices", shape=(total_nnz,), dtype=np.int32,
                                  chunks=True)
        d_iptr = xg.create_dataset("indptr", shape=(n_obs + 1,), dtype=np.int64,
                                   chunks=True)
        d_iptr[0] = 0
        row_off = 0
        nnz_off = 0
        for k, (npz_path, nnz_k, nrows_k) in enumerate(chunk_meta):
            Xk = sp.load_npz(npz_path)
            # append data/indices
            d_data[nnz_off:nnz_off + nnz_k] = Xk.data.astype(np.float32)
            d_idx[nnz_off:nnz_off + nnz_k] = Xk.indices.astype(np.int32)
            # indptr: this chunk's indptr is 0..nnz_k; shift by nnz_off, drop first
            d_iptr[row_off + 1:row_off + 1 + nrows_k] = Xk.indptr[1:].astype(np.int64) + nnz_off
            row_off += nrows_k
            nnz_off += nnz_k
            del Xk
            if k % 5 == 0:
                log(f"    merged chunk {k+1}/{len(chunk_meta)}")

    # Cleanup temp chunks
    for npz_path, _, _ in chunk_meta:
        try: npz_path.unlink()
        except Exception: pass
    try: tmp_dir.rmdir()
    except Exception: pass
    log(f"  wrote {out_path} ({out_path.stat().st_size/1e9:.1f} GB)")


def main() -> int:
    import anndata as ad
    import snapatac2 as snap

    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--peaks", required=True, type=Path, help="CATLAS peak BED file")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=16)
    args = ap.parse_args()

    file_gb = args.src.stat().st_size / 1e9
    use_backed = file_gb > 10
    log(f"Loading {args.src} (size {file_gb:.1f} GB, backed={use_backed})")
    adata = ad.read_h5ad(args.src, backed="r" if use_backed else None)
    log(f"  n_cells: {adata.n_obs}  n_vars: {adata.n_vars}")

    has_fragments = bool(adata.uns.get("fragments_bed_path"))

    if use_backed or not has_fragments:
        # Peak-matrix lift path (chunked). Works whether backed or not.
        log("Quantification via CATLAS peak lift (no fragments / huge query)")
        catlas_df, catlas_pr = _load_catlas_peaks(args.peaks)
        src_df, _ = _var_to_pyranges(adata)
        M = _build_lift_matrix(src_df, catlas_df, catlas_pr, overlap_frac=0.25)
        _lift_backed(adata, M, args.out)
        if use_backed:
            adata.file.close()
        return 0

    # Fragment-mode path: use snap.pp.make_peak_matrix (loads unbacked)
    log("Quantification via snap.pp.make_peak_matrix (fragment-mode)")
    snap.pp.make_peak_matrix(
        adata,
        peak_file=str(args.peaks),
        n_jobs=args.threads,
    )
    log(f"  n_peaks after requantification: {adata.n_vars}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    log(f"Wrote → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
