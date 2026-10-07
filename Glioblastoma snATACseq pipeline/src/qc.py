"""
Per-cohort QC. Thresholds from config/pipeline.yaml.

Use plain `anndata.read_h5ad` always (polars-backed snap.read panics on obs
indexing with heterogeneous dtypes). Fragment-based metrics (TSSe, frag size
distribution) only run when the cohort was imported from fragments by
snapatac2, which we detect via `uns['fragments_bed_path']`.
"""

from __future__ import annotations
import argparse
import sys
import yaml
from pathlib import Path


def log(msg): print(f"[qc] {msg}", flush=True)


def _has_fragments(adata) -> bool:
    return bool(adata.uns.get("fragments_bed_path"))


def main() -> int:
    import anndata as ad
    import numpy as np
    import pandas as pd

    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--report", required=True, type=Path)
    args = ap.parse_args()

    with open("config/pipeline.yaml") as f:
        cfg = yaml.safe_load(f)["qc"]

    # Use backed='r' for huge files (GBM-Space is 127 GB uncompressed) — only
    # load slices into memory on demand.
    file_size_gb = args.src.stat().st_size / 1e9
    use_backed = file_size_gb > 20
    log(f"Loading via anndata (file size {file_size_gb:.1f} GB, backed={use_backed}): {args.src}")
    adata = ad.read_h5ad(args.src, backed="r" if use_backed else None)
    log(f"  n_cells before QC: {adata.n_obs}")
    if "sample_id" in adata.obs.columns:
        log(f"  n_samples: {adata.obs['sample_id'].nunique()}")
    if "patient_id" in adata.obs.columns:
        log(f"  n_patients: {adata.obs['patient_id'].nunique()}")

    has_frags = _has_fragments(adata)
    log(f"  fragment-based QC available: {has_frags}")

    if has_frags:
        # Hand off to snapatac2 for TSSe + fragment-size + fragment-count filtering
        import snapatac2 as snap
        log("Fragment-mode QC via snapatac2 (TSSe + fragment size distribution)")
        # Convert to backed AnnDataSet so snapatac2 can operate; use a scratch path
        scratch = Path("/data/projects/atacseq/processed") / args.src.parent.name / "_snap_qc_scratch.h5ad"
        scratch.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(scratch)
        del adata
        adata = snap.read(scratch, backed="r+")
        snap.metrics.tsse(adata, snap.genome.hg38)
        snap.metrics.frag_size_distr(adata)
        log(f"Filtering: TSSe ≥ {cfg['min_tss_enrichment']}, "
            f"{cfg['min_fragments_per_cell']} ≤ n_fragment ≤ {cfg['max_fragments_per_cell']}")
        snap.pp.filter_cells(
            adata,
            min_counts=cfg["min_fragments_per_cell"],
            max_counts=cfg["max_fragments_per_cell"],
            min_tsse=cfg["min_tss_enrichment"],
        )
        log(f"  n_cells after QC: {adata.n_obs}")
        adata.write(args.out)
        try:
            scratch.unlink()
        except Exception:
            pass
    else:
        log(f"Matrix/annotated-mode cohort: count-based QC only")
        log(f"  cells with {cfg['min_fragments_per_cell']} ≤ counts ≤ {cfg['max_fragments_per_cell']}")

        # Prefer pre-computed count columns if present (GBM-Space ships 'nFrags',
        # ArchR calls it that; others may use 'n_fragment' / 'total_counts').
        totals = None
        for col in ("nFrags", "n_fragment", "total_counts", "n_counts", "ReadsInPeaks"):
            if col in adata.obs.columns:
                totals = np.asarray(adata.obs[col]).astype(float)
                log(f"  using precomputed obs['{col}'] for count filter "
                    f"(median={np.median(totals):.0f})")
                break
        if totals is None:
            # Fallback: compute from X — streaming if backed, else in memory
            if use_backed and hasattr(adata.X, "nnz"):
                log("  summing backed sparse X in chunks (10K cells per chunk)")
                totals = np.zeros(adata.n_obs, dtype=np.float64)
                chunk = 10000
                for i in range(0, adata.n_obs, chunk):
                    j = min(i + chunk, adata.n_obs)
                    block = adata.X[i:j]
                    totals[i:j] = np.asarray(block.sum(axis=1)).flatten()
            else:
                totals = np.asarray(adata.X.sum(axis=1)).flatten()
            log(f"  computed totals (median={np.median(totals):.0f})")

        keep = (totals >= cfg["min_fragments_per_cell"]) & (totals <= cfg["max_fragments_per_cell"])
        log(f"  cells before: {adata.n_obs}, keeping: {int(keep.sum())}")

        # Slice to kept cells. On backed data we materialize into memory ONLY the
        # kept cells; on in-memory data this is just a view.
        keep_idx = np.where(keep)[0]
        if use_backed:
            args.out.parent.mkdir(parents=True, exist_ok=True)
            frac_keep = len(keep_idx) / adata.n_obs
            if frac_keep >= 0.9999:
                # 100% (or near-100%) pass → hardlink input to output, patch obs.
                # Avoids materializing a 100+ GB sparse matrix in RAM (what
                # OOM'd the box on the 1M-cell GBM-Space cohort).
                import os as _os
                import h5py as _h5
                try:
                    from anndata.io import write_elem as _write_elem
                except ImportError:
                    from anndata.experimental import write_elem as _write_elem
                log(f"  100% of cells pass QC — using hardlink trick (zero RAM/disk cost)")
                src_path = str(adata.filename)
                adata.file.close()  # must close before hardlinking the backing file
                del adata
                if args.out.exists():
                    args.out.unlink()
                try:
                    _os.link(src_path, args.out)
                    log(f"  hardlinked {src_path} → {args.out}")
                except OSError as e:
                    log(f"  hardlink failed ({e}); falling back to shutil.copy")
                    import shutil as _sh
                    _sh.copy(src_path, args.out)
                # Patch obs to add qc_count_total column
                import anndata as ad_mod
                _tmp = ad_mod.read_h5ad(args.out, backed="r")
                _tmp_obs = _tmp.obs.copy()
                _tmp.file.close()
                _tmp_obs["qc_count_total"] = totals
                with _h5.File(args.out, "r+") as hf:
                    if "obs" in hf:
                        del hf["obs"]
                    _write_elem(hf, "obs", _tmp_obs)
                log(f"  wrote {args.out}  ({args.out.stat().st_size/1e9:.1f} GB on disk)")
                # Reload as backed for the report step below
                adata = ad_mod.read_h5ad(args.out, backed="r")
            else:
                # Partial filter: do actually materialize the kept subset.
                # IMPORTANT: close the backed input BEFORE building the subset
                # so we don't hold two giant h5py handles open simultaneously.
                log(f"  materializing kept subset in RAM (sparse) — {frac_keep:.1%} of cells kept")
                sub = adata[keep_idx].to_memory()
                sub.obs["qc_count_total"] = totals[keep_idx]
                adata.file.close()
                del adata
                log(f"  writing filtered h5ad (gzip)")
                sub.write_h5ad(args.out, compression="gzip", compression_opts=4)
                log(f"  wrote {args.out}  ({args.out.stat().st_size/1e9:.1f} GB)")
                adata = sub
        else:
            adata = adata[keep_idx].copy()
            adata.obs["qc_count_total"] = totals[keep_idx]
            log(f"  n_cells after QC: {adata.n_obs}")
            args.out.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(args.out)

    # Per-sample report, works on either path
    obs_pd = adata.obs if isinstance(adata.obs, pd.DataFrame) else adata.obs.to_pandas() \
        if hasattr(adata.obs, "to_pandas") else pd.DataFrame(dict(adata.obs))
    if "sample_id" in obs_pd.columns:
        per_sample = obs_pd.groupby("sample_id").size().reset_index(name="n_cells")
    else:
        per_sample = pd.DataFrame({"sample_id": ["(unknown)"], "n_cells": [adata.n_obs]})
    log("Per-sample cell counts post-QC:")
    for _, row in per_sample.iterrows():
        log(f"    {str(row['sample_id'])[:40]:40s} {int(row['n_cells']):>7d}")

    args.report.parent.mkdir(parents=True, exist_ok=True)
    per_sample.to_html(args.report, index=False)
    log(f"Wrote QC summary → {args.report}")
    log(f"Wrote filtered h5ad → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
