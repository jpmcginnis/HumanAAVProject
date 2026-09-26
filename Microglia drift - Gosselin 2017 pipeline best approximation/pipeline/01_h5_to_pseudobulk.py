#!/usr/bin/env python3
"""
Step 1 of the pipeline: 10x/Cell Ranger .h5 files → pseudobulk counts.

Given a set of snRNA-seq h5 files (one per donor-timepoint sample) and
a donor-pairing table, this script:

  1. Loads and concatenates the h5s.
  2. Runs standard scanpy QC (mito < 10 %, min_genes = 200,
     min_cells = 3 per gene).
  3. HVG selection, PCA, neighbors, leiden clustering at resolution 0.5.
  4. Identifies microglia by scoring each cell against a canonical
     marker panel (CX3CR1, P2RY12, TMEM119, C1QA, AIF1, CSF1R, ITGAM)
     using scanpy.tl.score_genes, then keeping any leiden cluster with
     mean score > 0.05. This multi-cluster call catches D14/D7 microglia
     that have drifted transcriptionally and clustered separately from
     D0 microglia — critical, because a naive top-cluster call misses
     them entirely.
  5. For each (donor × timepoint), sums raw counts across all
     microglial nuclei to produce one pseudobulk column.
  6. Writes pseudobulk counts, sample metadata, and a QC log.

Outputs match what Step 2 (`02_limma_voom.R`) expects:
  - pseudobulk_counts.csv  : gene × sample integer counts
  - sample_metadata.csv    : donor, condition (ex/iv) per sample
  - microglia_counts_per_sample.csv : nuclei counts for QC
  - qc_summary.txt         : per-step counts for auditability

Usage:
  python 01_h5_to_pseudobulk.py --config config.yaml --outdir output/

Config format (YAML) — see pipeline/config_example.yaml.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import yaml


# Canonical microglia marker panel used for cell-type identification.
MICROGLIA_MARKERS = ["CX3CR1", "P2RY12", "TMEM119", "C1QA", "AIF1",
                     "CSF1R", "ITGAM"]

# QC defaults.
QC_MIN_GENES = 200          # per-cell floor
QC_MIN_CELLS = 3            # per-gene floor
QC_MAX_PCT_MT = 10.0        # per-cell mito percent ceiling
LEIDEN_RESOLUTION = 0.5
N_PCA_COMPONENTS = 50
N_NEIGHBORS = 15
N_TOP_HVG = 2000
MICROGLIA_SCORE_THRESHOLD = 0.05    # cluster mean score to keep


def load_config(path: str) -> Dict:
    with open(path) as fh:
        return yaml.safe_load(fh)


def load_and_merge(pairing: pd.DataFrame) -> ad.AnnData:
    """Load each h5 in the pairing table, tag with donor+condition, concat."""
    parts = []
    for _, row in pairing.iterrows():
        for cond_label, path_key in [("ex", "d0_path"), ("iv", "d14_path")]:
            path = row[path_key]
            print(f"  loading {row['donor_id']} {cond_label}: "
                  f"{os.path.basename(path)}")
            a = sc.read_10x_h5(path)
            a.var_names_make_unique()
            a.obs["donor"] = row["donor_id"]
            a.obs["condition"] = cond_label
            a.obs["sample"] = f"{row['donor_id']}_{cond_label}"
            a.obs_names = [f"{row['donor_id']}_{cond_label}_{b}"
                            for b in a.obs_names]
            parts.append(a)
    merged = ad.concat(parts, axis=0, join="inner", merge="same",
                        index_unique=None)
    return merged


def run_qc(adata: ad.AnnData) -> ad.AnnData:
    """Standard scanpy QC. Preserves raw counts in layers['counts'] before
    normalization; scanpy .raw slot is set to the log-normalized data."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                log1p=False, inplace=True)
    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=QC_MIN_GENES)
    sc.pp.filter_genes(adata, min_cells=QC_MIN_CELLS)
    adata = adata[adata.obs["pct_counts_mt"] < QC_MAX_PCT_MT].copy()
    print(f"  cells: {n0} -> {adata.n_obs} after QC")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    return adata


def cluster_and_score_microglia(adata: ad.AnnData) -> ad.AnnData:
    """PCA → neighbors → leiden → microglia identity score → is_microglia."""
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_HVG, flavor="seurat")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=N_PCA_COMPONENTS, random_state=0)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=N_NEIGHBORS, random_state=0)
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, flavor="igraph",
                  n_iterations=2, directed=False, random_state=0)

    markers = [g for g in MICROGLIA_MARKERS if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=markers, score_name="microglia_score",
                       ctrl_size=50, n_bins=25, random_state=0)

    per_cluster = (adata.obs.groupby("leiden", observed=True)
                    .agg(n=("leiden", "size"),
                         microglia_score=("microglia_score", "mean"),
                         frac_ex=("condition", lambda s: (s == "ex").mean()))
                    .sort_values("microglia_score", ascending=False))
    print("  per-cluster microglia score (top 10):")
    print(per_cluster.head(10).to_string())

    keep = per_cluster[per_cluster["microglia_score"]
                       > MICROGLIA_SCORE_THRESHOLD].index.tolist()
    print(f"  -> microglia clusters (score > {MICROGLIA_SCORE_THRESHOLD}): "
          f"{keep}")
    adata.obs["is_microglia"] = adata.obs["leiden"].isin(keep)
    return adata


def pseudobulk_microglia(microglia: ad.AnnData, outdir: Path
                         ) -> pd.DataFrame:
    """Sum raw counts across nuclei per (donor, condition) → gene × sample."""
    from scipy.sparse import csr_matrix
    counts = microglia.layers["counts"]
    if not isinstance(counts, csr_matrix):
        counts = csr_matrix(counts)
    genes = microglia.var_names.to_numpy()
    samples = microglia.obs["sample"].to_numpy()
    order = sorted(pd.unique(samples))
    n_per_sample = {}
    pb = np.zeros((len(genes), len(order)), dtype=np.float64)
    for j, s in enumerate(order):
        mask = samples == s
        n_per_sample[s] = int(mask.sum())
        pb[:, j] = np.asarray(counts[mask].sum(axis=0)).ravel()
    df = pd.DataFrame(pb.astype(int), index=genes, columns=order)
    df.to_csv(outdir / "pseudobulk_counts.csv")
    # metadata
    meta_rows = []
    for s in order:
        donor, cond = s.rsplit("_", 1)
        meta_rows.append({"sample": s, "donor": donor, "condition": cond,
                           "n_microglia_nuclei": n_per_sample[s]})
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(outdir / "sample_metadata.csv", index=False)
    per_donor = (meta.pivot(index="donor", columns="condition",
                             values="n_microglia_nuclei")
                     .fillna(0).astype(int))
    per_donor.to_csv(outdir / "microglia_counts_per_sample.csv")
    print("  microglia nuclei per donor x condition:")
    print(per_donor.to_string())
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, help="YAML config path")
    ap.add_argument("--outdir", required=True, help="Output directory")
    args = ap.parse_args()

    cfg = load_config(args.config)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pairing = pd.DataFrame(cfg["donors"])
    required_cols = {"donor_id", "d0_path", "d14_path"}
    missing = required_cols - set(pairing.columns)
    if missing:
        sys.exit(f"config donors table missing columns: {missing}")

    print("Step 1a: load + merge")
    merged = load_and_merge(pairing)
    print(f"  merged: {merged.shape}")
    merged.write(outdir / "merged_raw.h5ad")

    print("\nStep 1b: QC + normalize")
    qc = run_qc(merged)
    qc.write(outdir / "merged_qc.h5ad")

    print("\nStep 1c: cluster + identify microglia")
    qc = cluster_and_score_microglia(qc)
    microglia = qc[qc.obs["is_microglia"]].copy()
    print(f"  microglia nuclei retained: {microglia.n_obs}")
    microglia.write(outdir / "microglia.h5ad")

    print("\nStep 1d: pseudobulk per (donor, condition)")
    pseudobulk_microglia(microglia, outdir)
    print(f"\nDone. Output in: {outdir}")


if __name__ == "__main__":
    main()
