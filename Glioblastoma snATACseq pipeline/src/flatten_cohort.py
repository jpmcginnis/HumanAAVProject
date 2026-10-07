"""
Flatten per-sample snapatac2 h5ads into one cohort-level ``catlas_quantified.h5ad``
in CATLAS peak space, with marker-peak-based cell type labels.

Why this script exists
----------------------
``src/ingest.py`` imports fragments.tsv.gz files via ``snapatac2.pp.import_fragments``
and bundles them into an ``snapatac2.AnnDataSet`` written to ``ingested.h5ad``.
AnnDataSet is a *lazy pointer bundle* — ``.X`` is unresolved until peaks are quantified.
Downstream ``src/qc.py`` loads with vanilla ``anndata.read_h5ad`` and sees ``n_vars=0``,
so the QC → CNV → labels → quantify chain breaks for fragments-mode cohorts.

This module replaces that broken chain with a single pass that:

1. Reads each per-sample snapatac2 h5ad (already produced by ``src/ingest.py``'s
   ``_import_fragments_set``).
2. Projects each sample onto the CATLAS 544 K-peak BED via
   ``snapatac2.pp.make_peak_matrix`` → sparse per-sample peak matrix.
3. Builds a plain ``anndata.AnnData`` per sample with harmonized obs columns.
4. ``anndata.concat`` across samples (outer join on peaks — identical since all
   project onto the same CATLAS peaks).
5. Scores cells against marker TSSes (``src/marker_peak_scoring.py``) to produce
   ``cell_type`` labels.
6. Writes the result as ``catlas_quantified.h5ad`` — the input shape
   ``src/matrix_build.py`` expects.

Usage
-----
::

    python -u src/flatten_cohort.py \\
        --per-sample-glob "/data/processed/<cohort>/sample_*.h5ad" \\
        --cohort tcga_scatac \\
        --catlas-bed /data/reference/catlas_peaks.bed \\
        --markers-tsv config/marker_tsses.tsv \\
        --out /data/processed/<cohort>/catlas_quantified.h5ad

Per-sample patient_id, region, modality, idh_status are either read from each
h5ad's obs (if ``src/ingest.py`` wrote them) or defaulted via ``--default-*`` flags.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[flatten] {msg}", flush=True)


def main() -> int:
    import hdf5plugin  # noqa: F401 — register HDF5 filters snapatac2 uses
    import snapatac2 as snap
    import anndata as ad
    import pandas as pd
    import numpy as np
    import scipy.sparse as sp

    # marker_peak_scoring is in the same src/ package
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from marker_peak_scoring import score_cells_matrix_mode

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--per-sample-glob",
        required=True,
        help="Glob pattern matching per-sample snapatac2 h5ads. "
        "Must NOT match _import.h5ad or _peakmat.h5ad helper files; use a narrow pattern.",
    )
    ap.add_argument("--cohort", required=True, help="Cohort slug (goes into obs['cohort'])")
    ap.add_argument("--catlas-bed", type=Path, required=True, help="CATLAS peaks BED (chrom start end)")
    ap.add_argument("--markers-tsv", type=Path, required=True, help="Marker TSS TSV for cell type scoring")
    ap.add_argument("--out", type=Path, required=True, help="Output catlas_quantified.h5ad path")
    ap.add_argument("--default-region", default="tumor")
    ap.add_argument("--default-modality", default="snATAC")
    ap.add_argument("--default-idh-status", default="wildtype")
    ap.add_argument("--default-primary-recurrent", default="primary")
    ap.add_argument(
        "--patient-id-source",
        choices=["sample_id", "obs"],
        default="sample_id",
        help="'sample_id' (default) uses the h5ad's filename stem; 'obs' uses obs['patient_id'] if present",
    )
    args = ap.parse_args()

    # Resolve per-sample files from the glob. Exclude helper suffixes defensively.
    glob_root = Path(args.per_sample_glob).parent
    glob_pat = Path(args.per_sample_glob).name
    per_sample = sorted([
        p for p in glob_root.glob(glob_pat)
        if "_peakmat" not in p.name and "_import" not in p.name
    ])
    log(f"found {len(per_sample)} per-sample h5ads matching {args.per_sample_glob}")
    if len(per_sample) == 0:
        log("ABORT: no per-sample h5ads found. Check --per-sample-glob.")
        return 2

    # Load CATLAS peak BED once — all per-sample peak matrices share this var.
    cat_df = pd.read_csv(args.catlas_bed, sep="\t", header=None, names=["chrom", "start", "end"])
    cat_df["chrom"] = cat_df["chrom"].astype(str)
    cat_df["start"] = cat_df["start"].astype(int)
    cat_df["end"] = cat_df["end"].astype(int)
    cat_df.index = [f"{c}:{s}-{e}" for c, s, e in zip(cat_df["chrom"], cat_df["start"], cat_df["end"])]
    log(f"CATLAS peaks: {len(cat_df):,}")

    adatas = []
    for p in per_sample:
        sample_id = p.stem
        log(f"processing {sample_id}")

        a = snap.read(str(p), backed=None)

        # Patient id: default to sample_id; or read from per-sample h5ad obs when requested.
        patient_id = sample_id
        if args.patient_id_source == "obs":
            try:
                vals = a.obs["patient_id"]
                patient_id = str(vals[0]) if len(vals) else sample_id
            except (KeyError, AttributeError):
                pass

        pm_file = p.with_name(p.stem + "_peakmat.h5ad")
        pm_file.unlink(missing_ok=True)
        pm = snap.pp.make_peak_matrix(a, peak_file=str(args.catlas_bed), file=str(pm_file))

        X = pm.X[:]
        if not sp.issparse(X):
            X = sp.csr_matrix(X)
        n_cells, n_peaks = X.shape
        if n_peaks != len(cat_df):
            raise RuntimeError(
                f"{sample_id}: make_peak_matrix produced {n_peaks} peaks, "
                f"CATLAS BED has {len(cat_df)}. BED mismatch."
            )
        log(f"  n_obs={n_cells}  n_vars={n_peaks}  nnz={X.nnz:,}")

        barcodes = [str(b) for b in pm.obs_names]
        obs = pd.DataFrame(
            {
                "sample_id": sample_id,
                "patient_id": patient_id,
                "region": args.default_region,
                "cohort": args.cohort,
                "modality": args.default_modality,
                "idh_status": args.default_idh_status,
                "primary_recurrent": args.default_primary_recurrent,
            },
            index=barcodes,
        )
        a_plain = ad.AnnData(X=X, obs=obs, var=cat_df.copy())
        adatas.append(a_plain)

        try: pm.close()
        except Exception: pass
        try: a.close()
        except Exception: pass

    log(f"concatenating {len(adatas)} per-sample AnnDatas")
    merged = ad.concat(adatas, join="outer", index_unique="-")
    log(f"merged: n_obs={merged.n_obs:,}  n_vars={merged.n_vars:,}")

    # Cell type labels via marker-peak scoring
    log("running marker-peak scoring")
    score_cells_matrix_mode(merged, args.markers_tsv)
    if "label_tme_gbmap" in merged.obs.columns:
        merged.obs["cell_type"] = merged.obs["label_tme_gbmap"].astype(str)

    # Downstream-compat defaults (matrix_build / scoring expect these)
    merged.obs.setdefault = None  # pandas DataFrame has no setdefault; use conditional set
    for col, default in [("neftel_state", "unresolved"),
                         ("neuron_subtype", "non_neuron"),
                         ("malignant_cnv", 0)]:
        if col not in merged.obs.columns:
            merged.obs[col] = default

    log("cell_type distribution:")
    for ct, n in merged.obs["cell_type"].value_counts().head(15).items():
        log(f"    {str(ct):22s} {int(n):>7d}")

    log("patient distribution:")
    for pid, n in merged.obs["patient_id"].value_counts().items():
        log(f"    {str(pid):60s} {int(n):>7d}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.unlink(missing_ok=True)
    merged.write_h5ad(args.out, compression="lzf")
    log(f"WROTE {args.out} ({args.out.stat().st_size / 1e9:.2f} GB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
