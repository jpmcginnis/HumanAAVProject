"""
Pan-malignant matrix build.

Rebuilds the (peak × patient × cohort) accessibility aggregate specifically for
ALL tumor cells, bypassing Neftel-subtype classification. Produced because the
standard ``src/matrix_build.py`` chain in this pipeline grouped cells by
``obs['cell_type']``, which was assigned by marker-peak scoring — a classifier
that often MIS-LABELS CNV-malignant cells as TAM / microglia / neuron / OPC
when their chromatin incidentally opens at a TME marker promoter.

The symptom: in the main `enhancer_candidate_matrix.parquet`, the
``malignant_unresolved`` pool max-cohort coverage was 3/8 cohorts and 21/51
patients. Not enough replication for an AAV cloning shortlist. Looking at the
per-cohort h5ads directly: GBM-Space had 438,994 cells with
``malignant_cnv == 1`` (the CNV-based tumor cell call), but ZERO of them
landed under the ``malignant_unresolved`` label — all were mis-labeled by
marker-peak scoring to TME types.

This module fixes that by aggregating along the CNV flag (plus a label
fallback for cohorts where CNV wasn't computed), yielding ~2.7× more
pan-malignant cells and up to 6/8 cohort coverage per peak.

Inputs
------
A list of ``catlas_quantified.h5ad`` files (one per cohort). Each file's
``obs`` must have:

- ``cohort`` (str), ``patient_id`` (str), ``cell_type`` (str)
- ``malignant_cnv`` (int 0/1) when CNV calling was run; else absent

A cell is included in the pan-malignant pool if
``malignant_cnv == 1 OR cell_type == "malignant_unresolved"``.

Outputs
-------
Two parquet files:

- ``pan_malignant_matrix.parquet`` — long-format rows at the
  ``(peak_id, patient_id, cohort)`` level with ``n_cells_malignant``,
  ``n_accessible``, ``frac_accessible``, ``mean_counts``.
- ``pan_malignant_scores.parquet`` — per-peak summary with ``strength_mean``,
  ``strength_median``, ``n_patients_accessible_10pct``, ``n_patients_detected``,
  ``n_cohorts``, ``consistency`` = fraction of detected patients ≥10%
  accessible.

Usage
-----
::

    python -u src/pan_malignant_matrix.py \\
        --inputs /data/processed/*/catlas_quantified.h5ad \\
        --matrix matrix/pan_malignant_matrix.parquet \\
        --scores matrix/pan_malignant_scores.parquet \\
        --consistency-threshold 0.1

Caveats
-------
This explicitly includes CNV-called malignant cells that may have been
scored as TME cells by marker-peak scoring. For fragment-mode cohorts whose
``src/flatten_cohort.py`` run did not populate ``malignant_cnv`` (currently
``tcga_scatac`` and ``gse276177_khan_astro``), only cells already labeled
``malignant_unresolved`` are included — i.e. the cohort contributes less than
it should. Running CNV calling in the flatten path (TODO) would recover
these.
"""
from __future__ import annotations
import argparse
import re
import sys
from pathlib import Path


def log(msg: str) -> None:
    print(f"[pan] {msg}", flush=True)


def main() -> int:
    import anndata as ad
    import numpy as np
    import pandas as pd
    import polars as pl

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--inputs", nargs="+", required=True, type=Path,
        help="Per-cohort catlas_quantified.h5ad files (one per cohort).",
    )
    ap.add_argument("--matrix", type=Path, required=True, help="Output long-format parquet")
    ap.add_argument("--scores", type=Path, required=True, help="Output per-peak scores parquet")
    ap.add_argument(
        "--consistency-threshold", type=float, default=0.1,
        help="A patient 'counts' toward consistency if frac_accessible >= this (default 0.1 = 10%% of malignant cells open)",
    )
    ap.add_argument(
        "--min-cells-per-patient", type=int, default=10,
        help="Skip patients with fewer than this many malignant cells",
    )
    args = ap.parse_args()

    all_dfs = []

    for src in args.inputs:
        log(f"=== {src.parent.name} ===")
        a = ad.read_h5ad(src, backed="r")
        cohort_name = str(a.obs["cohort"].iloc[0])
        modality = str(a.obs["modality"].iloc[0]) if "modality" in a.obs.columns else "unknown"
        obs = a.obs

        mal_cnv = (
            obs["malignant_cnv"].astype(int).values if "malignant_cnv" in obs.columns
            else np.zeros(a.n_obs, dtype=int)
        )
        mal_lab = obs["cell_type"].astype(str).values == "malignant_unresolved"
        pan_mal = (mal_cnv == 1) | mal_lab
        log(f"  pan-malignant cells: {int(pan_mal.sum()):,} / {a.n_obs:,}  "
            f"(cnv={int((mal_cnv==1).sum()):,}, label={int(mal_lab.sum()):,})")

        if pan_mal.sum() == 0:
            log(f"  SKIP {cohort_name} — no pan-malignant cells")
            try: a.file.close()
            except Exception: pass
            continue

        # Precompute var side (chrom/start/end + peak_id)
        n_peaks = a.n_vars
        if "chrom" in a.var.columns:
            chrom_all = a.var["chrom"].astype(str).values
            start_all = a.var["start"].astype(np.int64).values
            end_all = a.var["end"].astype(np.int64).values
        else:
            var_names = np.asarray(a.var_names, dtype=object)
            chrom_all = np.empty(n_peaks, dtype=object)
            start_all = np.zeros(n_peaks, dtype=np.int64)
            end_all = np.zeros(n_peaks, dtype=np.int64)
            pat = re.compile(r"^(chr[\dXYM]+):(\d+)-(\d+)$")
            for pi, name in enumerate(var_names):
                m = pat.match(str(name))
                if m:
                    chrom_all[pi] = m.group(1)
                    start_all[pi] = int(m.group(2))
                    end_all[pi] = int(m.group(3))
        peak_id_all = np.array(
            [f"{c}:{s}-{e}" for c, s, e in zip(chrom_all, start_all, end_all)],
            dtype=object,
        )

        # Per-patient aggregation (ignore region — fold over regions)
        for pid, pid_df in obs[pan_mal].groupby("patient_id", observed=True):
            idx = pid_df.index
            n_cells = len(idx)
            if n_cells < args.min_cells_per_patient:
                continue
            mat = a[idx, :].X
            n_access = np.asarray((mat > 0).sum(axis=0)).flatten()
            mean_cnt = np.asarray(mat.mean(axis=0)).flatten()
            keep = n_access > 0
            k = int(keep.sum())
            if k == 0:
                continue

            df_group = pl.DataFrame({
                "peak_id":          peak_id_all[keep].astype(str),
                "chrom":            chrom_all[keep].astype(str),
                "start":            start_all[keep].astype(np.int64),
                "end":              end_all[keep].astype(np.int64),
                "cohort":           [cohort_name] * k,
                "patient_id":       [str(pid)] * k,
                "modality":         [modality] * k,
                "n_cells_malignant": np.full(k, n_cells, dtype=np.int64),
                "n_accessible":     n_access[keep].astype(np.int64),
                "frac_accessible":  (n_access[keep].astype(np.float64) / n_cells).astype(np.float32),
                "mean_counts":      mean_cnt[keep].astype(np.float32),
            })
            all_dfs.append(df_group)
            log(f"  patient {pid}: {n_cells:,} malignant cells, {k:,} accessible peaks")

        try: a.file.close()
        except Exception: pass

    log(f"\nConcat {len(all_dfs)} per-patient DataFrames")
    big = pl.concat(all_dfs, how="vertical", rechunk=True)
    log(f"pan-malignant matrix: {big.height:,} rows x {big.width} cols")

    args.matrix.parent.mkdir(parents=True, exist_ok=True)
    big.write_parquet(args.matrix, compression="zstd")
    log(f"wrote {args.matrix} ({args.matrix.stat().st_size / 1e6:.1f} MB)")

    # Point-estimate scoring per peak
    log("\n=== per-peak scoring ===")
    scores = (
        big.group_by("peak_id")
        .agg([
            pl.col("chrom").first(),
            pl.col("start").first(),
            pl.col("end").first(),
            pl.col("frac_accessible").mean().alias("strength_mean"),
            pl.col("frac_accessible").median().alias("strength_median"),
            (pl.col("frac_accessible") >= args.consistency_threshold).sum().alias("n_patients_accessible"),
            pl.col("patient_id").n_unique().alias("n_patients_detected"),
            pl.col("cohort").n_unique().alias("n_cohorts"),
        ])
        .with_columns(
            (pl.col("n_patients_accessible").cast(pl.Float64) / pl.col("n_patients_detected"))
            .alias("consistency")
        )
        .sort("strength_mean", descending=True)
    )
    log(f"scored peaks: {scores.height:,}")

    args.scores.parent.mkdir(parents=True, exist_ok=True)
    scores.write_parquet(args.scores, compression="zstd")
    log(f"wrote {args.scores} ({args.scores.stat().st_size / 1e6:.1f} MB)")

    log("\n=== TOP 30 pan-malignant peaks (n_cohorts >= 3, sorted by n_cohorts then n_patients then strength) ===")
    top = (
        scores
        .filter(pl.col("n_cohorts") >= 3)
        .sort(
            ["n_cohorts", "n_patients_detected", "strength_mean"],
            descending=[True, True, True],
        )
        .head(30)
    )
    print(top.to_pandas().to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
