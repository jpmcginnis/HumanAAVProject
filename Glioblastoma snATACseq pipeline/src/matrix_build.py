"""
Build the persistent (peak × cell_type × cohort × patient) matrix.
This is JP's artifact to return to for every future GBM enhancer question.

Output schema (parquet, long format):
    peak_id          str   — CATLAS peak id (chrN:start-end)
    chrom, start, end        — GRCh38 coordinates
    cell_type        str   — harmonized vocabulary
    cohort           str   — one of the catalog entries
    patient_id       str   — standardized across cohorts
    modality         str   — "snATAC" or "multiome"
    region           str   — if multi-region
    n_cells          int   — cells in this (cohort, patient, cell_type) bin
    n_accessible     int   — cells with ≥1 fragment at this peak
    frac_accessible  float — n_accessible / n_cells
    mean_counts      float — mean fragment count per cell
    rna_mean         float — mean target-gene expression (Multiome only; null otherwise)
    rna_coordinated  float — frac of accessible cells also expressing target gene (Multiome only)
"""

from __future__ import annotations
import argparse, sys
from pathlib import Path


def log(msg): print(f"[matrix] {msg}", flush=True)


SCHEMA_MD = """# GBM Enhancer Candidate Matrix — schema

A long-format parquet keyed on `(peak_id, cell_type, cohort, patient_id)`.
Query with polars or duckdb for fast slice-and-dice.

## Columns

| column           | type   | meaning |
|------------------|--------|---------|
| peak_id          | str    | CATLAS peak id, format `chrN:start-end` |
| chrom            | str    | GRCh38 chromosome |
| start            | int    | GRCh38 start |
| end              | int    | GRCh38 end |
| cell_type        | str    | harmonized vocabulary from config/celltypes.yaml |
| cohort           | str    | source dataset (guilhamon, tcga_scatac, mathewson_lupien, gbm_space, ...) |
| patient_id       | str    | standardized patient id |
| modality         | str    | "snATAC" or "multiome" |
| region           | str    | multi-region label (gbm_space); else "unknown" |
| n_cells          | int    | total cells in this (cohort, patient, cell_type) bin |
| n_accessible     | int    | cells with ≥1 fragment at this peak |
| frac_accessible  | float  | n_accessible / n_cells |
| mean_counts      | float  | mean fragment count per cell |
| rna_mean         | float  | mean paired-RNA expression of nearest gene (Multiome only; else NaN) |
| rna_coordinated  | float  | fraction of accessible cells also expressing nearest gene (Multiome only) |

## Example queries

```python
import polars as pl
m = pl.read_parquet("matrix/enhancer_candidate_matrix.parquet")

# Top TAM-specific enhancers that replicate in >=3 cohorts
tam = (
    m.filter(pl.col("cell_type") == "TAM")
     .group_by("peak_id")
     .agg([
         pl.col("frac_accessible").mean().alias("mean_frac"),
         pl.col("cohort").n_unique().alias("n_cohorts"),
         pl.col("patient_id").n_unique().alias("n_patients"),
     ])
     .filter(pl.col("n_cohorts") >= 3)
     .sort("mean_frac", descending=True)
     .head(20)
)
```
"""


def main() -> int:
    import anndata as ad
    import numpy as np
    import polars as pl
    import re as _re

    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, type=Path)
    ap.add_argument("--matrix", required=True, type=Path)
    ap.add_argument("--provenance", required=True, type=Path)
    ap.add_argument("--schema", required=True, type=Path)
    ap.add_argument("--threads", type=int, default=32)
    args = ap.parse_args()

    # Collect per-group polars DataFrames; concat once at the end.
    # Vectorized construction (no per-peak Python loop) is ~50-100x faster than
    # the dict-append version this replaced.
    all_dfs = []
    provenance_rows = []

    for src in args.inputs:
        log(f"=== Aggregating {src} ===")
        adata = ad.read_h5ad(src, backed="r")
        cohort = str(adata.obs["cohort"].iloc[0])
        modality = str(adata.obs["modality"].iloc[0])
        n_cells = adata.n_obs
        n_peaks = adata.n_vars
        log(f"  cohort={cohort}  n_cells={n_cells}  n_peaks={n_peaks}")

        provenance_rows.append(dict(
            cohort=cohort,
            n_cells=n_cells,
            n_patients=adata.obs["patient_id"].nunique(),
            n_samples=adata.obs["sample_id"].nunique(),
            modality=modality,
            cell_types=sorted(adata.obs["cell_type"].unique().tolist()),
            n_peaks=n_peaks,
            source_file=str(src),
        ))

        # PRECOMPUTE var-side arrays ONCE per cohort (not per group).
        # var_names like 'chr1:12345-67890' are either coming from a column trio
        # (chrom/start/end) or need to be parsed out of the index.
        has_chrom_col = "chrom" in adata.var.columns
        if has_chrom_col:
            chrom_all = adata.var["chrom"].astype(str).values
            start_all = adata.var["start"].astype(np.int64).values
            end_all = adata.var["end"].astype(np.int64).values
            peak_id_all = np.array(
                [f"{c}:{s}-{e}" for c, s, e in zip(chrom_all, start_all, end_all)],
                dtype=object,
            )
            peak_valid = np.ones(n_peaks, dtype=bool)
        else:
            var_names = np.asarray(adata.var_names, dtype=object)
            peak_id_all = var_names
            chrom_all = np.empty(n_peaks, dtype=object)
            start_all = np.zeros(n_peaks, dtype=np.int64)
            end_all = np.zeros(n_peaks, dtype=np.int64)
            peak_valid = np.zeros(n_peaks, dtype=bool)
            pat = _re.compile(r"^(chr[\dXYM]+):(\d+)-(\d+)$")
            for pi, name in enumerate(var_names):
                m = pat.match(str(name))
                if m:
                    chrom_all[pi] = m.group(1)
                    start_all[pi] = int(m.group(2))
                    end_all[pi] = int(m.group(3))
                    peak_valid[pi] = True
            if not peak_valid.all():
                log(f"  WARNING: {int((~peak_valid).sum())}/{n_peaks} vars had unparseable peak_id; dropping")

        groups = adata.obs.groupby(["patient_id", "cell_type", "region"], dropna=False, observed=True)
        n_groups = groups.ngroups
        log(f"  {n_groups} (patient, cell_type, region) groups to aggregate")

        for i, ((pid, ct, region), df) in enumerate(groups, 1):
            idx = df.index.values
            n_cells_group = len(idx)
            if n_cells_group < 10:
                if i % 50 == 0:
                    log(f"    group {i}/{n_groups} ({ct}, {pid}): {n_cells_group} cells — skipped (<10)")
                continue
            mat = adata[idx, :].X  # sparse CSR, shape (n_cells_group, n_peaks)
            n_accessible = np.asarray((mat > 0).sum(axis=0)).flatten()
            mean_counts = np.asarray(mat.mean(axis=0)).flatten()
            keep = (n_accessible > 0) & peak_valid
            k = int(keep.sum())
            if k == 0:
                if i % 50 == 0:
                    log(f"    group {i}/{n_groups} ({ct}, {pid}): {n_cells_group} cells, 0 peaks")
                continue

            # VECTORIZED: build a polars DataFrame directly from numpy slices
            df_group = pl.DataFrame({
                "peak_id":        peak_id_all[keep].astype(str),
                "chrom":          chrom_all[keep].astype(str),
                "start":          start_all[keep].astype(np.int64),
                "end":            end_all[keep].astype(np.int64),
                "cell_type":      [str(ct)] * k,
                "cohort":         [cohort] * k,
                "patient_id":     [str(pid)] * k,
                "region":         [str(region)] * k,
                "modality":       [modality] * k,
                "n_cells":        np.full(k, n_cells_group, dtype=np.int64),
                "n_accessible":   n_accessible[keep].astype(np.int64),
                "frac_accessible": (n_accessible[keep].astype(np.float64) / n_cells_group).astype(np.float32),
                "mean_counts":    mean_counts[keep].astype(np.float32),
                "rna_mean":       np.full(k, np.nan, dtype=np.float32),
                "rna_coordinated": np.full(k, np.nan, dtype=np.float32),
            })
            all_dfs.append(df_group)
            if i % 50 == 0 or i == n_groups:
                log(f"    group {i}/{n_groups} ({ct}, {pid}): {n_cells_group} cells, {k:,} peaks")

        # Explicitly close the backed file so we can free its file handle before the next
        try:
            adata.file.close()
        except Exception:
            pass

    log(f"Concatenating {len(all_dfs):,} per-group DataFrames...")
    big = pl.concat(all_dfs, how="vertical", rechunk=True)
    log(f"Final long-format matrix: {big.height:,} rows x {big.width} cols")

    args.matrix.parent.mkdir(parents=True, exist_ok=True)
    big.write_parquet(args.matrix, compression="zstd")
    log(f"Wrote matrix → {args.matrix}")

    pl.DataFrame(provenance_rows).write_parquet(args.provenance, compression="zstd")
    log(f"Wrote provenance → {args.provenance}")

    args.schema.write_text(SCHEMA_MD)
    log(f"Wrote schema doc → {args.schema}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
