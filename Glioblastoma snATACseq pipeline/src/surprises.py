"""
Surface patterns worth knowing that weren't the question being asked.
Each detector returns a dataframe of candidate "surprises" with provenance.
"""

from __future__ import annotations
import argparse, sys, yaml
from pathlib import Path


def log(msg): print(f"[surprises] {msg}", flush=True)


def detect_dual_celltype_enhancers(m, min_strength=0.3, max_other=0.05):
    """Peaks accessible in exactly two cell types with high selectivity against all others."""
    import polars as pl
    log("  Dual-cell-type enhancer detector")
    agg = (m.group_by(["peak_id", "cell_type"])
             .agg(pl.col("frac_accessible").mean().alias("mean_frac")))
    pivot = agg.pivot(index="peak_id", on="cell_type", values="mean_frac").fill_null(0)
    cell_types = [c for c in pivot.columns if c != "peak_id"]
    hits = []
    for row in pivot.iter_rows(named=True):
        high = [c for c in cell_types if row[c] >= min_strength]
        low  = [c for c in cell_types if row[c] <= max_other]
        if len(high) == 2 and len(low) >= (len(cell_types) - 3):
            hits.append({"peak_id": row["peak_id"], "ct1": high[0], "ct2": high[1],
                         "frac1": row[high[0]], "frac2": row[high[1]]})
    log(f"    {len(hits)} dual-cell-type enhancers")
    return pl.DataFrame(hits) if hits else pl.DataFrame(schema={
        "peak_id": pl.Utf8, "ct1": pl.Utf8, "ct2": pl.Utf8, "frac1": pl.Float64, "frac2": pl.Float64})


def detect_tumor_vs_normal_divergence(m, min_delta=0.3):
    """Peaks open in homeostatic microglia but closed in TAM at the same patient — regulatory shutdown."""
    import polars as pl
    log("  Tumor-vs-normal divergence detector (microglia vs TAM)")
    sub = m.filter(pl.col("cell_type").is_in(["microglia", "TAM"]))
    pivot = (sub.group_by(["peak_id", "patient_id", "cell_type"])
                .agg(pl.col("frac_accessible").mean().alias("mean_frac"))
                .pivot(index=["peak_id", "patient_id"], on="cell_type", values="mean_frac")
                .fill_null(0))
    if "microglia" not in pivot.columns or "TAM" not in pivot.columns:
        log("    insufficient coverage for microglia/TAM comparison")
        return pl.DataFrame(schema={"peak_id": pl.Utf8, "patient_id": pl.Utf8, "delta": pl.Float64, "direction": pl.Utf8})
    pivot = pivot.with_columns((pl.col("microglia") - pl.col("TAM")).alias("delta"))
    hits = pivot.filter(pl.col("delta").abs() >= min_delta)
    hits = hits.with_columns(
        pl.when(pl.col("delta") > 0).then(pl.lit("closed_in_TAM")).otherwise(pl.lit("opened_in_TAM")).alias("direction")
    )
    log(f"    {hits.height} per-patient divergence hits")
    return hits.select(["peak_id", "patient_id", "delta", "direction"])


def detect_patient_subpopulation_specific(m, high_thr=0.5, low_thr=0.05, min_split=0.4):
    """Peaks accessible in a sharp subset of patients — potential subtype markers (MGMT, Verhaak)."""
    import polars as pl
    log("  Patient-subpopulation-specific detector")
    agg = (m.group_by(["peak_id", "cell_type", "patient_id"])
             .agg(pl.col("frac_accessible").mean().alias("mean_frac")))
    hits = []
    for (peak, ct), df in agg.group_by(["peak_id", "cell_type"]):
        if df.height < 10:
            continue
        high = (df["mean_frac"] >= high_thr).sum()
        low  = (df["mean_frac"] <= low_thr).sum()
        frac_high = high / df.height
        if min_split <= frac_high <= (1 - min_split) and low >= 3 and high >= 3:
            hits.append({"peak_id": peak, "cell_type": ct, "n_high": int(high),
                         "n_low": int(low), "frac_high": float(frac_high)})
    log(f"    {len(hits)} subpopulation-specific hits")
    return pl.DataFrame(hits) if hits else pl.DataFrame(schema={
        "peak_id": pl.Utf8, "cell_type": pl.Utf8, "n_high": pl.Int64, "n_low": pl.Int64, "frac_high": pl.Float64})


def detect_multiome_discordance(m, atac_thr=0.3, rna_thr=0.1):
    """Multiome: peaks accessible without expression (poised) or expressed without accessibility (compensatory)."""
    import polars as pl
    log("  Multiome accessibility/expression discordance detector")
    multi = m.filter((pl.col("modality") == "multiome") & pl.col("rna_mean").is_not_null())
    if multi.height == 0:
        log("    no Multiome rows with RNA data; skipping")
        return pl.DataFrame(schema={"peak_id": pl.Utf8, "cell_type": pl.Utf8, "pattern": pl.Utf8})
    poised = multi.filter((pl.col("frac_accessible") >= atac_thr) & (pl.col("rna_mean") < rna_thr))
    compensatory = multi.filter((pl.col("frac_accessible") < atac_thr) & (pl.col("rna_mean") >= rna_thr))
    log(f"    {poised.height} poised (open, not expressed); {compensatory.height} compensatory")
    poised = poised.with_columns(pl.lit("poised").alias("pattern")).select(["peak_id", "cell_type", "pattern"])
    comp = compensatory.with_columns(pl.lit("compensatory").alias("pattern")).select(["peak_id", "cell_type", "pattern"])
    return pl.concat([poised, comp])


def main() -> int:
    import polars as pl

    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    log(f"Loading matrix: {args.src}")
    m = pl.read_parquet(args.src)

    results = {
        "dual_celltype": detect_dual_celltype_enhancers(m),
        "tumor_vs_normal": detect_tumor_vs_normal_divergence(m),
        "subpopulation": detect_patient_subpopulation_specific(m),
        "multiome_discordance": detect_multiome_discordance(m),
    }

    frames = []
    for name, df in results.items():
        if df.height > 0:
            frames.append(df.with_columns(pl.lit(name).alias("detector")))
    combined = pl.concat(frames, how="diagonal_relaxed") if frames else pl.DataFrame(schema={"detector": pl.Utf8})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.write_parquet(args.out, compression="zstd")
    log(f"Wrote → {args.out}  ({combined.height:,} surprise rows across {len(frames)} detectors)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
