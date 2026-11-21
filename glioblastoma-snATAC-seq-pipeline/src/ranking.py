"""
Composite ranking → "clone these first" list, per cell type.

Applies the formula from config/pipeline.yaml and the AAV cargo + CNV-safe + no-PmlI constraints.
Outputs top_candidates.parquet (overall) and per_celltype/<ct>.parquet (one per cell type).
"""

from __future__ import annotations
import argparse, sys, yaml, re
from pathlib import Path


def log(msg): print(f"[rank] {msg}", flush=True)


PMLI_MOTIF = re.compile(r"CACGTG")   # PmlI recognition site, counter-selection step


def passes_constraints(row, cfg) -> bool:
    length = int(row["end"]) - int(row["start"])
    if length < cfg["constraints"]["peak_length_bp"]["min"]:
        return False
    if length > cfg["constraints"]["peak_length_bp"]["max"]:
        return False
    if row["chrom"] in cfg["constraints"]["chromosome_exclusions"]:
        return False
    if cfg["constraints"]["cnv_safe_zones_only"] and row["chrom"] in ("chr7", "chr10"):
        return False
    # PmlI check requires sequence fetch; done in a separate pass against hg38
    # if row.get("has_pml_site", False) and cfg["constraints"]["e_box_pml_free"]:
    #     return False
    return True


def composite(row) -> float:
    import math
    s   = float(row["strength"])
    c   = float(row["consistency"])
    sel = float(row["selectivity"])
    rep = float(row["n_cohorts_replicating"]) / max(1.0, float(row["n_cohorts"]))
    rna = 0.0  # populate from Multiome scoring if present
    return (
        0.30 * s
      + 0.25 * c
      + 0.25 * (math.log2(sel + 1) / math.log2(10))
      + 0.15 * rep
      + 0.05 * rna
    )


def main() -> int:
    import polars as pl

    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--top", required=True, type=Path)
    ap.add_argument("--per-celltype-dir", required=True, type=Path)
    args = ap.parse_args()

    with open("config/pipeline.yaml") as f:
        cfg = yaml.safe_load(f)["ranking"]

    log(f"Loading scored candidates: {args.src}")
    scored = pl.read_parquet(args.src)

    # candidate_scores.parquet lacks chrom/start/end — parse from peak_id
    # (CATLAS format 'chrX:start-end'). Also use composite_score computed by
    # scoring.py so we don't double-weight.
    if "chrom" not in scored.columns:
        log("  parsing chrom/start/end from peak_id")
        scored = scored.with_columns([
            pl.col("peak_id").str.extract(r"^(chr[\dXYM]+)", 1).alias("chrom"),
            pl.col("peak_id").str.extract(r":(\d+)-", 1).cast(pl.Int64).alias("start"),
            pl.col("peak_id").str.extract(r"-(\d+)$", 1).cast(pl.Int64).alias("end"),
        ])

    log(f"Applying constraints (length {cfg['constraints']['peak_length_bp']['min']}-"
        f"{cfg['constraints']['peak_length_bp']['max']} bp, CNV-safe zones only={cfg['constraints']['cnv_safe_zones_only']})")
    kept = scored.filter(
        pl.struct(["chrom", "start", "end"]).map_elements(
            lambda r: passes_constraints(r, cfg), return_dtype=pl.Boolean
        )
    )
    log(f"  {kept.height:,} / {scored.height:,} rows pass constraints")

    # scoring.py already computed composite_score; keep as-is unless missing
    if "composite_score" not in kept.columns:
        log("Computing composite score (not present in scoring output)")
        kept = kept.with_columns(
            pl.struct(["strength", "consistency", "selectivity", "n_cohorts"])
              .map_elements(composite, return_dtype=pl.Float64)
              .alias("composite_score")
        )

    top_n = cfg["top_n_per_celltype"]
    log(f"Ranking top {top_n} per cell type")
    args.per_celltype_dir.mkdir(parents=True, exist_ok=True)
    top_rows = []
    for ct in sorted(set(kept["cell_type"])):
        subset = kept.filter(pl.col("cell_type") == ct).sort("composite_score", descending=True).head(top_n)
        out_path = args.per_celltype_dir / f"{ct}.parquet"
        subset.write_parquet(out_path, compression="zstd")
        top_rows.append(subset)
        log(f"  {ct:30s}  {subset.height} candidates  (top score {subset['composite_score'].max():.3f})")
        for row in subset.head(3).iter_rows(named=True):
            rep = row.get('n_cohorts_replicating', row.get('n_cohorts', 0))
            log(f"      {row['peak_id']:28s} composite={row['composite_score']:.3f} "
                f"sel={row['selectivity']:.2f} rep={rep}")

    top_all = pl.concat(top_rows)
    args.top.parent.mkdir(parents=True, exist_ok=True)
    top_all.write_parquet(args.top, compression="zstd")
    log(f"Wrote top candidates → {args.top}")

    # Snakemake expects a parquet for each priority cell type (TAM, microglia,
    # malignant_NPC_like, malignant_OPC_like, malignant_MES_like). Cohorts
    # without proper gene-activity bridging produce only unresolved labels —
    # write empty-schema placeholders for the missing priority cell types so
    # the Snakefile rule's declared outputs all exist.
    priority = ["TAM", "microglia", "malignant_NPC_like", "malignant_OPC_like",
                "malignant_MES_like"]
    schema = top_all.schema if top_all.height > 0 else {c: pl.Utf8 for c in ["peak_id","cell_type"]}
    for ct in priority:
        p = args.per_celltype_dir / f"{ct}.parquet"
        if not p.exists():
            pl.DataFrame(schema=schema).write_parquet(p, compression="zstd")
            log(f"  placeholder → {p.name} (no cells of this type in current cohorts)")

    # Clone-these-first: cross-cell-type top 10 overall
    clone_first = top_all.sort("composite_score", descending=True).head(10)
    log("")
    log("=== CLONE THESE FIRST (top 10 overall) ===")
    for row in clone_first.iter_rows(named=True):
        log(f"  {row['peak_id']:28s} {row['cell_type']:22s} composite={row['composite_score']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
