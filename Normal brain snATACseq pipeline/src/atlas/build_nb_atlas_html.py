#!/usr/bin/env python3
"""Build the Normal Brain atlas v3 HTML report — the NB analog of GBM atlas_v4.html.

Reads the existing replicability_major_v2.0_scored.parquet (long format, peak × major_class)
and builds a self-contained HTML with:
  - Per-major_class tabs (16 cell types), top-1000 distal peaks per tab
  - Columns matching GBM v4 one-for-one, with scoring semantics ported:
      log2FC_top_vs_next     = log2(specificity_ratio)
      spec_pct               = percentile rank of log2FC across distal peaks
      strength_mean_evaluable_wn = width-norm mean fraction_open across evaluable
                                    donors (donors without signal count as 0; same
                                    denominator as GBM atlas) = (sum_frac /
                                    n_evaluable_class) × (500 / peak_width)
      strength_pct           = percentile rank of mean_frac across distal peaks
      primary_score          = (spec_pct + strength_pct) / 2
      n_cohorts              = n_studies (0-23 for NB vs 0-8 for GBM)
      patient_detection      = n_donors / n_donors_evaluable_in_class
  - Expansion panel: CPM bars across 16 cell types + cohort breakdown
  - UCSC Genome Browser deep links
  - Caveats box rewritten for NB (up to 25 cohorts, no CNV amplicons, 4DN+GSE322785 sparse-cell-call caveat)

Usage:
    python3 build_nb_atlas_v3_html.py --parquet <path> --out <path>
"""
import argparse
import json
import math
from pathlib import Path

import polars as pl


DEFAULT_PARQUET = "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/replicability_major_v2.5_gbm_strength.parquet"
DEFAULT_OUT = "/Users/jpmcginnis1/mcginnislab-data/normal-brain-atlas/top-1000.html"
TOP_N_PER_CT = 1000
DISTAL_THRESHOLD = 2000  # bp from TSS; matches GBM
MIN_PEAK_WIDTH = 100  # drop sub-100bp peaks (width-norm amplification artifact)
EPSILON = 1e-9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", default=DEFAULT_PARQUET)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--top-n", type=int, default=TOP_N_PER_CT)
    args = ap.parse_args()

    print(f"Loading {args.parquet}", flush=True)
    df = pl.read_parquet(args.parquet)
    print(f"  shape: {df.shape}", flush=True)

    # Compute derived columns per (peak, class)
    print("Computing derived columns (log2FC, distance_band, distal flag)", flush=True)
    df = df.with_columns([
        # log2FC = log2(specificity_ratio), guarding against zero/negative
        pl.when(pl.col("specificity_ratio") > 0)
          .then(pl.col("specificity_ratio").log(base=2.0))
          .otherwise(-10.0)  # peaks with zero specificity get a floor
          .alias("log2FC_top_vs_next"),
        pl.col("distance_to_tss").abs().alias("abs_dist"),
    ]).with_columns([
        pl.when(pl.col("abs_dist") < DISTAL_THRESHOLD).then(pl.lit("proximal"))
          .when(pl.col("abs_dist") <= 10_000).then(pl.lit("2-10kb"))
          .otherwise(pl.lit(">10kb"))
          .alias("distance_band"),
        (pl.col("abs_dist") >= DISTAL_THRESHOLD).alias("is_distal"),
    ])

    # Filter to distal + peak-width floor (sub-100bp get amplified by width-norm artifact)
    distal = df.filter(pl.col("is_distal") & (pl.col("peak_width") >= MIN_PEAK_WIDTH))
    print(f"  distal + peak_width>={MIN_PEAK_WIDTH}: {distal.shape[0]:,} rows", flush=True)

    # Compute per-class percentile ranks using GBM-style strength (strength_wn_gbm)
    print("Computing per-class percentile ranks (spec_pct, strength_pct, primary_score)", flush=True)
    distal = distal.with_columns([
        pl.col("log2FC_top_vs_next").rank(method="average").over("major_class").alias("_spec_rank"),
        pl.col("strength_wn_gbm").rank(method="average").over("major_class").alias("_strength_rank"),
        pl.len().over("major_class").alias("_class_n"),
    ]).with_columns([
        (pl.col("_spec_rank") / pl.col("_class_n")).alias("spec_pct"),
        (pl.col("_strength_rank") / pl.col("_class_n")).alias("strength_pct"),
    ]).with_columns([
        ((pl.col("spec_pct") + pl.col("strength_pct")) / 2.0).alias("primary_score"),
    ]).drop(["_spec_rank", "_strength_rank", "_class_n"])

    # List of cell types — ONLY classes with real evaluable direct quantification
    # (at least one donor with >=100 cells of this class in the nnz-filtered pool).
    # Classes with n_evaluable_fresh == 0 across all peaks had their stale rows carried
    # over from pre-nnz-filter scoring and have no real data in the current atlas.
    evaluable_classes = (
        distal.group_by("major_class")
              .agg(pl.col("n_evaluable_fresh").max().alias("max_eval"))
              .filter(pl.col("max_eval") > 0)["major_class"]
              .to_list()
    )
    distal = distal.filter(pl.col("major_class").is_in(evaluable_classes))
    ct_counts = distal.group_by("major_class").len().sort("len", descending=True)
    ct_list = ct_counts["major_class"].to_list()
    print(f"  {len(ct_list)} evaluable major_class values: {ct_list}", flush=True)
    dropped = sorted(set(distal.select(pl.col('major_class')).unique().to_series().to_list()) ^ set(ct_list))
    if dropped:
        print(f"  dropped (no evaluable donors): {dropped}", flush=True)

    # Top-N per cell type
    print(f"Selecting top-{args.top_n} per cell type by primary_score", flush=True)
    top_per_ct = {}
    for ct in ct_list:
        top = (
            distal.filter(pl.col("major_class") == ct)
                  .sort("primary_score", descending=True)
                  .head(args.top_n)
        )
        top_per_ct[ct] = top
        print(f"    {ct}: top {top.shape[0]:,} selected", flush=True)

    # Surgical-only alt ranking (GBM-style):
    # Only classes with surgical donors (n_surgical_evaluable_class > 0) get a surgical
    # top-1000. strength = surg_strength_wn_gbm (width-norm mean penetrance over the ~19
    # surgical donors; donors with no signal count as 0). spec = combined-donor log2FC
    # (no surgical-only off-class signal in current parquet; flagged in UI).
    print(f"Selecting surgical-only top-{args.top_n} per cell type (GBM-style)", flush=True)
    distal_surg = (
        distal.filter(pl.col("n_surgical_evaluable_class") > 0)
              .with_columns([
                  pl.col("log2FC_top_vs_next").rank(method="average").over("major_class").alias("_surg_spec_rank"),
                  pl.col("surg_strength_wn_gbm").rank(method="average").over("major_class").alias("_surg_strength_rank"),
                  pl.len().over("major_class").alias("_surg_class_n"),
              ])
              .with_columns([
                  (pl.col("_surg_spec_rank") / pl.col("_surg_class_n")).alias("surg_spec_pct"),
                  (pl.col("_surg_strength_rank") / pl.col("_surg_class_n")).alias("surg_strength_pct"),
              ])
              .with_columns([
                  ((pl.col("surg_spec_pct") + pl.col("surg_strength_pct")) / 2.0).alias("surgical_primary_score"),
              ])
              .drop(["_surg_spec_rank", "_surg_strength_rank", "_surg_class_n"])
    )
    top_per_ct_surgical = {}
    for ct in ct_list:
        sub = distal_surg.filter(pl.col("major_class") == ct)
        if sub.shape[0] == 0:
            top_per_ct_surgical[ct] = sub
            print(f"    {ct}: 0 (no surgical donors evaluable)", flush=True)
            continue
        top = sub.sort("surgical_primary_score", descending=True).head(args.top_n)
        top_per_ct_surgical[ct] = top
        print(f"    {ct}: top {top.shape[0]:,} selected (surgical)", flush=True)

    # Build TAB_DATA (per-ct rows) — GBM-style strength/detection
    print("Building TAB_DATA (default, GBM-style)", flush=True)
    tab_data = {}
    all_top_peak_ids = set()
    for ct, top in top_per_ct.items():
        rows = []
        for row in top.iter_rows(named=True):
            rows.append({
                "peak_id": row["peak_id"],
                "chrom": row["chrom"],
                "start": int(row["start"]),
                "end": int(row["end"]),
                "peak_width_bp": int(row["peak_width"]),
                "nearest_gene": row.get("nearest_gene") or "",
                "dist_to_tss_signed": int(row.get("distance_to_tss") or 0),
                "distance_band": row["distance_band"],
                "log2FC_top_vs_next": float(row["log2FC_top_vs_next"]),
                "spec_pct": float(row["spec_pct"]),
                "strength_mean_evaluable_wn": float(row["strength_wn_gbm"]),
                "strength_pct": float(row["strength_pct"]),
                "primary_score": float(row["primary_score"]),
                "n_cohorts": int(row["n_studies"]),
                "n_patients_detected": int(row["n_detected"]),
                "n_patients_evaluable": int(row["n_evaluable_class"]),
                "patient_detection_rate": float(row["patient_detection_rate_gbm"]),
            })
            all_top_peak_ids.add(row["peak_id"])
        tab_data[ct] = rows

    # Build TAB_DATA_SURGICAL — GBM-style surgical columns
    print("Building TAB_DATA_SURGICAL (GBM-style)", flush=True)
    tab_data_surgical = {}
    for ct, top in top_per_ct_surgical.items():
        rows = []
        if top.shape[0] > 0:
            for row in top.iter_rows(named=True):
                n_eval = int(row["n_surgical_evaluable_class"] or 0)
                n_det = int(row["n_detected_surg"] or 0)
                rate = (n_det / n_eval) if n_eval > 0 else 0.0
                rows.append({
                    "peak_id": row["peak_id"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "peak_width_bp": int(row["peak_width"]),
                    "nearest_gene": row.get("nearest_gene") or "",
                    "dist_to_tss_signed": int(row.get("distance_to_tss") or 0),
                    "distance_band": row["distance_band"],
                    "log2FC_top_vs_next": float(row["log2FC_top_vs_next"]),
                    "spec_pct": float(row["surg_spec_pct"]),
                    "strength_mean_evaluable_wn": float(row["surg_strength_wn_gbm"] or 0.0),
                    "strength_pct": float(row["surg_strength_pct"]),
                    "primary_score": float(row["surgical_primary_score"]),
                    "n_cohorts": int(row["n_studies"]),
                    "n_patients_detected": n_det,
                    "n_patients_evaluable": n_eval,
                    "patient_detection_rate": rate,
                })
                all_top_peak_ids.add(row["peak_id"])
        tab_data_surgical[ct] = rows

    # Build EXPANSION: for each top peak, pull strength_wn_gbm across all 16 classes
    print(f"Building EXPANSION data for {len(all_top_peak_ids):,} unique peaks", flush=True)
    pivot = (
        df.filter(pl.col("peak_id").is_in(list(all_top_peak_ids)))
          .select(["peak_id", "major_class", "strength_wn_gbm", "n_donors", "n_studies"])
    )
    # One row per peak_id with strength per class in CT_LIST order
    peak_class_cpms = {}
    for peak_id, group in pivot.group_by("peak_id"):
        cpm_map = {r["major_class"]: float(r["strength_wn_gbm"]) for r in group.iter_rows(named=True)}
        cpms = [cpm_map.get(ct, 0.0) for ct in ct_list]
        peak_class_cpms[peak_id[0]] = cpms

    # Cohort-level detection: n_studies tells us total cohorts; we'd need per-cohort breakdown
    # from direct_parquets to show chips. For MVP, show total count summary.
    expansion = {}
    for ct in ct_list:
        exp_for_ct = {}
        # include both default and surgical top peaks for this cell type
        pid_seen = set()
        for src in (top_per_ct.get(ct), top_per_ct_surgical.get(ct)):
            if src is None or src.shape[0] == 0:
                continue
            for row in src.iter_rows(named=True):
                pid = row["peak_id"]
                if pid in pid_seen:
                    continue
                pid_seen.add(pid)
                cpms = peak_class_cpms.get(pid, [0.0] * len(ct_list))
                exp_for_ct[pid] = {
                    "cpms": cpms,
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "n_studies_this_class": int(row["n_studies"]),
                    "n_donors_this_class": int(row["n_donors"]),
                    "surgical_pos": int(row.get("surgical_pos") or 0),
                    "surgical_total": int(row.get("surgical_total") or 0),
                    "postmortem_pos": int(row.get("postmortem_pos") or 0),
                    "postmortem_total": int(row.get("postmortem_total") or 0),
                    "specificity_delta": float(row.get("specificity_delta") or 0),
                    "specificity_ratio": float(row.get("specificity_ratio") or 0),
                }
        expansion[ct] = exp_for_ct

    # ATLAS_MAX_PER_CT — max strength_wn_gbm per class across the atlas. Each class is divided
    # by its OWN ceiling. Read vertically within a class ("how notable is this peak for this
    # class vs the strongest peak that class ever reaches"), NOT horizontally across classes
    # (small-ceiling classes like Dopaminergic will show large bars at modest absolute signal).
    # Useful for predicting AAV off-target leak: if a non-target class approaches its own
    # ceiling at this peak, that class is likely to drive the enhancer at a sizeable fraction
    # of its typical enhancer activity.
    atlas_max = {}
    for row in df.group_by("major_class").agg(pl.col("strength_wn_gbm").max()).iter_rows(named=True):
        atlas_max[row["major_class"]] = float(row["strength_wn_gbm"] or 0.0) or 1.0

    # Compute narrow Design B sub-peaks for each unique top peak (AAV-cloneable fragments
    # inside the wider v1.0.0 peak). Lookup table lives at ~/.cache/nb_atlas/narrow_union_B.bed
    # (47 MB, 2.06M peaks, median width 518 bp).
    # Union over default + surgical top-1000 sets so both modes have sub-peaks.
    tab_data_both = {ct: tab_data.get(ct, []) + tab_data_surgical.get(ct, []) for ct in ct_list}
    narrow_subpeaks = compute_narrow_subpeaks(tab_data_both, ct_list)

    # Classify tabs into primary (common) vs secondary (less common) based on cell counts
    primary_cts = ct_list[:10]
    secondary_cts = ct_list[10:]

    print(f"Writing HTML to {args.out}", flush=True)
    html = render_html(tab_data, tab_data_surgical, expansion, ct_list, atlas_max, narrow_subpeaks,
                       primary_cts, secondary_cts, top_n=args.top_n)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html)
    size_mb = out_path.stat().st_size / 1_048_576
    print(f"  wrote {out_path} ({size_mb:.1f} MB)", flush=True)


def compute_narrow_subpeaks(tab_data, ct_list):
    """For every unique peak in TAB_DATA, find overlapping narrow Design B peaks.

    Returns {peak_id: [{s, e, w}]} using short keys to minimize serialized JSON size
    (16k peaks × ~2-5 sub-peaks each would be ~1 MB with short keys).
    """
    import os
    from collections import defaultdict
    from pathlib import Path as _P

    bed_cache = _P.home() / ".cache" / "nb_atlas" / "narrow_union_B.bed"
    env_bed = os.environ.get("NARROW_BED")
    bed_path = _P(env_bed) if env_bed else bed_cache
    if not bed_path.exists():
        print(f"[narrow] WARN: {bed_path} not found; sub-peaks will be empty", flush=True)
        return {}
    print(f"[narrow] loading {bed_path}", flush=True)
    by_chrom = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            toks = line.rstrip().split("\t")
            if len(toks) < 3:
                continue
            try:
                by_chrom[toks[0]].append((int(toks[1]), int(toks[2])))
            except ValueError:
                continue
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    total = sum(len(v) for v in by_chrom.values())
    print(f"[narrow]   {total:,} narrow B peaks loaded across {len(by_chrom)} chroms", flush=True)

    unique = {}  # peak_id -> (chrom, start, end)
    for ct in ct_list:
        for row in tab_data.get(ct, []):
            pid = row["peak_id"]
            if pid not in unique:
                unique[pid] = (row["chrom"], int(row["start"]), int(row["end"]))

    out = {}
    missing_chrom = 0
    for pid, (chrom, start, end) in unique.items():
        if chrom not in by_chrom:
            missing_chrom += 1
            continue
        hits = []
        for b_start, b_end in by_chrom[chrom]:
            if b_start >= end:
                break
            if b_end > start:
                hits.append({"s": b_start, "e": b_end, "w": b_end - b_start})
        if hits:
            out[pid] = hits
    print(f"[narrow]   mapped {len(out):,} / {len(unique):,} peaks "
          f"({sum(len(v) for v in out.values()):,} total sub-peaks)", flush=True)
    if missing_chrom:
        print(f"[narrow]   {missing_chrom:,} peaks skipped (chrom not in BED)", flush=True)
    return out


def render_html(tab_data, tab_data_surgical, expansion, ct_list, atlas_max, narrow_subpeaks,
                primary_cts, secondary_cts, top_n=1000):
    # Serialize data as compact JSON
    tab_data_json = json.dumps(tab_data, separators=(",", ":"))
    tab_data_surgical_json = json.dumps(tab_data_surgical, separators=(",", ":"))
    expansion_json = json.dumps(expansion, separators=(",", ":"))
    ct_list_json = json.dumps(ct_list)
    atlas_max_json = json.dumps(atlas_max)
    narrow_json = json.dumps(narrow_subpeaks, separators=(",", ":"))

    # Build tab buttons
    tabs_html = []
    for ct in primary_cts:
        n = len(tab_data.get(ct, []))
        tabs_html.append(f'<button class="tab-btn primary" data-ct="{ct}">{ct} <span class="n">({n})</span></button>')
    for ct in secondary_cts:
        n = len(tab_data.get(ct, []))
        tabs_html.append(f'<button class="tab-btn secondary" data-ct="{ct}">{ct} <span class="n">({n})</span></button>')
    tabs_block = "\n".join(tabs_html)

    # Default tab = first primary cell type with data
    default_tab = next((ct for ct in primary_cts if tab_data.get(ct)), ct_list[0])

    return TEMPLATE.format(
        tabs_block=tabs_block,
        default_tab=default_tab,
        top_n=top_n,
        tab_data=tab_data_json,
        tab_data_surgical=tab_data_surgical_json,
        expansion=expansion_json,
        ct_list=ct_list_json,
        atlas_max=atlas_max_json,
        narrow_subpeaks=narrow_json,
    )


TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>Normal Brain snATAC Atlas v3 — Enhancer Candidates</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css">
<style>
  *,*::before,*::after {{ box-sizing: border-box; }}
  :root {{ --sticky-head-top: 173px; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
    margin: 0; padding: 0; color: #222; background: #f8f9fa; font-size: 13px; }}
  header {{ background:#1b1f24; color:#f8f9fa; padding:14px 24px;
    position:sticky; top:0; z-index:50; box-shadow:0 2px 4px rgba(0,0,0,0.15); }}
  header h1 {{ margin:0 0 4px 0; font-size:18px; font-weight:600; }}
  header .sub {{ font-size:11px; opacity:0.7; }}
  header .version-badge {{ background:#8b5cf6; color:white; font-weight:700;
    font-size:10px; padding:2px 8px; border-radius:4px; margin-left:10px; vertical-align:middle; }}
  .tabs {{ background:white; padding:10px 24px; border-bottom:1px solid #dee2e6;
    display:flex; flex-wrap:wrap; gap:6px; position:sticky; top:var(--header-height, 56px); z-index:40; }}
  .tab-btn {{ padding:6px 12px; font-size:12px; border:1px solid #ced4da; border-radius:4px;
    background:white; cursor:pointer; color:#495057; font-family:inherit; }}
  .tab-btn:hover {{ background:#f1f3f5; }}
  .tab-btn.active {{ background:#8b5cf6; color:white; border-color:#8b5cf6; }}
  .tab-btn.secondary {{ border-style:dashed; opacity:0.85; }}
  .tab-btn .n {{ font-weight:400; font-size:10px; opacity:0.7; }}
  .col-ctrl {{ background:white; padding:8px 24px; border-bottom:1px solid #dee2e6;
    display:flex; gap:12px; align-items:center; font-size:11px; color:#666; flex-wrap:wrap; }}
  .col-ctrl label {{ cursor:pointer; }}
  .col-ctrl .group-label {{ font-weight:600; color:#333; margin-right:4px; }}
  .col-ctrl select {{ padding:3px 8px; font-size:11px; border:1px solid #ced4da; border-radius:3px; background:white; }}
  body.surgical-mode .tab-btn.active {{ background:#047857; border-color:#047857; }}
  body.surgical-mode table.dataTable thead th {{ border-bottom-color:#10b981 !important; }}
  body.surgical-mode .primary-score-cell {{ background:#dcfce7 !important; }}
  body.surgical-mode #surgical-banner {{ display:flex; }}
  .tab-btn.empty-ct {{ opacity:0.35; text-decoration:line-through; cursor:not-allowed; }}
  #surgical-banner {{ display:none; background:#dcfce7; border:1px solid #10b981;
    border-radius:4px; padding:8px 12px; margin:12px 0 4px 0; font-size:11px; color:#064e3b;
    gap:8px; align-items:flex-start; }}
  #surgical-banner strong {{ color:#065f46; }}
  .main {{ padding:16px 24px; }}
  .caveats {{ background:#fff3cd; border:1px solid #ffe08a; border-radius:4px;
    padding:10px 16px; margin-bottom:16px; font-size:12px; color:#664d03; line-height:1.5; }}
  .caveats strong {{ color:#594100; }}
  table.dataTable {{ border-collapse: separate !important; border-spacing: 0 !important; font-size:11px; }}
  table.dataTable thead th {{ background:#343a40 !important; color:#f8f9fa !important; padding:7px 8px !important;
    font-weight:600; font-size:10px; text-transform:uppercase; letter-spacing:0.3px;
    cursor:help; border-bottom:2px solid #8b5cf6; white-space:nowrap;
    position: sticky !important; top: var(--sticky-head-top) !important; z-index: 45 !important;
    box-shadow: 0 2px 4px rgba(0,0,0,0.25); }}
  table.dataTable thead th[title]:hover {{ background:#495057 !important; }}
  table.dataTable tbody td {{ padding:5px 8px !important; white-space:nowrap; }}
  table.dataTable tbody tr:nth-child(even) {{ background:#fafbfc; }}
  table.dataTable tbody tr:hover {{ background:#f3edff; }}
  .num {{ font-family:ui-monospace, Menlo, monospace; text-align:right; }}
  .peakid {{ font-family:ui-monospace, Menlo, monospace; font-size:10px; }}
  .sparse-flag {{ background:#fef3c7; color:#78350f; padding:1px 6px; border-radius:3px; font-size:10px; }}
  .primary-score-cell {{ font-weight:700; background:#f3edff; }}
  .expander {{ cursor:pointer; color:#8b5cf6; user-select:none; font-weight:600; }}
  .no-expand {{ padding:12px 16px; color:#999; font-size:11px; background:#f8f9fa; }}
  .expansion {{ background:#f8f9fa; padding:12px 16px;
    max-width: calc(100vw - 48px); box-sizing: border-box; overflow: hidden; }}
  .expansion .ucsc-link, .expansion .ucsc-ctx, .expansion .aav-note {{ white-space: normal; }}
  .expansion h4 {{ margin:0 0 6px 0; font-size:11px; text-transform:uppercase; color:#495057;
    letter-spacing:0.3px; }}
  .expansion .ucsc-link {{ padding:10px 12px; background:white; border:1px solid #8b5cf6;
    border-radius:4px; margin-bottom:12px; font-size:11px; }}
  .expansion .ucsc-link a {{ color:#8b5cf6; font-weight:600; text-decoration:none; }}
  .expansion .ucsc-link a:hover {{ text-decoration:underline; }}
  .expansion .ucsc-ctx {{ color:#666; font-size:10px; }}
  .expansion .aav-note {{ display:inline-block; margin-top:4px; font-size:11px; color:#495057; }}
  .expansion .ct-bar-pair {{ display:grid; grid-template-columns: 1fr 1fr; gap:24px; margin-bottom:12px; }}
  .expansion .ct-bars {{ display:grid; grid-template-columns:140px 1fr 60px; gap:4px 10px;
    align-items:center; font-size:10px; }}
  .expansion .ct-bars .bar-cont {{ background:#e9ecef; height:12px; border-radius:2px; position:relative; }}
  .expansion .ct-bars .bar-fill {{ background:#8b5cf6; height:100%; border-radius:2px; }}
  .expansion .ct-bars .bar-fill.abs-scale {{ background:#f59e0b; }}
  .expansion .ct-bars .ct-name.highlighted {{ font-weight:700; color:#8b5cf6; }}
  .expansion .meta-grid {{ display:grid; grid-template-columns:repeat(auto-fill, minmax(180px, 1fr));
    gap:8px; font-size:11px; margin-bottom:12px; }}
  .expansion .meta-card {{ background:white; padding:8px 12px; border-radius:4px; border:1px solid #e9ecef; }}
  .expansion .meta-card .lab {{ font-size:10px; color:#666; text-transform:uppercase; letter-spacing:0.3px; }}
  .expansion .meta-card .val {{ font-size:14px; font-weight:600; color:#222; }}
  .expansion .narrow-box {{ background:white; border:1px solid #10b981; border-radius:4px;
    padding:10px 12px; margin-bottom:12px; font-size:11px; }}
  .expansion .narrow-box h4 {{ color:#047857; margin-bottom:6px; }}
  .expansion .narrow-box.narrow-empty {{ border-color:#dc143c; }}
  .expansion .narrow-box.narrow-empty h4 {{ color:#991b1b; }}
  .expansion .narrow-ctx {{ color:#666; font-size:10px; margin-bottom:8px; }}
  .expansion .narrow-cli {{ color:#666; font-size:10px; margin-top:8px; padding-top:6px;
    border-top:1px dashed #e9ecef; }}
  .expansion .narrow-cli code {{ background:#f3f4f6; padding:1px 5px; border-radius:3px;
    font-family:ui-monospace, Menlo, monospace; color:#047857; }}
  .expansion .sub-list {{ display:flex; flex-direction:column; gap:4px; }}
  .expansion .sub-row {{ display:grid;
    grid-template-columns: minmax(190px, auto) 70px 100px 60px 70px;
    gap:12px; align-items:center; padding:4px 8px; background:#f0fdf4;
    border-radius:3px; font-size:11px; }}
  .expansion .sub-coord {{ font-family:ui-monospace, Menlo, monospace; color:#064e3b; font-weight:600; }}
  .expansion .sub-width {{ font-family:ui-monospace, Menlo, monospace; color:#333; text-align:right; }}
  .expansion .sub-row a {{ color:#8b5cf6; text-decoration:none; font-size:10px; }}
  .expansion .sub-row a:hover {{ text-decoration:underline; }}
  .expansion .aav-ok {{ color:#047857; font-weight:600; font-size:10px; }}
  .expansion .aav-mid {{ color:#b45309; font-weight:600; font-size:10px; }}
  .expansion .aav-no {{ color:#b91c1c; font-weight:600; font-size:10px; }}
  .footer {{ font-size:10px; color:#666; padding:20px 24px; text-align:center; }}
</style>
</head><body>
<header>
  <h1>Normal Brain snATAC Atlas v3 — Enhancer Candidates <span class="version-badge">v3 · NB · rank-percentile sum</span></h1>
  <div class="sub">
    Up to 830 donors × 25 cohorts in abundant cell types (OPC 830 / Oligo 745 / Inhibitory 704 / Astrocyte 704 / Microglia 661 / Excitatory 528 / Endothelial 362); rare types (Meninges, T_cell, Serotonergic, Dopaminergic) have single-digit to ~20 donors. 41M barcodes (6M post-nnz≥1k QC) · 990K union peaks (v1.0.0 peak set) · 16 major cell-type classes.
    Top-{top_n} DISTAL peaks (|dist to TSS| ≥ 2 kb) per cell type,
    ranked by <b>primary_score = (spec_pct + strength_pct) / 2</b>.
    Specificity = log2(mean_frac target / max off-target). Strength = mean per-donor fraction of cells with ≥1 fragment in peak.
    Click row to expand per-cell-type CPM breakdown + per-sample-type (surgical/postmortem) split.
  </div>
</header>

<div class="tabs">{tabs_block}
</div>

<div class="col-ctrl">
  <span class="group-label">Distance band:</span>
  <select id="dist-filter">
    <option value="all">All distal (≥ 2 kb)</option>
    <option value="2-10kb">2–10 kb only</option>
    <option value=">10kb">&gt; 10 kb only</option>
  </select>
  <label style="margin-left:12px;" title="Switch to a surgical-only top-{top_n} per cell type, re-ranked using surgical-donor replication. See caveat box below for the exact definition.">
    <input type="checkbox" id="surgical-only"> <strong>Surgical only</strong>
    <span style="color:#666; font-size:10px;">(re-ranked top-{top_n})</span>
  </label>
  <span class="group-label" style="margin-left:12px;">Toggle columns:</span>
  <label><input type="checkbox" class="col-toggle" data-col="1"> chrom</label>
  <label><input type="checkbox" class="col-toggle" data-col="2"> start</label>
  <label><input type="checkbox" class="col-toggle" data-col="3"> end</label>
</div>

<div class="main">
  <div id="surgical-banner">
    <span>🩺</span>
    <div>
      <strong>Surgical only — top-{top_n} re-ranked (GBM-style).</strong> Only shows cell types with surgical donors (classes with ≥1 surgical donor meeting the ≥100-cells evaluability floor). <b>Strength</b> = <code>surg_strength_wn_gbm</code> = width-normalized mean penetrance across ALL surgical evaluable donors (donors without signal count as 0, same convention as GBM). <b>Patient detection</b> = n_detected / n_surgical_evaluable; denominator is CONSTANT per class (19 for most abundant types, 18 Excitatory, 14 Endothelial). <b>Specificity</b> still uses combined-donor <code>log2FC_top_vs_next</code> — no surgical-only off-class signal in the current parquet. Rare cell types with no surgical donors (Meninges, T_cell, Serotonergic, Dopaminergic, Immature_neuron, Neuron_other, Immune_other, Progenitor) are grayed-out tabs.
    </div>
  </div>
  <div class="caveats">
    <strong>v3 primary ranking (rank-percentile sum):</strong> Top-{top_n} per cell type ranked by
    <code>primary_score = (spec_pct + strength_pct) / 2</code>, where <code>spec_pct</code> is the percentile rank of
    <code>log2FC_top_vs_next = log2(mean_frac_target / off_target_max)</code> and <code>strength_pct</code> is the percentile rank of
    <code>mean_frac</code>, both computed within distal peaks of each cell type.
    <br>
    <strong>v3 strength (GBM-matched, width-normalized):</strong>
    <code>strength_wn_gbm = mean over ALL EVALUABLE donors (≥100 cells of class) of fraction_open × (500 / peak_width_bp)</code>, where fraction_open = cells accessing peak / total cells of that class per donor. Donors with no signal at a peak count as 0, not as unevaluable — denominator is CONSTANT per class (Excitatory/Inhibitory 225, Astrocyte/Oligo 212, Microglia 205, OPC 186, Mural 54, Endothelial 44, Progenitor 32, Dopaminergic 11, Serotonergic 2). Width-norm to 500 bp reference matches GBM v4 convention.
    <br><br>
    <strong>Distal filter:</strong> |dist_to_tss| ≥ 2 kb applied before top-N selection. Mirrors GBM v4 convention — promoter-proximal peaks aren't useful AAV enhancer scaffolds and would otherwise dominate the top ranks.
    <br><br>
    <strong>Peak set is v1.0.0 (pre-Design-B).</strong> This report is scored against the 990K pre-Design-B union peak set, median width 895 bp (mean 2,434 bp), which carries wider peaks than the current Design-B narrow union (2.06M peaks, median 518 bp). Chrom/start/end coordinates are correct; coord lookups in UCSC work. For AAV cloning, use the sister script <code>narrow_lookup.py</code> to resolve a top-ranked wide peak to the narrow Design B sub-peak(s) with the actual cloneable sequence.
    <br><br>
    <strong>Sparse-cell-call caveat.</strong> The 4DN hippocampus cohort (25.7M barcodes) and GEO GSE322785 (9.7M barcodes) shipped with permissive upstream barcode-calling; median nnz per barcode is 3-4 in those cohorts vs 5,000-15,000 in the rest. The upstream v2.0 pipeline inherited whatever filtering was applied by the author deposits. If a peak's apparent cell-type signal is dominated by patient counts from 4DN or GSE322785 and no other cohort replicates, treat the specificity call with caution; cross-check <code>n_cohorts</code> (should be ≥ 3) and surgical/postmortem agreement.
    <br><br>
    <strong>Nearest-gene is informational only.</strong> AAV scaffold efficacy depends on cell-type-selective driving strength (what primary_score ranks), not nearest-gene identity. Median mammalian enhancer-promoter distance is ~50-100 kb with a long tail beyond 1 Mb.
    <br><br>
    <strong>Hover any column header for its definition and formula.</strong>
  </div>

  <div id="table-container">
    <table id="atlas-table" class="display compact" style="width:100%">
      <thead>
        <tr>
          <th title="GRCh38 peak coordinates (chr:start-end).">peak_id</th>
          <th title="Chromosome, GRCh38. Hidden by default.">chrom</th>
          <th title="Peak start coordinate (0-based). Hidden by default.">start</th>
          <th title="Peak end coordinate (0-based). Hidden by default.">end</th>
          <th title="Peak width in bp (end − start). Note peaks on v1.0.0 peak set are wider than the Design B narrow union; use narrow_lookup.py to resolve to AAV-cloneable width.">width (bp)</th>
          <th title="HGNC symbol of nearest annotated TSS. CONTEXT ONLY — does NOT identify the regulatory target.">nearest gene</th>
          <th title="Signed distance from peak midpoint to the nearest TSS. Negative = upstream on gene's sense strand.">dist to TSS</th>
          <th title="Distance band: '2-10kb' or '>10kb'. Use the Distance band dropdown to filter.">band</th>
          <th title="Specificity: log2 fold-change of mean_frac in the target cell type vs the max among off-target cell types (across all 16 major classes).">log2FC top vs next</th>
          <th title="Percentile rank of log2FC_top_vs_next across all distal peaks within this cell type. 0 = lowest, 1 = highest specificity.">spec pct</th>
          <th title="STRENGTH (strength_mean_evaluable_wn, same formula as GBM atlas): mean over all evaluable donors (donors with ≥100 cells of this class) of the fraction of cells of this class with ≥1 Tn5 fragment in the peak, multiplied by (500 / peak_width_bp) to normalize for peak width.">strength (mean frac)</th>
          <th title="Percentile rank of strength across all distal peaks within this cell type. 0 = lowest, 1 = highest.">strength pct</th>
          <th title="DEFAULT RANKING. (spec_pct + strength_pct) / 2. Equal-weighted sum of nonparametric percentile ranks of specificity and strength. 0-1, higher is better.">primary score</th>
          <th title="Count of distinct studies (0-23) in which this peak × cell type is detected. Measures cross-cohort replication across the full NB atlas.">n cohorts</th>
          <th title="Detected in N of M evaluable donors. Evaluable = donor contributed enough cells of this class to be scored. Detected = donor's fraction_open crossed the pipeline's internal threshold.">donor detection</th>
        </tr>
      </thead>
      <tbody></tbody>
    </table>
  </div>
</div>

<div class="footer">
  Normal Brain atlas v3 · scored on v1.0.0 peak set · 16 major cell classes · up to 25 cohorts · up to 830 donors (abundant types) · {top_n} per tab
</div>

<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>
<script>
const TAB_DATA_DEFAULT = {tab_data};
const TAB_DATA_SURGICAL = {tab_data_surgical};
let TAB_DATA = TAB_DATA_DEFAULT;
const EXPANSION = {expansion};
const CT_LIST = {ct_list};
const ATLAS_MAX_PER_CT = {atlas_max};
const NARROW_SUBPEAKS = {narrow_subpeaks};

function parsePeak(p) {{
  const m = p.match(/^(chr[^:]+):(\d+)-(\d+)$/);
  return m ? {{ chrom: m[1], start: parseInt(m[2]), end: parseInt(m[3]) }} : null;
}}
function fmtDist(d) {{
  const sign = d < 0 ? '-' : '+';
  const abs = Math.abs(d);
  if (abs < 1000) return `${{sign}}${{abs}} bp`;
  if (abs < 1e6)  return `${{sign}}${{(abs/1000).toFixed(1)}} kb`;
  return `${{sign}}${{(abs/1e6).toFixed(2)}} Mb`;
}}
function fmtDetection(r) {{
  const pct = (100 * r.patient_detection_rate).toFixed(0);
  return `${{pct}}% (${{r.n_patients_detected}}/${{r.n_patients_evaluable}})`;
}}

let dt, currentCT = null, currentDistFilter = 'all';

function renderExpansion(peak_id, tab_ct) {{
  const data = EXPANSION[tab_ct] && EXPANSION[tab_ct][peak_id];
  if (!data) {{
    return `<div class="no-expand">(expansion data available for top rows per tab only)</div>`;
  }}
  const max_cpm = Math.max(...data.cpms, 1e-6);
  let ctBars = '';
  let ctBarsAbs = '';
  for (let i = 0; i < CT_LIST.length; i++) {{
    const ct = CT_LIST[i];
    const cpm = data.cpms[i];
    const w = (100 * cpm / max_cpm).toFixed(1);
    const atlasMax = ATLAS_MAX_PER_CT[ct] || 1;
    const wAbs = Math.min(100, (100 * cpm / atlasMax)).toFixed(1);
    const hl = (ct === tab_ct) ? 'highlighted' : '';
    ctBars += `<div class="ct-name ${{hl}}">${{ct}}</div>`
           + `<div class="bar-cont"><div class="bar-fill" style="width:${{w}}%"></div></div>`
           + `<div class="num">${{cpm.toFixed(4)}}</div>`;
    ctBarsAbs += `<div class="ct-name ${{hl}}">${{ct}}</div>`
              + `<div class="bar-cont"><div class="bar-fill abs-scale" style="width:${{wAbs}}%" title="${{cpm.toFixed(4)}} / ${{atlasMax.toFixed(4)}} atlas-max"></div></div>`
              + `<div class="num">${{wAbs}}%</div>`;
  }}
  // Meta cards: cohorts, donors, surgical/postmortem split, delta/ratio
  const surgRate = data.surgical_total > 0 ? (100 * data.surgical_pos / data.surgical_total).toFixed(0) : '-';
  const pmRate = data.postmortem_total > 0 ? (100 * data.postmortem_pos / data.postmortem_total).toFixed(0) : '-';
  const metaCards = `
    <div class="meta-card">
      <div class="lab">cohorts replicating</div>
      <div class="val">${{data.n_studies_this_class}} <small>cohorts (of up to 25 per class)</small></div>
    </div>
    <div class="meta-card">
      <div class="lab">donors in class</div>
      <div class="val">${{data.n_donors_this_class}}</div>
    </div>
    <div class="meta-card">
      <div class="lab">surgical agreement</div>
      <div class="val">${{data.surgical_pos}} / ${{data.surgical_total}} <small>(${{surgRate}}%)</small></div>
    </div>
    <div class="meta-card">
      <div class="lab">postmortem agreement</div>
      <div class="val">${{data.postmortem_pos}} / ${{data.postmortem_total}} <small>(${{pmRate}}%)</small></div>
    </div>
    <div class="meta-card">
      <div class="lab">specificity_delta</div>
      <div class="val">${{data.specificity_delta.toFixed(4)}}</div>
    </div>
    <div class="meta-card">
      <div class="lab">specificity_ratio</div>
      <div class="val">${{data.specificity_ratio.toFixed(2)}}×</div>
    </div>
  `;
  // Use chrom/start/end from the expansion data directly (NB peak_ids are peak_NNNN style, not chr:start-end)
  const pp = (data.chrom && data.start != null && data.end != null)
    ? {{ chrom: data.chrom, start: data.start, end: data.end }}
    : parsePeak(peak_id);
  let ucscLink = '';
  if (pp) {{
    const flank = 500;
    const ucscPos = `${{pp.chrom}}:${{pp.start - flank}}-${{pp.end + flank}}`;
    const peakWidth = pp.end - pp.start;
    const bedText = `track name="NB atlas peak" description="${{peak_id}} (${{tab_ct}})" color=139,92,246 visibility=pack itemRgb=on\n${{pp.chrom}}\t${{pp.start}}\t${{pp.end}}\tv3_${{tab_ct}}_peak\t1000\t.\t${{pp.start}}\t${{pp.end}}\t139,92,246`;
    const customText = encodeURIComponent(bedText);
    const ucscParams = [
      'db=hg38',
      `position=${{ucscPos}}`,
      `hgct_customText=${{customText}}`,
      'hgt.labelWidth=20',
      'encRegTfbsClusteredV3=dense',
      'encRegDnaseClustered=dense',
      'wgEncodeRegMarkH3k27ac=full',
      'phastCons100way=full',
      'jaspar2024=pack',
      'rmsk=dense'
    ].join('&');
    const ucscUrl = `https://genome.ucsc.edu/cgi-bin/hgTracks?${{ucscParams}}`;
    const aavNote = peakWidth > 700
      ? `<span class="aav-note">⚠ Peak is ${{peakWidth}} bp — larger than a typical ~500 bp AAV enhancer insert. See the Design B narrow sub-peaks below for AAV-cloneable fragments inside this region.</span>`
      : `<span class="aav-note">✓ Peak is ${{peakWidth}} bp — fits a standard ~500 bp AAV enhancer insert directly.</span>`;
    ucscLink = `<div class="ucsc-link">
      <a href="${{ucscUrl}}" target="_blank" rel="noopener noreferrer">🧬 View in UCSC Genome Browser</a>
      <span class="ucsc-ctx">— opens ${{ucscPos}} (peak + ±500 bp flanks). Our atlas peak is overlaid as a <span style="color:#8b5cf6;font-weight:600">purple block</span>. ENCODE / H3K27ac / phastCons / JASPAR / RepeatMasker pre-enabled.</span>
      <br>${{aavNote}}
    </div>`;
  }}

  // Design B narrow sub-peaks (AAV-cloneable fragments inside this wider v1.0.0 peak)
  let narrowBlock = '';
  const subs = NARROW_SUBPEAKS[peak_id];
  if (pp && subs && subs.length > 0) {{
    const parentW = pp.end - pp.start;
    const rows = subs.map(s => {{
      const w = s.w;
      const fit = w <= 700
        ? '<span class="aav-ok">✓ AAV-fit</span>'
        : (w <= 1500 ? '<span class="aav-mid">◐ trim to ~500</span>' : '<span class="aav-no">⚠ too wide</span>');
      const subUrl = `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${{pp.chrom}}:${{s.s - 100}}-${{s.e + 100}}&wgEncodeRegMarkH3k27ac=full&phastCons100way=full&jaspar2024=pack`;
      const dasUrl = `https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=${{pp.chrom}};start=${{s.s}};end=${{s.e}}`;
      return `<div class="sub-row">
        <span class="sub-coord">${{pp.chrom}}:${{s.s}}-${{s.e}}</span>
        <span class="sub-width">${{w}} bp</span>
        ${{fit}}
        <a href="${{subUrl}}" target="_blank" rel="noopener noreferrer">UCSC</a>
        <a href="${{dasUrl}}" target="_blank" rel="noopener noreferrer" title="UCSC REST API returns the hg38 sequence as JSON">get seq</a>
      </div>`;
    }}).join('');
    narrowBlock = `<div class="narrow-box">
      <h4>Design B narrow sub-peaks — AAV-cloneable fragments inside this ${{parentW}} bp wide A peak</h4>
      <div class="narrow-ctx">${{subs.length}} narrow MACS3 sub-peak(s) from the 2.06M-peak Design B union BED. Pick a ✓ AAV-fit sub-peak (≤700 bp) as a direct AAV insert, or trim a ◐ sub-peak to the ~500 bp core using the UCSC view (follow phastCons peaks + JASPAR TFBS clusters within our purple block).</div>
      <div class="sub-list">${{rows}}</div>
      <div class="narrow-cli">For full hg38 sequence of all sub-peaks + FASTA: <code>python3 scripts/narrow_lookup.py ${{pp.chrom}}:${{pp.start}}-${{pp.end}} --fasta</code></div>
    </div>`;
  }} else if (pp) {{
    narrowBlock = `<div class="narrow-box narrow-empty">
      <h4>Design B narrow sub-peaks</h4>
      <div class="narrow-ctx">No overlapping narrow Design B peak was found inside this A peak (unusual — the A peak may fall in a region the Design B MACS3 pass did not call). The wider A peak above is still your only AAV handle.</div>
    </div>`;
  }}

  return `<div class="expansion">
    ${{ucscLink}}
    ${{narrowBlock}}
    <div class="meta-grid">${{metaCards}}</div>
    <div class="ct-bar-pair">
      <div>
        <h4>Specificity — scaled to this peak's max (cell-type selectivity)</h4>
        <div class="ct-bars">${{ctBars}}</div>
      </div>
      <div>
        <h4>Absolute strength — % of this class's atlas-wide-max strength</h4>
        <div class="ct-bars">${{ctBarsAbs}}</div>
      </div>
    </div>
  </div>`;
}}

function renderTab(ct) {{
  currentCT = ct;
  const isSurgical = document.body.classList.contains('surgical-mode');
  const rows = TAB_DATA[ct] || [];
  const filtered = currentDistFilter === 'all'
    ? rows
    : rows.filter(r => r.distance_band === currentDistFilter);
  if (dt) {{ dt.destroy(); $('#atlas-table tbody').empty(); }}
  const tbody = document.getElementById('atlas-table').getElementsByTagName('tbody')[0];
  filtered.forEach(r => {{
    const tr = document.createElement('tr');
    tr.setAttribute('data-peak', r.peak_id);
    const strengthStr = r.strength_mean_evaluable_wn.toFixed(4);
    const detectionOrder = r.n_patients_detected;
    const detectionTitle = isSurgical
      ? `${{r.n_patients_detected}} of ${{r.n_patients_evaluable}} surgical evaluable donors (≥100 cells of this class) have ≥1% penetrance at this peak. Denominator is CONSTANT per class; donors with no signal count as detected=0, not unevaluable.`
      : `${{r.n_patients_detected}} of ${{r.n_patients_evaluable}} atlas donors (≥100 cells of this class) have ≥1% penetrance at this peak. GBM-style stable denominator per class.`;
    tr.innerHTML = `
      <td class="peakid"><span class="expander" title="Click to expand">▸</span> ${{r.peak_id}}</td>
      <td>${{r.chrom}}</td>
      <td class="num">${{r.start}}</td>
      <td class="num">${{r.end}}</td>
      <td class="num">${{r.peak_width_bp}}</td>
      <td>${{r.nearest_gene}}</td>
      <td class="num">${{fmtDist(r.dist_to_tss_signed)}}</td>
      <td>${{r.distance_band}}</td>
      <td class="num">${{r.log2FC_top_vs_next.toFixed(2)}}</td>
      <td class="num">${{r.spec_pct.toFixed(3)}}</td>
      <td class="num" title="strength_wn_gbm = mean over ALL EVALUABLE donors (≥100 cells of class) of fraction_open × (500 / peak_width_bp). Donors with no signal count as 0, not unevaluable. GBM v4-matched.">${{strengthStr}}</td>
      <td class="num">${{r.strength_pct.toFixed(3)}}</td>
      <td class="num primary-score-cell">${{r.primary_score.toFixed(3)}}</td>
      <td class="num">${{r.n_cohorts}}</td>
      <td class="num" title="${{detectionTitle}}" data-order="${{detectionOrder}}">${{fmtDetection(r)}}</td>`;
    tbody.appendChild(tr);
  }});
  dt = $('#atlas-table').DataTable({{
    paging: true, pageLength: 100,
    lengthMenu: [25, 50, 100, 250, 1000],
    order: [[12, 'desc']],
    dom: '<"top"lf>rt<"bottom"ip>',
    columnDefs: [
      {{ orderable: true, targets: '_all' }},
      {{ visible: false, targets: [1,2,3] }},
      {{ className: 'num', targets: [2,3,4,6,8,9,10,11,12,13,14] }}
    ]
  }});
  $('#atlas-table tbody').off('click').on('click', 'tr', function(e) {{
    const tr = $(this);
    const peak_id = tr.data('peak');
    const row = dt.row(tr);
    if (row.child.isShown()) {{
      row.child.hide();
      tr.find('.expander').text('▸');
    }} else {{
      row.child(renderExpansion(peak_id, currentCT)).show();
      tr.find('.expander').text('▾');
    }}
  }});
  $('.col-toggle').each(function() {{
    const idx = parseInt(this.dataset.col);
    this.checked = dt.column(idx).visible();
  }});
  document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
  document.querySelector(`.tab-btn[data-ct="${{ct}}"]`).classList.add('active');
}}

document.querySelectorAll('.tab-btn').forEach(btn => {{
  btn.addEventListener('click', () => renderTab(btn.dataset.ct));
}});
$(document).on('change', '.col-toggle', function() {{
  const col = dt.column(parseInt(this.dataset.col));
  col.visible(this.checked);
}});
document.getElementById('dist-filter').addEventListener('change', function() {{
  currentDistFilter = this.value;
  renderTab(currentCT);
}});

function updateTabCounts() {{
  document.querySelectorAll('.tab-btn').forEach(btn => {{
    const ct = btn.dataset.ct;
    const n = (TAB_DATA[ct] || []).length;
    const span = btn.querySelector('.n');
    if (span) span.textContent = `(${{n}})`;
    btn.classList.toggle('empty-ct', n === 0);
  }});
}}

document.getElementById('surgical-only').addEventListener('change', function() {{
  TAB_DATA = this.checked ? TAB_DATA_SURGICAL : TAB_DATA_DEFAULT;
  document.body.classList.toggle('surgical-mode', this.checked);
  updateTabCounts();
  // Keep current cell type if it still has rows; else fall back to first non-empty
  const stillHasRows = (TAB_DATA[currentCT] || []).length > 0;
  const fallback = CT_LIST.find(ct => (TAB_DATA[ct] || []).length > 0) || currentCT;
  renderTab(stillHasRows ? currentCT : fallback);
}});

function updateStickyHeadTop() {{
  const header = document.querySelector('header');
  const tabs = document.querySelector('.tabs');
  if (!header || !tabs) return;
  const headerH = header.offsetHeight;
  const tabsH = tabs.offsetHeight;
  document.documentElement.style.setProperty('--header-height', headerH + 'px');
  document.documentElement.style.setProperty('--sticky-head-top', (headerH + tabsH) + 'px');
}}
window.addEventListener('resize', updateStickyHeadTop);
window.addEventListener('load', updateStickyHeadTop);
updateStickyHeadTop();

renderTab('{default_tab}');
setTimeout(updateStickyHeadTop, 100);
</script>
</body></html>
"""


if __name__ == "__main__":
    main()
