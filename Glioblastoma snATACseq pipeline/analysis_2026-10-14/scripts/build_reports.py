"""Build annotated reports with patient-level + off-target metrics.

Data dependency: expects the per-report scored parquets under reports/detail/
and reports/detail/per_celltype/. Download them (along with the pseudobulk
matrix) from the Zenodo deposit: https://doi.org/10.5281/zenodo.23252324

For each per-cell-type/pan report in reports/detail/:
  - Annotate with nearest_gene, dist_to_tss_signed, distance_category
  - Write annotated.csv (all peaks sorted by daigle_z)
  - top100.csv (any-distance)
  - top50_distal.csv (|dist_to_tss| > 2 kb)
  - top30_distal.bed
  - .html report with sortable columns including the new metrics

Also build:
  - top50_combined.{csv,parquet,xlsx} across all reports
  - top50_distal_combined.{csv,parquet,xlsx}
  - index.html gallery

Then move reports/detail/ + the annotated outputs into reports/ (replacing v5 reports).
"""
import re
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "reports/detail"
OUT = ROOT / "reports" / "detail"
PER_CT = OUT / "per_celltype"
OUT.mkdir(parents=True, exist_ok=True)
PER_CT.mkdir(parents=True, exist_ok=True)
REFGENE_BED = ROOT / "reference" / "refgene_hg38.bed"

MIN_PATIENTS = 2  # replication filter: require ≥2 patients
TOP_N_FULL = 100
TOP_N_DISTAL = 50
TOP_N_BED = 30
PROMOTER_WINDOW = 2000
PROXIMAL_WINDOW = 10000
DISTAL_WINDOW = 100000


def load_tss_per_chrom():
    g = pd.read_csv(REFGENE_BED, sep="\t", header=None,
                    names=["chrom", "start", "end", "gene", "score", "strand"])
    g["tss"] = np.where(g["strand"] == "+", g["start"], g["end"])
    g = g[["chrom", "tss", "gene"]].drop_duplicates().sort_values(["chrom", "tss"]).reset_index(drop=True)
    return {ch: group.reset_index(drop=True) for ch, group in g.groupby("chrom")}


def nearest_tss(peak_df, tss_by_chrom):
    out_gene = np.full(len(peak_df), "", dtype=object)
    out_dist = np.zeros(len(peak_df), dtype=np.int64)
    peak_df = peak_df.reset_index(drop=True)
    for ch, grp in peak_df.groupby("chrom"):
        if ch not in tss_by_chrom:
            continue
        tss = tss_by_chrom[ch]
        tss_arr = tss["tss"].values
        gene_arr = tss["gene"].values
        centers = ((grp["start"].values + grp["end"].values) // 2).astype(np.int64)
        orig_i = grp.index.values
        ins = np.searchsorted(tss_arr, centers)
        prev_i = np.clip(ins - 1, 0, len(tss_arr) - 1)
        next_i = np.clip(ins, 0, len(tss_arr) - 1)
        d_prev = centers - tss_arr[prev_i]
        d_next = centers - tss_arr[next_i]
        pick_prev = np.abs(d_prev) <= np.abs(d_next)
        chosen_d = np.where(pick_prev, d_prev, d_next)
        chosen_i = np.where(pick_prev, prev_i, next_i)
        out_dist[orig_i] = chosen_d
        out_gene[orig_i] = gene_arr[chosen_i]
    peak_df = peak_df.copy()
    peak_df["nearest_gene"] = out_gene
    peak_df["dist_to_tss_signed"] = out_dist
    abs_d = np.abs(peak_df["dist_to_tss_signed"].values)
    peak_df["distance_category"] = np.where(abs_d <= PROMOTER_WINDOW, "promoter",
                                   np.where(abs_d <= PROXIMAL_WINDOW, "proximal",
                                   np.where(abs_d <= DISTAL_WINDOW, "distal", "far-distal")))
    return peak_df


def fmt_dist(d):
    if abs(d) >= 1000:
        return f"{d/1000:+.1f}kb"
    return f"{d:+d}"


def write_bed(df, path):
    s = df["daigle_z_arm_adjusted"].values
    if len(s):
        s = np.clip(((s - s.min()) / max(s.max() - s.min(), 1e-9) * 1000).astype(int), 0, 1000)
    else:
        s = np.zeros(len(df), dtype=int)
    name = [f"{g}_{fmt_dist(d)}_{p}" for g, d, p in
            zip(df["nearest_gene"], df["dist_to_tss_signed"], df["peak_id"])]
    bed = pd.DataFrame({
        "chrom": df["chrom"], "start": df["start"], "end": df["end"],
        "name": name, "score": s, "strand": "."
    })
    bed.to_csv(path, sep="\t", header=False, index=False)


HTML_TEMPLATE = """<!doctype html>
<html><head><meta charset="utf-8"><title>v7 — {title}</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:1800px;margin:1em auto;padding:0 1em;color:#222}}
h1{{color:#333;margin-bottom:0.2em}}
.meta{{background:#f4f4f4;padding:0.8em;border-radius:4px;margin:1em 0;font-size:11pt}}
table{{border-collapse:collapse;width:100%;font-size:10pt}}
th{{background:#333;color:white;padding:5px 6px;text-align:left;cursor:pointer;position:sticky;top:0;font-size:9pt}}
td{{padding:3px 6px;border-bottom:1px solid #eee;font-family:ui-monospace,monospace;font-size:9pt}}
tr:hover{{background:#f8f8f8}}
.promoter{{color:#999}}
.distal{{color:#0a7;font-weight:bold}}
.gene{{color:#06c;font-weight:bold}}
.batch{{background:#ffcc99}}
.leak{{background:#ffd0d0}}
</style></head><body>
<h1>v7 — {title}</h1>
<div class="meta">
<b>Primary replication:</b> n_target_patients ≥ {min_patients} filter applied.<br>
<b>Specificity:</b> daigle_z_arm_adjusted = (target_mean_cpm − bg_mean_cpm) / bg_std, with CPM floor 0.5 before ranking.<br>
<b>AAV safety:</b> offtarget_max_cpm = worst non-target cell type signal; aav_safety_ratio = target / (offtarget_max + 0.1); offtarget_top_cell_type = which non-target lineage leaks most.<br>
<b>Batch flags:</b> cohort_balance close to 1 means one cohort dominates target signal (highlighted orange ≥ 0.9).<br>
Showing top 500 by Daigle Z. Click column headers to re-sort.
</div>
<table id="t"><thead>
<tr><th>rank</th><th>peak_id</th><th>gene</th><th>dist</th><th>cat</th>
<th>z_arm_adj</th><th>z_raw</th><th>target_cpm</th>
<th>n_pat</th><th>n_cohorts</th><th>cohort_bal</th>
<th>off_max</th><th>off_top</th><th>aav_safety</th><th>n_off>1</th>
<th>cnv_class</th><th>arm</th><th>%pat_chr7_elev</th></tr>
</thead><tbody>
{rows}
</tbody></table>
<script>
document.querySelectorAll('th').forEach((th,i)=>{{th.addEventListener('click',()=>{{
  const rows=Array.from(document.querySelectorAll('#t tbody tr'));
  const asc=th.dataset.sort!=='asc';th.dataset.sort=asc?'asc':'desc';
  rows.sort((a,b)=>{{const av=a.children[i].innerText.replace(/,/g,'');const bv=b.children[i].innerText.replace(/,/g,'');
  const an=parseFloat(av),bn=parseFloat(bv);
  if(!isNaN(an)&&!isNaN(bn))return asc?an-bn:bn-an;return asc?av.localeCompare(bv):bv.localeCompare(av);}});
  const tb=document.querySelector('#t tbody');rows.forEach(r=>tb.appendChild(r));
}})}});
</script></body></html>"""


def html_rows(df):
    rows = []
    for i, r in df.iterrows():
        cat = r["distance_category"]
        css_cat = "promoter" if cat == "promoter" else ("distal" if cat in ("distal","far-distal") else "")
        bal = r["cohort_balance"]
        css_bal = ' class="batch"' if bal >= 0.9 else ""
        leak = r["offtarget_max_cpm"]
        css_leak = ' class="leak"' if leak > 5 else ""
        pfc7 = r.get("patient_fraction_chr7_elevated")
        pfc7_txt = "" if pd.isna(pfc7) else f"{pfc7:.0%}"
        rows.append(
            f'<tr><td>{i+1}</td>'
            f'<td>{r["peak_id"]}</td>'
            f'<td class="gene">{r["nearest_gene"]}</td>'
            f'<td>{fmt_dist(int(r["dist_to_tss_signed"]))}</td>'
            f'<td class="{css_cat}">{cat}</td>'
            f'<td>{r["daigle_z_arm_adjusted"]:.2f}</td>'
            f'<td>{r.get("daigle_z_raw", 0):.1f}</td>'
            f'<td>{r["target_mean_cpm"]:.2f}</td>'
            f'<td>{int(r["n_target_patients"])}</td>'
            f'<td>{int(r["n_target_cohorts"])}</td>'
            f'<td{css_bal}>{bal:.2f}</td>'
            f'<td{css_leak}>{leak:.2f}</td>'
            f'<td>{r["offtarget_top_cell_type"]}</td>'
            f'<td>{r["aav_safety_ratio"]:.1f}</td>'
            f'<td>{int(r["n_offtarget_ct_above_1"])}</td>'
            f'<td>{r.get("chromosome_cnv_class", "")}</td>'
            f'<td>{r.get("chromosome_arm", "")}</td>'
            f'<td>{pfc7_txt}</td>'
            f'</tr>'
        )
    return "\n".join(rows)


def process_one(parquet_path, tss_by_chrom):
    df = pd.read_parquet(parquet_path)
    n_total = len(df)
    df = df[df["n_target_patients"] >= MIN_PATIENTS].copy()
    df = df.sort_values("daigle_z_arm_adjusted", ascending=False).reset_index(drop=True)
    df = nearest_tss(df, tss_by_chrom)

    stem = parquet_path.stem
    is_pan = stem.startswith("pan_")
    out_dir = OUT if is_pan else PER_CT
    out_dir.mkdir(exist_ok=True)

    df.to_csv(out_dir / f"{stem}_annotated.csv", index=False)
    df.head(TOP_N_FULL).to_csv(out_dir / f"{stem}_top100.csv", index=False)
    distal = df[df["distance_category"] != "promoter"].reset_index(drop=True)
    distal.head(TOP_N_DISTAL).to_csv(out_dir / f"{stem}_top{TOP_N_DISTAL}_distal.csv", index=False)
    write_bed(distal.head(TOP_N_BED), out_dir / f"{stem}_top{TOP_N_BED}_distal.bed")

    html_df = df.head(500).reset_index(drop=True)
    html = HTML_TEMPLATE.format(title=stem, min_patients=MIN_PATIENTS, rows=html_rows(html_df))
    (out_dir / f"{stem}.html").write_text(html)

    print(f"  {stem}: {n_total:,} total → {len(df):,} after n_patients≥{MIN_PATIENTS}; "
          f"top-{TOP_N_BED} distal: {TOP_N_BED} non-promoter")
    # Copy parquet to canonical location only if source differs from destination
    dest = out_dir / f"{stem}.parquet"
    if parquet_path.resolve() != dest.resolve():
        shutil.copy(parquet_path, dest)
    return df


def build_combined(all_dfs, out_name, filter_distal=False):
    rows = []
    for stem, df in all_dfs.items():
        label = stem.replace("per_celltype_", "")
        if filter_distal:
            df = df[df["distance_category"] != "promoter"]
        top = df.head(50).copy()
        top.insert(0, "report", label)
        top.insert(1, "rank", range(1, len(top) + 1))
        rows.append(top)
    combined = pd.concat(rows, ignore_index=True)
    combined.to_parquet(OUT / f"{out_name}.parquet", index=False)
    combined.to_csv(OUT / f"{out_name}.csv", index=False)
    with pd.ExcelWriter(OUT / f"{out_name}.xlsx", engine="openpyxl") as xw:
        summary = combined.groupby("report").agg(
            n=("rank", "count"),
            median_z=("daigle_z_arm_adjusted", "median"),
            median_n_patients=("n_target_patients", "median"),
            median_cohort_balance=("cohort_balance", "median"),
            median_offtarget_max=("offtarget_max_cpm", "median"),
            median_aav_safety=("aav_safety_ratio", "median"),
        ).reset_index()
        summary.to_excel(xw, sheet_name="summary", index=False)
        for label in combined["report"].unique():
            combined[combined["report"] == label].to_excel(xw, sheet_name=label[:31], index=False)
    print(f"wrote {OUT / out_name}.{{csv,parquet,xlsx}} — {combined.shape}")


def build_index(all_dfs):
    rows = []
    for stem, df in all_dfs.items():
        is_pan = stem.startswith("pan_")
        sub = "." if is_pan else "per_celltype"
        top5 = df.head(5)[["nearest_gene", "daigle_z_arm_adjusted", "n_target_patients", "n_target_cohorts",
                           "cohort_balance", "offtarget_top_cell_type", "aav_safety_ratio"]]
        top5_html = top5.to_html(index=False, border=0, float_format=lambda x: f"{x:.2f}")
        rows.append(f'<h2><a href="{sub}/{stem}.html">{stem}</a></h2>{top5_html}')
    html = f"""<!doctype html><html><head><meta charset="utf-8"><title>v7 reports</title>
<style>body{{font-family:system-ui;max-width:1300px;margin:1em auto;padding:0 1em}}table{{border-collapse:collapse}}td,th{{padding:3px 7px;border-bottom:1px solid #eee;font-size:9pt}}</style></head>
<body><h1>GBM Enhancer Atlas v7 — Reports</h1>
<p>Peak set: 1.96M union. Scoring: Daigle Z + patient-level replication + per-non-target-cell-type contamination.
Filter: n_target_patients ≥ {MIN_PATIENTS}, target_mean_cpm ≥ 0.5.</p>
<p><b>What to look at per peak:</b> daigle_z (specificity), n_target_patients (biological replication),
cohort_balance (batch sanity, 0.9+ = one cohort dominates), offtarget_top_cell_type (worst non-target leaker),
aav_safety_ratio (target / offtarget_max).</p>
{"".join(rows)}
</body></html>"""
    (OUT / "index.html").write_text(html)
    print(f"wrote {OUT / 'index.html'}")


def main():
    print(f"loading refGene from {REFGENE_BED}")
    tss_by_chrom = load_tss_per_chrom()
    print(f"  {sum(len(v) for v in tss_by_chrom.values()):,} TSS entries")

    all_dfs = {}
    for parquet in sorted(PER_CT.glob("per_celltype_*.parquet")) + sorted(RAW.glob("pan_*.parquet")):
        df = process_one(parquet, tss_by_chrom)
        all_dfs[parquet.stem] = df

    print()
    build_combined(all_dfs, "top50_combined", filter_distal=False)
    build_combined(all_dfs, "top50_distal_combined", filter_distal=True)
    build_index(all_dfs)
    print("done.")


if __name__ == "__main__":
    sys.exit(main() or 0)
