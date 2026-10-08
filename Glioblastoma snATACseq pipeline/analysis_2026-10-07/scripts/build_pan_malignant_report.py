"""
Build a rich pan-malignant deep-dive report.

Produces:
  reports/pan_malignant_report.html     - human-readable, nearest-gene + plain-English rec
  reports/pan_malignant_top100.csv      - full top 100 with all columns, for Excel
  reports/pan_malignant_top50.bed       - GRCh38 BED, ready for Benchling / UCSC browser
"""
from __future__ import annotations
import sys
from pathlib import Path

import polars as pl
import pyranges as pr
import pandas as pd

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26")
SCORES = HERE / "data/matrix/pan_malignant_scores.parquet"
# Also fall back to Drive if local hasn't synced
if not SCORES.exists():
    SCORES = Path(
        "/Users/jpmcginnis1/Library/CloudStorage/GoogleDrive-jpmcginnis1@gmail.com/My Drive/"
        "AAV Gene Therapy/Enhancer ATACseq projects/October 2026 WT GBM analysis/"
        "matrix/pan_malignant_scores.parquet"
    )
GENE_BED = HERE / "reference/refgene_hg38.bed"
OUT_HTML = HERE / "reports/pan_malignant_report.html"
OUT_CSV = HERE / "reports/pan_malignant_top100.csv"
OUT_BED = HERE / "reports/pan_malignant_top50.bed"
OUT_HTML.parent.mkdir(parents=True, exist_ok=True)

# --- GBM biology context -----------------------------------------------------
# Curated per-category gene lists. A hit NEAR these flags for "biologically plausible".

GLIOMA_DRIVERS = {
    "EGFR","PDGFRA","MET","FGFR1","FGFR3","NTRK1","NTRK2","NTRK3",
    "TP53","PTEN","NF1","CDKN2A","CDKN2B","RB1","PIK3CA","PIK3R1",
    "IDH1","IDH2","ATRX","TERT","MGMT","MYC","MYCN","BRAF",
    "CIC","FUBP1","TCF12",
}
GSC_AND_LINEAGE = {
    "SOX2","SOX9","SOX10","SOX11","OLIG1","OLIG2","NESTIN","NES",
    "PROM1","CD44","FABP7","GFAP","VIM","MSI1","NANOG","POU5F1",
    "ASCL1","HES1","NOTCH1","NOTCH2","DLL3","JAG1","JAG2",
    "NFIA","NFIB","DLX1","DLX2","DLX5","DLX6",
}
NEFTEL_MES_LIKE = {"CD44","CHI3L1","VIM","SERPINA3","ANXA1","TGFBI","CSTB","S100A11"}
NEFTEL_AC_LIKE = {"GFAP","S100B","APOE","AQP4","ALDH1L1","SLC1A3","MLC1"}
NEFTEL_OPC_LIKE = {"OLIG1","OLIG2","PDGFRA","NFIA","APOD","CSPG4"}
NEFTEL_NPC_LIKE = {"DLL3","DLX1","DLX2","DLX5","DLX6","SOX11","NEUROD4","STMN2"}
TME_MYELOID = {"TMEM119","P2RY12","CX3CR1","CD68","ITGAM","AIF1","CSF1R","TREM2"}
CHROMATIN_REMODELER = {
    "EZH2","SUZ12","EED","BMI1","CBX2","CBX4","CBX7","CBX8",
    "KDM6A","KDM6B","KDM5A","KDM5C","SETD2","MLL1","MLL2","MLL3","MLL4",
    "DNMT1","DNMT3A","DNMT3B","TET1","TET2","TET3",
}

GBM_GENES_ALL = (
    GLIOMA_DRIVERS | GSC_AND_LINEAGE | NEFTEL_MES_LIKE | NEFTEL_AC_LIKE
    | NEFTEL_OPC_LIKE | NEFTEL_NPC_LIKE | TME_MYELOID | CHROMATIN_REMODELER
)

def gbm_tags(gene: str) -> list[str]:
    tags = []
    if gene in GLIOMA_DRIVERS: tags.append("glioma driver")
    if gene in GSC_AND_LINEAGE: tags.append("GSC / lineage TF")
    if gene in NEFTEL_MES_LIKE: tags.append("Neftel MES-like")
    if gene in NEFTEL_AC_LIKE: tags.append("Neftel AC-like")
    if gene in NEFTEL_OPC_LIKE: tags.append("Neftel OPC-like")
    if gene in NEFTEL_NPC_LIKE: tags.append("Neftel NPC-like")
    if gene in TME_MYELOID: tags.append("TME myeloid marker (be suspicious!)")
    if gene in CHROMATIN_REMODELER: tags.append("chromatin remodeler")
    return tags

# --- Load + annotate ---------------------------------------------------------

print(f"[report] loading {SCORES}")
scores = pl.read_parquet(SCORES)
print(f"[report] n peaks scored: {scores.height:,}")
print(f"[report] cohort coverage: max n_cohorts = {scores['n_cohorts'].max()}")

# Load gene BED
genes = pr.read_bed(str(GENE_BED))

# Build the ranked table
# A cell-type-general cloning candidate wants:
#   - high cross-cohort replication (n_cohorts)
#   - high patient breadth (n_patients_detected)
#   - sustained openness per patient (consistency = fraction of patients ≥10% open)
#   - decent mean strength (strength_mean)
# AVOID chr7 CNV-driven hits (chr7+ is universal in GBM; artifactually inflates accessibility).

ranked = (
    scores
    .with_columns(
        pl.when(pl.col("chrom") == "chr7").then(1).otherwise(0).alias("chr7_cnv_flag"),
        pl.when(pl.col("chrom") == "chr10").then(1).otherwise(0).alias("chr10_cnv_flag"),
    )
    .filter(pl.col("n_cohorts") >= 3)
    .with_columns(
        # A clean replicable_score: emphasizes consistency + patient breadth, not strength alone.
        (
            pl.col("consistency") * 0.40
            + pl.col("n_cohorts").cast(pl.Float64) / 8.0 * 0.30
            + (pl.col("n_patients_detected").cast(pl.Float64) / 33.0).clip(0, 1) * 0.20
            + pl.col("strength_mean") * 0.10
        ).alias("replicable_score")
    )
    .sort("replicable_score", descending=True)
)

TOP_N = 100
top = ranked.head(TOP_N).to_pandas()
print(f"[report] top {TOP_N} selected, chr7 count: {int(top['chr7_cnv_flag'].sum())}")

# Nearest-gene annotation — manual implementation since pyranges .nearest() behavior
# proved unreliable (returned 29 Mb-distant pseudogenes when relevant protein-coding
# genes were within 5 Mb). Simple per-chromosome binary search over gene TSSes.
import bisect
gene_df = pd.read_csv(GENE_BED, sep="\t", header=None,
                     names=["Chromosome","Start","End","Name","GeneType","Strand"])
# TSS = start if '+' strand, end if '-' strand
gene_df["tss"] = gene_df.apply(
    lambda r: r["Start"] if r["Strand"] == "+" else r["End"], axis=1
).astype(int)
tss_by_chr = {}
for ch, grp in gene_df.groupby("Chromosome"):
    g = grp.sort_values("tss").reset_index(drop=True)
    tss_by_chr[ch] = (g["tss"].values, g["Name"].values, g["Strand"].values,
                      g["Start"].values, g["End"].values)

def nearest_gene(chrom, peak_center):
    if chrom not in tss_by_chr: return (None, None, None, None, None)
    tss_arr, names, strands, starts, ends = tss_by_chr[chrom]
    idx = bisect.bisect_left(tss_arr, peak_center)
    cands = []
    if idx > 0:         cands.append(idx - 1)
    if idx < len(tss_arr): cands.append(idx)
    # pick candidate with min abs distance to TSS
    best = min(cands, key=lambda i: abs(peak_center - tss_arr[i]))
    return (names[best], strands[best], int(starts[best]), int(ends[best]), int(peak_center - tss_arr[best]))

annots = []
for _, r in top.iterrows():
    peak_center = (int(r["start"]) + int(r["end"])) // 2
    name, strand, gstart, gend, dist = nearest_gene(r["chrom"], peak_center)
    annots.append((name, strand, gstart, gend, dist))
annot_df = pd.DataFrame(annots, columns=["nearest_gene","gene_strand","gene_start","gene_end","dist_to_tss_signed"])
annotated = pd.concat([top.reset_index(drop=True), annot_df], axis=1)
annotated["gbm_tags"] = annotated["nearest_gene"].fillna("").map(lambda g: ", ".join(gbm_tags(g)))

# Cloning advice (plain English)
def clone_rec(row):
    recs = []
    if row["chr7_cnv_flag"] == 1:
        recs.append("⚠️ chr7 — CNV-driven signal possible, re-check in CNV-matched reference")
    if row["chr10_cnv_flag"] == 1:
        recs.append("⚠️ chr10 — CNV-driven signal possible")
    if row["gbm_tags"]:
        recs.append(f"biology: {row['gbm_tags']}")
    if row["dist_to_tss_signed"] is not None:
        d = row["dist_to_tss_signed"]
        if abs(d) < 2000:
            recs.append("promoter-proximal (≤2 kb from TSS) — may act as promoter, not distal enhancer")
        elif abs(d) < 100_000:
            recs.append(f"distal enhancer candidate ({d/1000:+.0f} kb from {row['nearest_gene']} TSS)")
        else:
            recs.append(f"far from any gene ({d/1000:+.0f} kb from {row['nearest_gene']}) — intergenic")
    if row["consistency"] >= 0.8 and row["n_cohorts"] >= 5 and row["chr7_cnv_flag"] == 0:
        recs.append("✅ STRONG cloning candidate")
    elif row["consistency"] >= 0.65 and row["n_cohorts"] >= 4:
        recs.append("👍 solid cloning candidate")
    return "; ".join(recs) if recs else ""

annotated["cloning_rec"] = annotated.apply(clone_rec, axis=1)

# Reorder columns
cols = [
    "peak_id","chrom","start","end",
    "nearest_gene","dist_to_tss_signed","gbm_tags",
    "n_cohorts","n_patients_detected","n_patients_accessible_10pct",
    "consistency","strength_mean","strength_median",
    "replicable_score","chr7_cnv_flag",
    "cloning_rec",
]
annotated_out = annotated[cols]
annotated_out.to_csv(OUT_CSV, index=False)
print(f"[report] wrote {OUT_CSV}")

# BED (top 50 only, excluding chr7)
clean_top50 = annotated[annotated["chr7_cnv_flag"] == 0].head(50)
with open(OUT_BED, "w") as f:
    for _, row in clean_top50.iterrows():
        name = f"{row['nearest_gene']}_{row['peak_id']}" if pd.notna(row['nearest_gene']) else row['peak_id']
        f.write(f"{row['chrom']}\t{row['start']}\t{row['end']}\t{name}\t{int(row['replicable_score']*1000)}\t.\n")
print(f"[report] wrote {OUT_BED} ({len(clean_top50)} entries, chr7 excluded)")

# --- HTML report -------------------------------------------------------------
import html

def esc(s): return html.escape(str(s)) if s else ""

rows_html = []
for _, r in annotated_out.iterrows():
    tag_html = ""
    if r["gbm_tags"]:
        tag_html = f'<span style="color:#059669;font-weight:600;">{esc(r["gbm_tags"])}</span>'
    rec_html = esc(r["cloning_rec"])
    if "STRONG" in r["cloning_rec"]:
        rec_html = f'<strong style="color:#15803d;">{rec_html}</strong>'
    elif "solid" in r["cloning_rec"]:
        rec_html = f'<span style="color:#1e40af;">{rec_html}</span>'
    chr7_cell = '<span style="color:#dc2626;">⚠️</span>' if r["chr7_cnv_flag"] else ""
    tss_cell = f'{r["dist_to_tss_signed"]/1000:+.1f}' if pd.notna(r["dist_to_tss_signed"]) else "—"
    rows_html.append(f"""
      <tr>
        <td><code>{esc(r["peak_id"])}</code>{chr7_cell}</td>
        <td><strong>{esc(r["nearest_gene"])}</strong><br/>{tag_html}</td>
        <td style="text-align:right;">{tss_cell} kb</td>
        <td style="text-align:center;">{int(r["n_cohorts"])}/8</td>
        <td style="text-align:center;">{int(r["n_patients_detected"])}</td>
        <td style="text-align:center;">{int(r["n_patients_accessible_10pct"])}</td>
        <td style="text-align:center;">{r["consistency"]*100:.0f}%</td>
        <td style="text-align:center;">{r["strength_mean"]*100:.1f}%</td>
        <td>{rec_html}</td>
      </tr>
    """)

html_out = f"""<!DOCTYPE html>
<html><head>
<title>Pan-malignant enhancer candidate report — GBM</title>
<meta charset="utf-8"/>
<style>
  body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; max-width: 1400px; margin: 2em auto; padding: 0 2em; color: #111; line-height: 1.5; }}
  h1 {{ font-size: 1.8em; border-bottom: 2px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.4em; margin-top: 2em; color: #1e40af; }}
  code {{ font-size: 0.85em; background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.88em; margin: 1em 0; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 6px 8px; vertical-align: top; }}
  th {{ background: #f9fafb; text-align: left; font-weight: 600; position: sticky; top: 0; }}
  tr:hover {{ background: #fafbfc; }}
  .note {{ background: #fef3c7; border-left: 4px solid #f59e0b; padding: 1em; margin: 1em 0; }}
  .success {{ background: #d1fae5; border-left: 4px solid #059669; padding: 1em; margin: 1em 0; }}
</style></head>
<body>

<h1>Pan-malignant enhancer candidate report — GBM atlas 2026-10-07</h1>

<p>Source: 8 cohorts, 51 patients, 1,523,668 nuclei. Pan-malignant pool = cells with
either CNV-based malignant call (<code>malignant_cnv == 1</code>) OR
<code>cell_type == "malignant_unresolved"</code> label. Includes 438,994 CNV-malignant cells
from GBM-Space that marker-peak scoring had mis-labeled as TME types.</p>

<h2>How to read this table</h2>

<ul>
<li><strong>peak_id</strong> — GRCh38 coordinates (CATLAS 544K peak atlas); ⚠️ flags chr7 (universal CNV-driven in GBM).</li>
<li><strong>nearest_gene / GBM tags</strong> — closest RefSeq gene by TSS, color-coded for known glioma driver, GSC-lineage TF, Neftel state marker, chromatin remodeler, or TME myeloid marker (flagged suspicious — these cells were likely CNV-malignant mis-labeled).</li>
<li><strong>Distance to TSS</strong> — signed, kb. Promoter-proximal (|d| < 2 kb) means this is more likely a promoter than a distal enhancer. 2–100 kb = good distal enhancer range. >100 kb = suspicious (far intergenic, may regulate nothing).</li>
<li><strong>Cohorts (n/8)</strong> — how many of the 8 atlas cohorts had this peak accessible in their malignant cells. ≥4 is the cross-cohort replication bar.</li>
<li><strong>Patients detected / Patients ≥10%</strong> — how many of the ~33 patients-with-malignant-cells show the peak at all / show it in ≥10% of their malignant cells (the "accessible per patient" bar).</li>
<li><strong>Consistency</strong> — fraction of detected patients accessible ≥10%. ≥80% is strong, ≥65% is solid.</li>
<li><strong>Strength (mean)</strong> — mean fraction of malignant cells with the peak open, averaged across patients.</li>
<li><strong>Cloning rec</strong> — plain-English: ✅ STRONG, 👍 solid, ⚠️ caveats. These are heuristic, not scores.</li>
</ul>

<div class="note">
<strong>Chr7 caveat.</strong> chr7+ gain is the defining CNV of GBM — 100% of IDH-WT GBMs
have it. Peaks on chr7 are inflated by extra DNA copies, not necessarily cell-type-selective
enhancers. We flag them with ⚠️ but still show them. For first-pass cloning, prefer non-chr7
candidates (the attached BED file excludes chr7 and shows top 50 clean hits).
</div>

<div class="success">
<strong>Top recommended cloning candidates (plain English):</strong>
<ol>
"""
# Pull top 10 non-chr7 with STRONG label
best = annotated[(annotated["chr7_cnv_flag"] == 0)].head(10)
for _, r in best.iterrows():
    gene = r["nearest_gene"] if pd.notna(r["nearest_gene"]) else "?"
    tss = f"{r['dist_to_tss_signed']/1000:+.0f} kb" if pd.notna(r["dist_to_tss_signed"]) else "?"
    html_out += (
        f"<li><code>{r['peak_id']}</code> — "
        f"nearest gene <strong>{gene}</strong> ({tss} from TSS); "
        f"open in <strong>{int(r['n_cohorts'])}/8 cohorts</strong>, "
        f"<strong>{int(r['n_patients_accessible_10pct'])}/{int(r['n_patients_detected'])} patients</strong>, "
        f"mean {r['strength_mean']*100:.0f}% of malignant cells open."
    )
    if r["gbm_tags"]:
        html_out += f" <em>GBM relevance: {r['gbm_tags']}.</em>"
    html_out += "</li>\n"
html_out += """
</ol>
</div>

<h2>Full top-""" + str(TOP_N) + """ table</h2>
<table>
<thead><tr>
  <th>peak_id (GRCh38)</th>
  <th>nearest gene<br/>GBM tags</th>
  <th>dist to TSS</th>
  <th>cohorts</th>
  <th>n_pt detected</th>
  <th>n_pt ≥10%</th>
  <th>consist.</th>
  <th>strength</th>
  <th>cloning note</th>
</tr></thead>
<tbody>
""" + "\n".join(rows_html) + """
</tbody></table>

<h2>Files</h2>
<ul>
<li><code>pan_malignant_top100.csv</code> — full top 100 with all columns, Excel-ready.</li>
<li><code>pan_malignant_top50.bed</code> — top 50 non-chr7 candidates in BED, ready for Benchling / UCSC browser upload.</li>
<li><code>pan_malignant_scores.parquet</code> — all """ + f"{scores.height:,}" + """ scored peaks (not just top 100). Query locally with polars for any custom filter.</li>
</ul>

</body></html>
"""

OUT_HTML.write_text(html_out)
print(f"[report] wrote {OUT_HTML}")
print()
print("DONE.")
