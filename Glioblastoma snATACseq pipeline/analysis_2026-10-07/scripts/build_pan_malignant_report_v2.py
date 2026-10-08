"""
Pan-malignant enhancer candidate report v2.

Adds:
  - distal filter (>=2 kb from any TSS; promoters excluded)
  - selectivity filter (malignant strength / max non-malignant strength from candidate_scores.parquet)
  - "clean" ranked table that's actually usable for AAV cell-type-selective design

Produces four files in reports/:
  pan_malignant_report_v2.html
  pan_malignant_top50_distal_selective.csv   -- the AAV cloning shortlist
  pan_malignant_top100_all.csv               -- top 100 including promoters/non-selective (for ref)
  pan_malignant_top30_distal_selective.bed   -- GRCh38 BED, Benchling-ready
"""
from __future__ import annotations
import bisect, html, sys
from pathlib import Path

import pandas as pd
import polars as pl

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26")
GDRIVE = Path("/Users/jpmcginnis1/Library/CloudStorage/GoogleDrive-jpmcginnis1@gmail.com/"
              "My Drive/AAV Gene Therapy/Enhancer ATACseq projects/October 2026 WT GBM analysis")

PAN_SCORES = GDRIVE / "matrix/pan_malignant_scores.parquet"
# candidate_scores on Drive is still mid-sync (multipart .xxx chunks); use local complete copy
CELLTYPE_SCORES = HERE / "data/matrix/candidate_scores.parquet"
GENE_BED = HERE / "reference/refgene_hg38.bed"
OUT_HTML = HERE / "reports/pan_malignant_report_v2.html"
OUT_CSV_CLEAN = HERE / "reports/pan_malignant_top50_distal_selective.csv"
OUT_CSV_ALL = HERE / "reports/pan_malignant_top100_all.csv"
OUT_BED = HERE / "reports/pan_malignant_top30_distal_selective.bed"
OUT_HTML.parent.mkdir(parents=True, exist_ok=True)


# --- GBM biology context -----------------------------------------------------
GLIOMA_DRIVERS = {"EGFR","PDGFRA","MET","FGFR1","FGFR3","NTRK1","NTRK2","NTRK3","TP53","PTEN","NF1","CDKN2A","CDKN2B","RB1","PIK3CA","PIK3R1","IDH1","IDH2","ATRX","TERT","MGMT","MYC","MYCN","BRAF","CIC","FUBP1","TCF12"}
GSC_LINEAGE = {"SOX2","SOX9","SOX10","SOX11","OLIG1","OLIG2","NES","PROM1","CD44","FABP7","GFAP","VIM","MSI1","ASCL1","HES1","HES5","NOTCH1","NOTCH2","DLL3","JAG1","NFIA","NFIB","DLX1","DLX2","DLX5","DLX6","POU3F2","POU3F3"}
NEFTEL_MES = {"CD44","CHI3L1","VIM","SERPINA3","ANXA1","TGFBI","CSTB","S100A11"}
NEFTEL_AC = {"GFAP","S100B","APOE","AQP4","ALDH1L1","SLC1A3","MLC1","FGFR3","SPARCL1"}
NEFTEL_OPC = {"OLIG1","OLIG2","PDGFRA","NFIA","APOD","CSPG4","SOX10"}
NEFTEL_NPC = {"DLL3","DLX1","DLX2","DLX5","DLX6","SOX11","NEUROD4","STMN2","TUBB3"}
TME_MYELOID = {"TMEM119","P2RY12","CX3CR1","CD68","ITGAM","AIF1","CSF1R","TREM2","IBA1"}
CHROMATIN = {"EZH2","SUZ12","EED","BMI1","CBX2","CBX4","CBX7","CBX8","KDM6A","KDM6B","KDM5A","KDM5C","SETD2","KMT2A","KMT2D","DNMT1","DNMT3A","DNMT3B","TET1","TET2","TET3"}

HOUSEKEEPING = {"ACTB","GAPDH","TUBB","B2M","RPLP0","RPL13A","HPRT1","PPIA","UBC",
                "VCP","MCM3","MCM4","MCM5","CCT1","CCT2","CCT3","CCT4","CCT5","CCT6A","CCT7","CCT8",
                "HSPA5","HSPA8","CALR","PDI","PPIB","RAN","RAB1A",
                "HMGCS1","SREBF1","FASN",  # metabolism
                "POLR2A","POLR2B","POLR2C",
                "RPS6","RPL7","EIF4A1","EIF4E"}

def gbm_tags(gene: str) -> list[str]:
    tags = []
    if gene in GLIOMA_DRIVERS: tags.append("glioma driver")
    if gene in GSC_LINEAGE: tags.append("GSC/lineage TF")
    if gene in NEFTEL_MES: tags.append("Neftel MES-like")
    if gene in NEFTEL_AC: tags.append("Neftel AC-like")
    if gene in NEFTEL_OPC: tags.append("Neftel OPC-like")
    if gene in NEFTEL_NPC: tags.append("Neftel NPC-like")
    if gene in TME_MYELOID: tags.append("⚠️ TME myeloid")
    if gene in CHROMATIN: tags.append("chromatin remodeler")
    if gene in HOUSEKEEPING: tags.append("⚠️ housekeeping")
    return tags


# --- Load scores -------------------------------------------------------------
print(f"[v2] loading pan-malignant scores: {PAN_SCORES}")
pan = pl.read_parquet(PAN_SCORES).to_pandas()
print(f"[v2] n peaks: {len(pan):,}")

print(f"[v2] loading all-cell-type scores: {CELLTYPE_SCORES}")
ct = pl.read_parquet(CELLTYPE_SCORES).to_pandas()
# Compute per-peak MAX NON-MALIGNANT strength for selectivity
ct_nonmal = ct[~ct["cell_type"].isin(["malignant_unresolved", "unassigned"])]
nonmal_strength = (
    ct_nonmal.groupby("peak_id")["strength"].max()
    .rename("strength_max_nonmalignant")
    .reset_index()
)
pan_enriched = pan.merge(nonmal_strength, on="peak_id", how="left")
pan_enriched["strength_max_nonmalignant"] = pan_enriched["strength_max_nonmalignant"].fillna(0.0)
pan_enriched["selectivity_vs_nonmal"] = pan_enriched["strength_mean"] / pan_enriched["strength_max_nonmalignant"].clip(lower=0.01)
print(f"[v2] joined non-malignant strengths; median selectivity = {pan_enriched['selectivity_vs_nonmal'].median():.2f}")


# --- Nearest-gene annotation (binary-search over TSSes) ----------------------
gene_df = pd.read_csv(GENE_BED, sep="\t", header=None,
                     names=["Chromosome","Start","End","Name","GeneType","Strand"])
gene_df["tss"] = gene_df.apply(
    lambda r: r["Start"] if r["Strand"] == "+" else r["End"], axis=1
).astype(int)
tss_by_chr = {}
for ch, grp in gene_df.groupby("Chromosome"):
    g = grp.sort_values("tss").reset_index(drop=True)
    tss_by_chr[ch] = (g["tss"].values, g["Name"].values, g["Strand"].values)

def nearest_gene(chrom, peak_center):
    if chrom not in tss_by_chr: return (None, None, None)
    tss_arr, names, strands = tss_by_chr[chrom]
    idx = bisect.bisect_left(tss_arr, peak_center)
    cands = []
    if idx > 0:            cands.append(idx - 1)
    if idx < len(tss_arr): cands.append(idx)
    best = min(cands, key=lambda i: abs(peak_center - tss_arr[i]))
    return (names[best], strands[best], int(peak_center - tss_arr[best]))


# --- Build annotated top table -----------------------------------------------
# First: candidate pool = n_cohorts >= 4 (stricter than before)
pool = pan_enriched[pan_enriched["n_cohorts"] >= 4].copy()
print(f"[v2] n_cohorts >= 4 pool: {len(pool):,}")

# Annotate every candidate with nearest gene
annots = []
for _, r in pool.iterrows():
    center = (int(r["start"]) + int(r["end"])) // 2
    name, strand, dist = nearest_gene(r["chrom"], center)
    annots.append((name, strand, dist))
pool[["nearest_gene","gene_strand","dist_to_tss_signed"]] = pd.DataFrame(annots, index=pool.index)
pool["dist_to_tss_abs_kb"] = (pool["dist_to_tss_signed"].abs() / 1000).round(1)
pool["gbm_tags"] = pool["nearest_gene"].fillna("").map(lambda g: ", ".join(gbm_tags(g)))

# Clean/scored for AAV cloning: distal (>=2 kb from TSS) + selective (>=2x non-malignant)
# + non-chr7 (CNV-safe) + housekeeping-free
HOUSEKEEPING_PAT = pool["gbm_tags"].str.contains("housekeeping", na=False)
TME_PAT = pool["gbm_tags"].str.contains("TME myeloid", na=False)
clean = pool[
    (pool["dist_to_tss_abs_kb"] >= 2.0)
    & (pool["dist_to_tss_abs_kb"] <= 100.0)
    & (pool["selectivity_vs_nonmal"] >= 2.0)
    & (pool["chrom"] != "chr7")
    & (~HOUSEKEEPING_PAT)
].copy()
print(f"[v2] AAV-clean pool (distal + selective ≥2x + non-chr7 + non-housekeeping): {len(clean):,}")

# Rank the clean pool by a plain, interpretable composite
clean["aav_score"] = (
    clean["consistency"] * 0.35
    + clean["n_cohorts"].astype(float) / 8.0 * 0.25
    + (clean["n_patients_detected"].astype(float) / 33.0).clip(upper=1) * 0.15
    + (clean["selectivity_vs_nonmal"].clip(upper=10) / 10) * 0.15
    + clean["strength_mean"].clip(upper=0.5) / 0.5 * 0.10
)
clean = clean.sort_values("aav_score", ascending=False).reset_index(drop=True)

# Full top 100 (promoters + non-selective too, for reference)
all100 = pool.sort_values(
    ["n_cohorts","n_patients_detected","consistency","selectivity_vs_nonmal"],
    ascending=[False, False, False, False],
).head(100)

# Output columns
cols_clean = [
    "peak_id","chrom","start","end",
    "nearest_gene","dist_to_tss_signed","gbm_tags",
    "n_cohorts","n_patients_detected","n_patients_accessible_10pct",
    "consistency","strength_mean","strength_max_nonmalignant","selectivity_vs_nonmal",
    "aav_score",
]
clean[cols_clean].head(50).to_csv(OUT_CSV_CLEAN, index=False)
print(f"[v2] wrote {OUT_CSV_CLEAN}")

all100[["peak_id","chrom","start","end","nearest_gene","dist_to_tss_signed","gbm_tags",
        "n_cohorts","n_patients_detected","n_patients_accessible_10pct",
        "consistency","strength_mean","strength_max_nonmalignant","selectivity_vs_nonmal"]
      ].to_csv(OUT_CSV_ALL, index=False)
print(f"[v2] wrote {OUT_CSV_ALL}")

# BED
with open(OUT_BED, "w") as f:
    for _, r in clean.head(30).iterrows():
        name = f"{r['nearest_gene']}_{int(r['dist_to_tss_signed']/1000):+d}kb_{r['peak_id']}" if pd.notna(r["nearest_gene"]) else r["peak_id"]
        f.write(f"{r['chrom']}\t{int(r['start'])}\t{int(r['end'])}\t{name}\t{int(r['aav_score']*1000)}\t.\n")
print(f"[v2] wrote {OUT_BED}")


# --- HTML report --------------------------------------------------------------
def esc(s): return html.escape(str(s)) if pd.notna(s) else ""

def rec_blurb(r):
    tags = r["gbm_tags"]
    bits = []
    if pd.notna(r["nearest_gene"]):
        d_kb = r["dist_to_tss_signed"] / 1000
        if abs(d_kb) < 2:
            bits.append(f"promoter of <strong>{r['nearest_gene']}</strong>")
        elif abs(d_kb) < 100:
            bits.append(f"<strong>{int(d_kb):+d} kb</strong> from <strong>{r['nearest_gene']}</strong> TSS")
        else:
            bits.append(f"intergenic (<strong>{int(d_kb):+d} kb</strong> from {r['nearest_gene']})")
    sel = r["selectivity_vs_nonmal"]
    if sel >= 5:
        bits.append(f"<strong style='color:#15803d'>{sel:.1f}× more open in malignant than any non-malignant cell type</strong>")
    elif sel >= 2:
        bits.append(f"<span style='color:#059669'>{sel:.1f}× more open in malignant vs non-malignant</span>")
    bits.append(f"replicated in {int(r['n_cohorts'])}/8 cohorts, {int(r['n_patients_accessible_10pct'])}/{int(r['n_patients_detected'])} patients ≥10% open")
    if tags:
        bits.append(f"<em style='color:#059669'>gene context: {esc(tags)}</em>")
    return " · ".join(bits)

def table_rows(df):
    rows = []
    for _, r in df.iterrows():
        d_kb = r["dist_to_tss_signed"] / 1000 if pd.notna(r["dist_to_tss_signed"]) else None
        d_str = f"{d_kb:+.1f} kb" if d_kb is not None else "—"
        tag_html = f'<span style="color:#059669;font-weight:600;">{esc(r["gbm_tags"])}</span>' if r["gbm_tags"] else ""
        chr7_flag = '<span style="color:#dc2626;" title="chr7 — CNV-driven signal possible">⚠️</span>' if r["chrom"] == "chr7" else ""
        rows.append(f"""
          <tr>
            <td><code>{esc(r['peak_id'])}</code> {chr7_flag}</td>
            <td><strong>{esc(r['nearest_gene'])}</strong><br/>{tag_html}</td>
            <td style="text-align:right;">{d_str}</td>
            <td style="text-align:center;">{int(r['n_cohorts'])}/8</td>
            <td style="text-align:center;">{int(r['n_patients_detected'])}</td>
            <td style="text-align:center;">{int(r['n_patients_accessible_10pct'])}</td>
            <td style="text-align:center;">{r['consistency']*100:.0f}%</td>
            <td style="text-align:center;">{r['strength_mean']*100:.1f}%</td>
            <td style="text-align:center;">{r['strength_max_nonmalignant']*100:.1f}%</td>
            <td style="text-align:center;"><strong>{r['selectivity_vs_nonmal']:.2f}×</strong></td>
          </tr>
        """)
    return "\n".join(rows)


html_doc = f"""<!DOCTYPE html>
<html><head>
<title>Pan-malignant AAV enhancer candidates — v2</title>
<meta charset="utf-8"/>
<style>
  body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; max-width: 1500px; margin: 2em auto; padding: 0 2em; color: #111; line-height: 1.5; }}
  h1 {{ font-size: 1.8em; border-bottom: 2px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.4em; margin-top: 2em; color: #1e40af; border-bottom: 1px solid #e5e7eb; padding-bottom: 0.2em; }}
  h3 {{ font-size: 1.1em; color: #334155; }}
  code {{ font-size: 0.85em; background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85em; margin: 1em 0; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 5px 8px; vertical-align: top; }}
  th {{ background: #f9fafb; text-align: left; font-weight: 600; }}
  tr:nth-child(even) {{ background: #fafbfc; }}
  tr:hover {{ background: #eff6ff; }}
  .note {{ background: #fef3c7; border-left: 4px solid #f59e0b; padding: 1em; margin: 1em 0; }}
  .success {{ background: #d1fae5; border-left: 4px solid #059669; padding: 1em; margin: 1em 0; }}
  .pick {{ background: #e0f2fe; border-left: 4px solid #0284c7; padding: 0.75em; margin: 0.5em 0; }}
  ol li, ul li {{ margin: 0.3em 0; }}
</style></head><body>

<h1>Pan-malignant AAV enhancer candidates — GBM atlas v2 (2026-10-07)</h1>

<p><strong>Source atlas:</strong> 8 cohorts, 51 patients, 1,523,668 nuclei. Pan-malignant pool
built from <code>malignant_cnv == 1 OR cell_type == malignant_unresolved</code>, which
recovered 438,994 CNV-malignant cells in GBM-Space that marker-peak scoring had
mis-assigned to TME types.</p>

<h2>What changed from v1</h2>
<ol>
<li><strong>Fixed nearest-gene bug</strong> — pyranges was returning 29 Mb-distant pseudogenes
when protein-coding genes were within kb. Rewrote with manual TSS binary search.</li>
<li><strong>Added selectivity-vs-non-malignant</strong> — now reports how many × more open
each peak is in malignant cells vs. the best-matched non-malignant cell type from the
main cell-type matrix. Critical for AAV cell-type selectivity.</li>
<li><strong>Added distal + housekeeping filters</strong> — the AAV cloning shortlist now
excludes promoters (|d|<2 kb from TSS, which drive broad expression) and housekeeping
genes (MCM3, VCP, CCT8 etc., open in every proliferating cell).</li>
</ol>

<div class="note">
<strong>Important caveat.</strong> The v1 top-15 was dominated by promoter-proximal
peaks (PRKCSH, FEM1A, UBAP2, MCM3, VCP, STAT2…). These are <strong>housekeeping
promoters</strong>, open in every proliferating cell — exactly what you DON'T want
for AAV cell-type-selective targeting. v2 filters them out so the shortlist actually
reflects distal, malignant-selective enhancer candidates.
</div>

<div class="success">
<strong>Headline numbers (v2):</strong>
<ul>
<li>Candidate pool (n_cohorts ≥ 4): <strong>{len(pool):,}</strong> peaks</li>
<li>AAV-clean pool (distal 2-100 kb + selective ≥2× + non-chr7 + non-housekeeping): <strong>{len(clean):,}</strong> peaks</li>
<li>Highest aav_score in the clean pool: <strong>{clean['aav_score'].iloc[0]:.3f}</strong></li>
<li>Max cross-cohort replication achieved: <strong>{int(clean['n_cohorts'].max())}/8 cohorts</strong></li>
</ul>
</div>

<h2>Top 10 AAV-cloning candidates (distal, selective, replicated)</h2>

<p>Reading each entry: a short English blurb, then the full numeric row below.</p>
"""

for i, (_, r) in enumerate(clean.head(10).iterrows(), start=1):
    html_doc += f'<div class="pick"><strong>#{i}: <code>{esc(r["peak_id"])}</code></strong><br/>{rec_blurb(r)}</div>'

html_doc += f"""

<h2>Top 50 AAV-cloning candidates — full table</h2>

<p>Sorted by <code>aav_score</code> (consistency × 0.35 + n_cohorts/8 × 0.25 + n_patients/33 × 0.15 + selectivity/10 × 0.15 + strength × 0.10).
Non-chr7, distal (2-100 kb from TSS), ≥2× more open in malignant than any non-malignant cell type.</p>

<table>
<thead><tr>
  <th>peak_id (GRCh38)</th>
  <th>nearest gene<br/>GBM tags</th>
  <th>dist to TSS</th>
  <th>cohorts</th>
  <th>n_pt det.</th>
  <th>n_pt ≥10%</th>
  <th>consist.</th>
  <th>strength (mal)</th>
  <th>strength (max non-mal)</th>
  <th>selectivity</th>
</tr></thead>
<tbody>
{table_rows(clean.head(50))}
</tbody></table>

<h2>Reference: top 100 pan-malignant (unfiltered)</h2>

<p>This includes promoters, chr7 hits, housekeeping genes etc. Shown for reference —
<strong>don't clone from this table</strong>, clone from the filtered one above.</p>

<table>
<thead><tr>
  <th>peak_id</th>
  <th>nearest gene<br/>GBM tags</th>
  <th>dist TSS</th>
  <th>cohorts</th>
  <th>n_pt det.</th>
  <th>n_pt ≥10%</th>
  <th>consist.</th>
  <th>str (mal)</th>
  <th>str (max non-mal)</th>
  <th>sel</th>
</tr></thead>
<tbody>
{table_rows(all100)}
</tbody></table>

<h2>Downstream files</h2>
<ul>
<li><code>pan_malignant_top50_distal_selective.csv</code> — the AAV cloning shortlist, Excel-ready.</li>
<li><code>pan_malignant_top100_all.csv</code> — top 100 pan-malignant for reference (includes promoters etc.).</li>
<li><code>pan_malignant_top30_distal_selective.bed</code> — GRCh38 BED, 30 best candidates, Benchling / UCSC browser upload.</li>
<li><code>pan_malignant_scores.parquet</code> — all {len(pan):,} scored peaks; query locally with polars for any custom filter.</li>
</ul>

<h2>Caveats / limits of this analysis</h2>
<ul>
<li><strong>chr7 excluded on principle.</strong> chr7+ gain is universal in GBM, so chr7 peaks have inflated strength from CNV, not cell-type selectivity. Even if a chr7 peak is in a biologically perfect spot, re-ranking against a CNV-matched reference is required before trusting it.</li>
<li><strong>Selectivity reference is coarse.</strong> We compare malignant strength to the <em>max</em> strength across non-malignant cell types (TAM, neuron, OPC, microglia, etc.). A peak that's <em>very</em> open in e.g. neurons will score poorly on selectivity even if it's clearly tumor-driving in a non-neuron context. Future iteration: compare malignant vs each specific non-malignant type.</li>
<li><strong>Promoter cutoff is 2 kb.</strong> Some bona fide enhancers sit within 2 kb of a TSS. If a candidate looks great on all other metrics but gets excluded by the distal filter, pull it from <code>pan_malignant_top100_all.csv</code> and inspect manually.</li>
<li><strong>Marker-peak mis-labeling residual.</strong> Even with the pan-malignant rebuild, cohorts where CNV wasn't computed (tcga_scatac, gse276177_khan_astro) contribute less than they should. 2/8 cohorts max achievable there.</li>
</ul>

</body></html>
"""

OUT_HTML.write_text(html_doc)
print(f"[v2] wrote {OUT_HTML}")
print("\nDONE.")
