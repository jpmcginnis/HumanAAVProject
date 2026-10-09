"""
Pan-malignant enhancer candidate report v3.

What's new in v3 vs v2:
  - Uses 45-patient scores (v4 pan_malignant run on 2026-10-07) — tcga_scatac
    and gse276177_khan_astro are now included (CNV fix in flatten_cohort.py)
  - Cohort QC caveats baked in: hra004942 low-cell platform, gse165037 sci-ATAC,
    gse276177 is actually Sojka et al. (not Khan), guilhamon hg38-verified
  - n_patients denominator updated 33 → 45

Produces four files in reports/:
  pan_malignant_report_v3.html
  pan_malignant_top50_distal_selective_v3.csv   -- the AAV cloning shortlist
  pan_malignant_top100_all_v3.csv               -- top 100 including promoters/non-selective (for ref)
  pan_malignant_top30_distal_selective_v3.bed   -- GRCh38 BED, Benchling-ready
"""
from __future__ import annotations
import bisect, html, sys
from pathlib import Path

import pandas as pd
import polars as pl

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results")
TOP_ROOT = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26")

# v4 scores are on local disk (freshly scp'd from EC2); fall back to Drive if needed.
PAN_SCORES = HERE / "data/matrix/pan_malignant_scores_v4.parquet"
CELLTYPE_SCORES = TOP_ROOT / "2026-10-6 dataset and results/data/matrix/candidate_scores.parquet"
GENE_BED = HERE / "reference/refgene_hg38.bed"
OUT_HTML = HERE / "reports/pan_malignant_report_v3.html"
OUT_CSV_CLEAN = HERE / "reports/pan_malignant_top50_distal_selective_v3.csv"
OUT_CSV_ALL = HERE / "reports/pan_malignant_top100_all_v3.csv"
OUT_BED = HERE / "reports/pan_malignant_top30_distal_selective_v3.bed"
OUT_HTML.parent.mkdir(parents=True, exist_ok=True)

# Updated patient totals from v4 run (2026-10-07)
N_PATIENTS_TOTAL = 51
N_PATIENTS_IN_POOL = 45  # patients contributing ≥10 cells to pan-malignant pool
N_MALIGNANT_CELLS = 575_971  # sum across all 8 cohorts in v4 log
N_COHORTS = 8


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

# Daigle-style Z-score (Allen Institute Armamentarium vocabulary):
# Z = (target_strength - mean_other_strengths) / sd_other_strengths
# Lets us report "Z > 2 by Daigle criterion" in the R01, the field-standard
# threshold for cell-type-specific peak nomination.
# "Other" = non-malignant, non-unassigned cell types (as above).
nonmal_stats = (
    ct_nonmal.groupby("peak_id")["strength"]
    .agg(nonmal_mean="mean", nonmal_sd="std")
    .reset_index()
)
pan_enriched = pan_enriched.merge(nonmal_stats, on="peak_id", how="left")
pan_enriched["nonmal_mean"] = pan_enriched["nonmal_mean"].fillna(0.0)
# Floor sd at a small epsilon so peaks absent from most cell types don't get infinite Z
pan_enriched["z_daigle"] = (
    (pan_enriched["strength_mean"] - pan_enriched["nonmal_mean"])
    / pan_enriched["nonmal_sd"].clip(lower=0.001)
)
pan_enriched["passes_daigle_z2"] = pan_enriched["z_daigle"] >= 2.0
print(f"[v2] Daigle Z >= 2 peaks (atlas-wide): {int(pan_enriched['passes_daigle_z2'].sum()):,}")
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
# Annotation (not a filter): where this peak sits relative to the nearest gene.
# For AAV cloning the enhancer sequence is extracted from genomic context, so
# this is context for the reader, not a reason to drop any row.
def _dist_cat(kb):
    if pd.isna(kb): return "unknown"
    if kb < 2:    return "promoter"       # excluded by the 2 kb floor
    if kb < 10:   return "near"            # 2-10 kb
    if kb < 100:  return "distal"          # classic 10-100 kb enhancer range
    if kb < 500:  return "far-distal"      # common for cell-type-specific enhancers
    return "gene-desert"                   # >500 kb — often enhancer-dense regions
pool["distance_category"] = pool["dist_to_tss_abs_kb"].map(_dist_cat)
pool["gbm_tags"] = pool["nearest_gene"].fillna("").map(lambda g: ", ".join(gbm_tags(g)))

# Clean/scored for AAV cloning: distal (>=2 kb from TSS) + selective (>=2x non-malignant)
# + non-chr7 (CNV-safe) + housekeeping-free
HOUSEKEEPING_PAT = pool["gbm_tags"].str.contains("housekeeping", na=False)
TME_PAT = pool["gbm_tags"].str.contains("TME myeloid", na=False)
clean = pool[
    (pool["dist_to_tss_abs_kb"] >= 2.0)
    # NOTE 2026-10-08: dropped the 100 kb upper cap. For AAV, the enhancer is
    # extracted from its native context and placed next to a minimal promoter,
    # so native-genome distance to the "regulated" gene is irrelevant. TF binding
    # is a sequence property, not a location property. Keeping the 2 kb floor to
    # exclude promoter-proximal peaks (not enhancers in the regulatory sense).
    # Armamentarium precedent: Mich et al. + Hooks striatum papers both used
    # ~500 kb search windows around marker-gene loci; several Hunker hits sit
    # 200-400 kb from their marker genes, AiE0387m sits in a gene desert.
    # See README.md "distance filter" note.
    & (pool["selectivity_vs_nonmal"] >= 2.0)
    & (pool["chrom"] != "chr7")
    & (~HOUSEKEEPING_PAT)
].copy()
print(f"[v2] AAV-clean pool (distal + selective ≥2x + non-chr7 + non-housekeeping): {len(clean):,}")

# Rank the clean pool by a plain, interpretable composite
clean["aav_score"] = (
    clean["consistency"] * 0.35
    + clean["n_cohorts"].astype(float) / 8.0 * 0.25
    + (clean["n_patients_detected"].astype(float) / float(N_PATIENTS_IN_POOL)).clip(upper=1) * 0.15
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
    "nearest_gene","dist_to_tss_signed","distance_category","gbm_tags",
    "n_cohorts","n_patients_detected","n_patients_accessible",
    "consistency","strength_mean","strength_max_nonmalignant","selectivity_vs_nonmal",
    "nonmal_mean","nonmal_sd","z_daigle","passes_daigle_z2",
    "aav_score",
]
clean[cols_clean].head(50).to_csv(OUT_CSV_CLEAN, index=False)
print(f"[v2] wrote {OUT_CSV_CLEAN}")

all100[["peak_id","chrom","start","end","nearest_gene","dist_to_tss_signed","distance_category","gbm_tags","z_daigle","passes_daigle_z2",
        "n_cohorts","n_patients_detected","n_patients_accessible",
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
    bits.append(f"replicated in {int(r['n_cohorts'])}/8 cohorts, {int(r['n_patients_accessible'])}/{int(r['n_patients_detected'])} patients ≥10% open")
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
            <td style="text-align:center;">{int(r['n_patients_accessible'])}</td>
            <td style="text-align:center;">{r['consistency']*100:.0f}%</td>
            <td style="text-align:center;">{r['strength_mean']*100:.1f}%</td>
            <td style="text-align:center;">{r['strength_max_nonmalignant']*100:.1f}%</td>
            <td style="text-align:center;"><strong>{r['selectivity_vs_nonmal']:.2f}×</strong></td>
            <td style="text-align:center;">{('<strong style="color:#059669;">' if r['z_daigle'] >= 2 else '')}{r['z_daigle']:.2f}{'</strong>' if r['z_daigle'] >= 2 else ''}</td>
          </tr>
        """)
    return "\n".join(rows)


html_doc = f"""<!DOCTYPE html>
<html><head>
<title>Pan-malignant AAV enhancer candidates — v3</title>
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

<h1>Pan-malignant AAV enhancer candidates — GBM atlas v3 (2026-10-07)</h1>

<p><strong>Source atlas:</strong> {N_COHORTS} cohorts, {N_PATIENTS_TOTAL} patients
(45 contributing ≥10 cells to pan-malignant pool), 1,523,668 nuclei total,
{N_MALIGNANT_CELLS:,} pan-malignant cells after CNV + label join. Pool rule:
<code>malignant_cnv == 1 OR cell_type == malignant_unresolved</code>.</p>

<h2>What changed from v2</h2>
<ol>
<li><strong>Fixed CNV missing for fragments-mode cohorts</strong> — v2 silently defaulted
<code>malignant_cnv=0</code> for tcga_scatac and gse276177 because <code>src/flatten_cohort.py</code>
skipped the CNV calling path. Patched to call <code>_parallel_chr7_10_ratio</code> in-process;
recovered <strong>12 patients</strong> and <strong>127,045 malignant cells</strong> for the pool.
Candidate pool size up accordingly.</li>
<li><strong>Patient denominator updated 33 → 45</strong> for n_patients_detected normalization
in <code>aav_score</code>.</li>
<li><strong>Cohort QC caveats incorporated</strong> — see Methods / Caveats section below.
Short version: hra004942 is low-cell-per-patient (20-250 malignant/pt), gse165037 is sci-ATAC
(not 10x), gse276177 is really Sojka et al. 2025 (slug misnomer), guilhamon confirmed hg38
(2,332 malignant cells across 4 pts match CATLAS peaks).</li>
</ol>

<div class="note">
<strong>Important caveat (preserved from v2).</strong> The v1 top-15 was dominated by
promoter-proximal peaks (PRKCSH, FEM1A, UBAP2, MCM3, VCP, STAT2…). These are
<strong>housekeeping promoters</strong>, open in every proliferating cell — exactly
what you DON'T want for AAV cell-type-selective targeting. v2+ filters them out so
the shortlist reflects distal, malignant-selective enhancer candidates.
</div>

<div class="success">
<strong>Headline numbers (v2):</strong>
<ul>
<li>Candidate pool (n_cohorts ≥ 4): <strong>{len(pool):,}</strong> peaks</li>
<li>AAV-clean pool (|dist to TSS| ≥ 2 kb, no upper cap + selective ≥2× + non-chr7 + non-housekeeping): <strong>{len(clean):,}</strong> peaks</li>
<li><strong>Pass Daigle Z ≥ 2</strong> (Allen Institute Armamentarium criterion for cell-type specificity): <strong>{int(clean['passes_daigle_z2'].sum()):,}</strong> peaks within the AAV-clean pool, <strong>{int(pan_enriched['passes_daigle_z2'].sum()):,}</strong> atlas-wide.</li>
<li>Highest aav_score in the clean pool: <strong>{clean['aav_score'].iloc[0]:.3f}</strong></li>
<li>Max cross-cohort replication achieved: <strong>{int(clean['n_cohorts'].max())}/8 cohorts</strong></li>
</ul>
</div>

<div class="note">
<strong>Daigle Z-score</strong> — Z = (strength<sub>pan-mal</sub> − mean<sub>other</sub>) / sd<sub>other</sub>,
where "other" = non-malignant / non-unassigned cell types (astrocyte, oligo, OPC, neuron,
GABA_neuron, microglia, TAM, endothelial, T_cell). Z ≥ 2 is the standard Armamentarium
threshold for calling a peak "cell-type specific." Columns in the output: <code>z_daigle</code>,
<code>passes_daigle_z2</code>, <code>nonmal_mean</code>, <code>nonmal_sd</code>.
</div>

<h2>Top 10 AAV-cloning candidates (distal, selective, replicated)</h2>

<p>Reading each entry: a short English blurb, then the full numeric row below.</p>
"""

for i, (_, r) in enumerate(clean.head(10).iterrows(), start=1):
    html_doc += f'<div class="pick"><strong>#{i}: <code>{esc(r["peak_id"])}</code></strong><br/>{rec_blurb(r)}</div>'

html_doc += f"""

<h2>Top 50 AAV-cloning candidates — full table</h2>

<p>Sorted by <code>aav_score</code> (consistency × 0.35 + n_cohorts/8 × 0.25 + n_patients/33 × 0.15 + selectivity/10 × 0.15 + strength × 0.10).
Non-chr7, |dist to TSS| ≥ 2 kb (no upper cap — AAV extracts the enhancer from genomic context,
so native distance to the nominally regulated gene is irrelevant), ≥2× more open in malignant than
any non-malignant cell type. See <code>distance_category</code> column for context
(near 2-10 kb / distal 10-100 kb / far-distal 100-500 kb / gene-desert >500 kb).</p>

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
  <th>Z<sub>Daigle</sub></th>
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
  <th>Z<sub>Daigle</sub></th>
</tr></thead>
<tbody>
{table_rows(all100)}
</tbody></table>

<h2>Downstream files</h2>
<ul>
<li><code>pan_malignant_top50_distal_selective_v3.csv</code> — the AAV cloning shortlist, Excel-ready.</li>
<li><code>pan_malignant_top100_all_v3.csv</code> — top 100 pan-malignant for reference (includes promoters etc.).</li>
<li><code>pan_malignant_top30_distal_selective_v3.bed</code> — GRCh38 BED, 30 best candidates, Benchling / UCSC browser upload.</li>
<li><code>data/matrix/pan_malignant_scores_v4.parquet</code> — all {len(pan):,} scored peaks; query locally with polars for any custom filter.</li>
<li><code>data/matrix/pan_malignant_matrix_v4.parquet</code> — 2.9 M-row (peak × patient × cohort) long-form matrix for custom aggregations.</li>
</ul>

<h2>Caveats / limits of this analysis</h2>
<ul>
<li><strong>chr7 excluded on principle.</strong> chr7+ gain is universal in GBM, so chr7 peaks have inflated strength from CNV, not cell-type selectivity. Even if a chr7 peak is in a biologically perfect spot, re-ranking against a CNV-matched reference is required before trusting it.</li>
<li><strong>Selectivity reference is coarse.</strong> We compare malignant strength to the <em>max</em> strength across non-malignant cell types (TAM, neuron, OPC, microglia, etc.). A peak that's <em>very</em> open in e.g. neurons will score poorly on selectivity even if it's clearly tumor-driving in a non-neuron context. Future iteration: compare malignant vs each specific non-malignant type.</li>
<li><strong>Promoter cutoff is 2 kb.</strong> Some bona fide enhancers sit within 2 kb of a TSS. If a candidate looks great on all other metrics but gets excluded by the distal filter, pull it from <code>pan_malignant_top100_all_v3.csv</code> and inspect manually.</li>
</ul>

<h2>Cohort-level methods caveats (new in v3)</h2>
<p>Per the methods comparison report (reports/methods_comparison.md). These impact how
much weight to put on cross-cohort replication for a given peak:</p>
<ul>
<li><strong>hra004942</strong> (gbm_tme_atlas) — only 20-250 malignant cells per patient
(vs 10K-80K in gbm_space). Platform may be low-throughput snATAC or a different chemistry
than 10x. A peak replicated here contributes less information than one replicated in a
10K-cell patient. Peaks relying heavily on this cohort for replication should be viewed cautiously.</li>
<li><strong>gse165037</strong> — sci-ATAC (combinatorial indexing, Sinnamon et al.-style),
not 10x droplet. Peaks per cell much shallower; co-opening patterns can diverge from
droplet. Treat as orthogonal technology; agreement here is extra-strong evidence.</li>
<li><strong>gse276177 (slug "khan_astro")</strong> — actually Sojka et al. 2025 Nature Cell Biol.
Rename TODO; data content and analysis are correct. Three IDH-WT GBM patients, 10x snATAC.</li>
<li><strong>guilhamon</strong> — pan-malignant pool matches CATLAS peaks (2,332 cells, 10-25K
accessible peaks per patient) which confirms hg38 build. 4 patients, 500-700 malignant cells each.</li>
<li><strong>tcga_scatac (sundaram_gbm)</strong> — recovered in v3 after CNV fix. 9 patients,
53,632 malignant cells. CATLAS peak coverage is 540K+ per patient (near-saturated), unusually
deep.</li>
</ul>

</body></html>
"""

OUT_HTML.write_text(html_doc)
print(f"[v2] wrote {OUT_HTML}")
print("\nDONE.")
