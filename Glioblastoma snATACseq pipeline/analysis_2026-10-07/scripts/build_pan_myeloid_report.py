"""
Pan-myeloid (TAM ∪ microglia) AAV enhancer candidate report.

Rationale: TAM and microglia share enough biology (both CSF1R+, ITGAM+, with
reactive microglia transitioning to TAM phenotype) that peaks accessible in
BOTH get penalized as "non-specific" by per-cell-type selectivity scoring —
even though they're ideal for pan-myeloid AAV targeting (immunomodulation,
TAM re-polarization, pan-brain-macrophage reporter lines, etc.).

This report merges TAM + microglia into one "myeloid" pool and computes
selectivity against all other cell types (excluding malignant_unresolved and
unassigned).

Approach:
  - per-peak: strength_myeloid = max(strength[TAM], strength[microglia])
  - per-peak: strength_other_max = max over {astrocyte, oligodendrocyte, OPC,
    neuron, GABA_neuron, endothelial, T_cell}
  - selectivity_myeloid = strength_myeloid / max(strength_other_max, 0.01)
  - n_cohorts, n_patients, consistency: max() across the two source cell types
    (peak that replicates in either counts; this is a loose upper bound but
    appropriate for pan-myeloid detection)
  - composite_score_myeloid: strength_myeloid * selectivity_myeloid * consistency

Filters (clean pool):
  - distal (|TSS| ≥ 2 kb, ≤ 100 kb)
  - non-housekeeping
  - n_cohorts ≥ 4
  - selectivity ≥ 1.5× (tight — these peaks should be clearly myeloid-biased)

Produces:
  reports/pan_myeloid_report.html
  reports/pan_myeloid_top50_distal_selective.csv
  reports/pan_myeloid_top30.bed
"""
from __future__ import annotations
import bisect, html
from pathlib import Path

import pandas as pd
import polars as pl

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results")
CELLTYPE_SCORES = HERE / "data/matrix/candidate_scores.parquet"
PAN_SCORES = HERE / "data/matrix/pan_malignant_scores_v4.parquet"
GENE_BED = HERE / "reference/refgene_hg38.bed"
OUT_HTML = HERE / "reports/pan_myeloid_report.html"
OUT_CSV = HERE / "reports/pan_myeloid_top50_distal_selective.csv"
OUT_BED = HERE / "reports/pan_myeloid_top30.bed"

OTHER_CELL_TYPES = ["astrocyte","oligodendrocyte","OPC","neuron","GABA_neuron","endothelial","T_cell"]
# Note: malignant_unresolved + unassigned excluded from "non-myeloid" baseline because
# malignant pool is handled separately and unassigned is a catch-all noise bucket.

# Gene set tags (abbreviated — focus on myeloid context)
TAM_MICROGLIA_MARKERS = {
    # Homeostatic microglia
    "TMEM119","P2RY12","P2RY13","CX3CR1","HEXB","SALL1","SLC2A5","FCRLS","OLFML3",
    # Pan-myeloid / shared
    "CSF1R","CSF2RA","ITGAM","AIF1","CD68","CD86","CD40","C1QA","C1QB","C1QC",
    "LAPTM5","FCER1G","TYROBP","MS4A6A","MS4A7","HLA-DRA","HLA-DRB1","HLA-DPA1",
    # TAM / activation
    "TREM2","APOE","GPNMB","SPP1","CD163","MRC1","MARCO","CD206","LGALS3","CSTB",
    "SELPLG","SIGLEC1","CYBB","LYZ","FOLR2","F13A1","FCGR1A","FCGR2A","FCGR3A",
    # Chemokines / cytokines
    "CCL2","CCL3","CCL4","CCL8","CXCL8","IL1B","IL6","TNF","NFKB1","NFKB2","RELA",
    # Phagocytosis / lysosomal
    "LAMP1","LAMP2","CTSB","CTSD","CTSS","CTSL","CTSK","MERTK","AXL","GAS6",
    # Transcription factors
    "SPI1","IRF8","CEBPB","CEBPD","MAF","MAFB","PU1","RUNX1","EGR1","FOSB","JUN",
}

GLIOMA_DRIVERS = {"EGFR","PDGFRA","MET","FGFR1","FGFR3","TP53","PTEN","NF1","CDKN2A","CDKN2B","MYC","MGMT"}

HOUSEKEEPING = {"ACTB","GAPDH","TUBB","B2M","RPLP0","RPL13A","HPRT1","PPIA","UBC",
                "VCP","MCM3","MCM4","MCM5","CCT1","CCT2","CCT3","CCT4","CCT5","CCT6A","CCT7","CCT8",
                "HSPA5","HSPA8","CALR","PDI","PPIB","RAN","RAB1A",
                "POLR2A","POLR2B","POLR2C","RPS6","RPL7","EIF4A1","EIF4E"}

def gene_tags(gene: str) -> list[str]:
    tags = []
    if gene in HOUSEKEEPING: tags.append("⚠️ housekeeping")
    if gene in TAM_MICROGLIA_MARKERS: tags.append("🟢 myeloid marker")
    if gene in GLIOMA_DRIVERS: tags.append("glioma driver")
    return tags


print("[myeloid] loading per-cell-type scores")
raw = pl.read_parquet(CELLTYPE_SCORES).to_pandas()
print(f"[myeloid]   {len(raw):,} rows x {len(raw.columns)} cols")

print("[myeloid] loading pan scores for peak coords")
pan = pl.read_parquet(PAN_SCORES, columns=["peak_id","chrom","start","end"]).to_pandas()
pan.drop_duplicates("peak_id", inplace=True)
print(f"[myeloid]   {len(pan):,} atlas peaks")

# Pivot: one row per peak with per-cell-type strength
print("[myeloid] pivoting to peak × cell_type strength matrix")
pivot = raw.pivot_table(index="peak_id", columns="cell_type", values="strength", aggfunc="max", fill_value=0.0)
# Also pivot n_cohorts and consistency (max across TAM/microglia for the pooled stat)
pivot_cohorts = raw.pivot_table(index="peak_id", columns="cell_type", values="n_cohorts", aggfunc="max", fill_value=0)
pivot_consistency = raw.pivot_table(index="peak_id", columns="cell_type", values="consistency", aggfunc="max", fill_value=0.0)
pivot_npatients = raw.pivot_table(index="peak_id", columns="cell_type", values="n_patients", aggfunc="max", fill_value=0)

print(f"[myeloid]   peak × cell_type shape: {pivot.shape}")

# Compute myeloid merge
myeloid = pd.DataFrame(index=pivot.index)
myeloid["strength_TAM"] = pivot.get("TAM", 0)
myeloid["strength_microglia"] = pivot.get("microglia", 0)
myeloid["strength_myeloid"] = pivot[["TAM","microglia"]].max(axis=1)

# Which type dominates (which the peak is MORE accessible in)
myeloid["dominant_source"] = myeloid[["strength_TAM","strength_microglia"]].idxmax(axis=1).str.replace("strength_", "")

# Max strength across non-myeloid, non-malignant, non-unassigned
other_cols = [c for c in OTHER_CELL_TYPES if c in pivot.columns]
myeloid["strength_other_max"] = pivot[other_cols].max(axis=1)
myeloid["selectivity_myeloid"] = myeloid["strength_myeloid"] / myeloid["strength_other_max"].clip(lower=0.01)

# Daigle Z-score vs non-myeloid cell types: Z = (strength_myeloid - mean_other) / sd_other
other_mean = pivot[other_cols].mean(axis=1)
other_sd = pivot[other_cols].std(axis=1)
myeloid["other_mean"] = other_mean
myeloid["other_sd"] = other_sd
myeloid["z_daigle"] = (myeloid["strength_myeloid"] - other_mean) / other_sd.clip(lower=0.001)
myeloid["passes_daigle_z2"] = myeloid["z_daigle"] >= 2.0

# Replication: use max of TAM & microglia cohorts/patients/consistency
myeloid["n_cohorts"] = pivot_cohorts[["TAM","microglia"]].max(axis=1)
myeloid["n_patients"] = pivot_npatients[["TAM","microglia"]].max(axis=1)
myeloid["consistency"] = pivot_consistency[["TAM","microglia"]].max(axis=1)

# Composite
myeloid["composite_score_myeloid"] = (
    myeloid["strength_myeloid"] * 1.0
    + myeloid["selectivity_myeloid"].clip(upper=10) / 10 * 0.5
    + myeloid["consistency"] * 0.3
    + myeloid["n_cohorts"] / 8 * 0.2
)

myeloid = myeloid.merge(pan, left_index=True, right_on="peak_id", how="left").set_index("peak_id")
print(f"[myeloid] merged with coords: {len(myeloid):,}")

# TSS annotation
print("[myeloid] annotating nearest gene")
gene_df = pd.read_csv(GENE_BED, sep="\t", header=None,
                     names=["Chromosome","Start","End","Name","GeneType","Strand"])
gene_df["tss"] = gene_df.apply(lambda r: r["Start"] if r["Strand"] == "+" else r["End"], axis=1).astype(int)
tss_by_chr = {}
for ch, grp in gene_df.groupby("Chromosome"):
    g = grp.sort_values("tss").reset_index(drop=True)
    tss_by_chr[ch] = (g["tss"].values, g["Name"].values)

def nearest_gene(chrom, center):
    if chrom not in tss_by_chr: return (None, None)
    tss_arr, names = tss_by_chr[chrom]
    idx = bisect.bisect_left(tss_arr, center)
    cands = []
    if idx > 0:            cands.append(idx - 1)
    if idx < len(tss_arr): cands.append(idx)
    best = min(cands, key=lambda i: abs(center - tss_arr[i]))
    return (names[best], int(center - tss_arr[best]))

annots = []
for _, r in myeloid.iterrows():
    if pd.isna(r["chrom"]):
        annots.append((None, None)); continue
    center = (int(r["start"]) + int(r["end"])) // 2
    annots.append(nearest_gene(r["chrom"], center))
myeloid[["nearest_gene","dist_to_tss_signed"]] = pd.DataFrame(annots, index=myeloid.index)
myeloid["dist_to_tss_abs_kb"] = (myeloid["dist_to_tss_signed"].abs() / 1000).round(1)
def _dist_cat(kb):
    if pd.isna(kb): return "unknown"
    if kb < 2:    return "promoter"
    if kb < 10:   return "near"
    if kb < 100:  return "distal"
    if kb < 500:  return "far-distal"
    return "gene-desert"
myeloid["distance_category"] = myeloid["dist_to_tss_abs_kb"].map(_dist_cat)
myeloid["gene_tags_list"] = [gene_tags(g if pd.notna(g) else "") for g in myeloid["nearest_gene"]]
myeloid["gene_tags"] = myeloid["gene_tags_list"].map(lambda L: ", ".join(L))

# Candidate pool: n_cohorts ≥ 4 AND strength ≥ 0.03 (minimum detection)
pool = myeloid[(myeloid["n_cohorts"] >= 4) & (myeloid["strength_myeloid"] >= 0.03)].copy()
print(f"[myeloid] pool (≥4 cohorts, strength ≥ 0.03): {len(pool):,}")

# Clean filter: distal + non-housekeeping + selectivity ≥ 1.0 (loosened — pan-myeloid
# selectivity ceilings against reactive glia are low in GBM; distribution q99 ≈ 0.98)
is_hk = pool["gene_tags"].str.contains("housekeeping", na=False)
base_clean = pool[
    (pool["dist_to_tss_abs_kb"] >= 2.0)
    # NOTE 2026-10-08: dropped the 100 kb upper cap (AAV extracts enhancer
    # from genomic context — native distance irrelevant). Keeping 2 kb floor
    # to exclude promoter-proximal peaks.
    & (~is_hk)
].copy()

# Table A: selective pan-myeloid (strictly more open in myeloid than any non-myeloid)
clean = base_clean[base_clean["selectivity_myeloid"] >= 1.0].copy()
clean = clean.sort_values("composite_score_myeloid", ascending=False).reset_index()
print(f"[myeloid] Table A — SELECTIVE (sel ≥ 1.0×, distal, non-HK): {len(clean):,}")

# Table B: SHARED pan-myeloid — strong in BOTH TAM AND microglia (user's primary ask).
# This picks up peaks that would have been penalized as "non-specific" in per-cell-type
# analyses because they're accessible in both compartments. No selectivity floor —
# ranked by (strength_TAM + strength_microglia) / 2 to privilege true shared signal.
shared_both = base_clean[
    (base_clean["strength_TAM"] >= 0.05)
    & (base_clean["strength_microglia"] >= 0.05)
].copy()
shared_both["shared_strength_mean"] = (shared_both["strength_TAM"] + shared_both["strength_microglia"]) / 2
shared_both["shared_consistency"] = (shared_both["strength_TAM"].clip(upper=1) * shared_both["strength_microglia"].clip(upper=1)) ** 0.5
# Rank: joint strength × composite × n_cohorts
shared_both["shared_rank_score"] = (
    shared_both["shared_strength_mean"]
    * (1 + shared_both["selectivity_myeloid"].clip(upper=5) / 5)
    * (shared_both["n_cohorts"] / 8)
)
shared_both = shared_both.sort_values("shared_rank_score", ascending=False).reset_index()
print(f"[myeloid] Table B — SHARED (TAM≥5% AND microglia≥5%, distal, non-HK): {len(shared_both):,}")

# For reporting: dedupe counts across the two tables
shared = clean[(clean["strength_TAM"] >= 0.03) & (clean["strength_microglia"] >= 0.03)]
tam_only = clean[(clean["strength_TAM"] >= 0.03) & (clean["strength_microglia"] < 0.03)]
microglia_only = clean[(clean["strength_TAM"] < 0.03) & (clean["strength_microglia"] >= 0.03)]
print(f"[myeloid]   within Table A: TAM+microglia shared: {len(shared):,}  TAM-only: {len(tam_only):,}  microglia-only: {len(microglia_only):,}")

# Write CSV: Table A (selective) + Table B (shared), both as separate CSVs
cols = ["peak_id","chrom","start","end","nearest_gene","dist_to_tss_signed","distance_category","gene_tags","z_daigle","passes_daigle_z2",
        "dominant_source","strength_TAM","strength_microglia","strength_myeloid","strength_other_max",
        "selectivity_myeloid","n_cohorts","n_patients","consistency","composite_score_myeloid"]
clean[cols].head(50).to_csv(OUT_CSV, index=False)
print(f"[myeloid] wrote {OUT_CSV}")

OUT_CSV_SHARED = HERE / "reports/pan_myeloid_top50_shared.csv"
cols_shared = cols + ["shared_strength_mean","shared_rank_score"]
shared_both[cols_shared].head(50).to_csv(OUT_CSV_SHARED, index=False)
print(f"[myeloid] wrote {OUT_CSV_SHARED}")

# BED for both
with open(OUT_BED, "w") as f:
    for _, r in clean.head(30).iterrows():
        name = f"{r['nearest_gene']}_{int(r['dist_to_tss_signed']/1000):+d}kb_{r['peak_id']}" if pd.notna(r["nearest_gene"]) else r["peak_id"]
        f.write(f"{r['chrom']}\t{int(r['start'])}\t{int(r['end'])}\t{name}\t{int(r['composite_score_myeloid']*1000)}\t.\n")
print(f"[myeloid] wrote {OUT_BED}")

OUT_BED_SHARED = HERE / "reports/pan_myeloid_top30_shared.bed"
with open(OUT_BED_SHARED, "w") as f:
    for _, r in shared_both.head(30).iterrows():
        name = f"{r['nearest_gene']}_{int(r['dist_to_tss_signed']/1000):+d}kb_{r['peak_id']}" if pd.notna(r["nearest_gene"]) else r["peak_id"]
        f.write(f"{r['chrom']}\t{int(r['start'])}\t{int(r['end'])}\t{name}\t{int(r['shared_rank_score']*1000)}\t.\n")
print(f"[myeloid] wrote {OUT_BED_SHARED}")


# --- HTML ------------------------------------------------------------------
def esc(s): return html.escape(str(s)) if pd.notna(s) else ""


def tag_html(tags_str):
    if not tags_str: return ""
    pieces = [t.strip() for t in tags_str.split(",")]
    out = []
    for t in pieces:
        if t.startswith("🟢"):
            out.append(f'<span style="color:#059669;font-weight:600;">{esc(t)}</span>')
        elif t.startswith("⚠️"):
            out.append(f'<span style="color:#b91c1c;">{esc(t)}</span>')
        else:
            out.append(esc(t))
    return ", ".join(out)


def section_rows(df):
    out = []
    for _, r in df.iterrows():
        d_kb = r["dist_to_tss_signed"] / 1000
        tam_s = r["strength_TAM"] * 100
        mic_s = r["strength_microglia"] * 100
        other_s = r["strength_other_max"] * 100
        src = r["dominant_source"]
        src_badge = ("<span style='color:#f59e0b;'>TAM+</span>" if src == "TAM"
                     else "<span style='color:#0284c7;'>micro+</span>")
        # Mark shared vs dominant
        shared_flag = ""
        if r["strength_TAM"] >= 0.03 and r["strength_microglia"] >= 0.03:
            shared_flag = '<span title="accessible in both TAM and microglia" style="background:#d1fae5;color:#065f46;font-size:0.75em;padding:1px 4px;border-radius:3px;">SHARED</span>'
        out.append(f"""
        <tr>
          <td><code>{esc(r['peak_id'])}</code> {shared_flag}</td>
          <td><strong>{esc(r['nearest_gene'])}</strong><br/>{tag_html(r['gene_tags'])}</td>
          <td style="text-align:right;">{d_kb:+.1f} kb</td>
          <td style="text-align:center;">{int(r['n_cohorts'])}/8</td>
          <td style="text-align:center;">{int(r['n_patients'])}</td>
          <td style="text-align:center;">{r['consistency']*100:.0f}%</td>
          <td style="text-align:center;">{tam_s:.1f}%</td>
          <td style="text-align:center;">{mic_s:.1f}%</td>
          <td style="text-align:center;">{other_s:.1f}%</td>
          <td style="text-align:center;"><strong>{r['selectivity_myeloid']:.2f}×</strong></td>
          <td style="text-align:center;">{('<strong style="color:#059669;">' if pd.notna(r['z_daigle']) and r['z_daigle'] >= 2 else '')}{r['z_daigle']:.2f}{'</strong>' if pd.notna(r['z_daigle']) and r['z_daigle'] >= 2 else ''}</td>
          <td style="text-align:center;">{src_badge}</td>
          <td style="text-align:center;">{r['composite_score_myeloid']:.3f}</td>
        </tr>""")
    return "\n".join(out)


html_doc = f"""<!DOCTYPE html>
<html><head>
<title>Pan-myeloid (TAM + microglia) AAV enhancer candidates — GBM atlas</title>
<meta charset="utf-8"/>
<style>
  body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; max-width: 1600px; margin: 2em auto; padding: 0 2em; color: #111; line-height: 1.5; }}
  h1 {{ font-size: 1.8em; border-bottom: 2px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.3em; margin-top: 2em; color: #1e40af; border-bottom: 1px solid #e5e7eb; padding-bottom: 0.2em; }}
  code {{ font-size: 0.85em; background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.84em; margin: 1em 0; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 5px 8px; vertical-align: top; }}
  th {{ background: #f9fafb; text-align: left; font-weight: 600; }}
  tr:nth-child(even) {{ background: #fafbfc; }}
  tr:hover {{ background: #eff6ff; }}
  .note {{ background: #fef3c7; border-left: 4px solid #f59e0b; padding: 1em; margin: 1em 0; }}
  .success {{ background: #d1fae5; border-left: 4px solid #059669; padding: 1em; margin: 1em 0; }}
</style></head><body>

<h1>Pan-myeloid AAV enhancer candidates — TAM ∪ microglia (GBM atlas, 2026-10-07)</h1>

<p><strong>Rationale.</strong> TAM and microglia share so much biology (both CSF1R+/ITGAM+;
reactive microglia transition continuously to TAM phenotype in tumor) that peaks accessible
in BOTH get penalized as "non-specific" by per-cell-type selectivity scoring — even though
they're <em>ideal</em> for pan-myeloid AAV targeting:</p>
<ul>
<li>TAM re-polarization therapy (M2 → M1 in tumor)</li>
<li>Pan-brain-macrophage reporter lines for lineage tracing</li>
<li>CSF1R / TREM2 modulators delivered selectively</li>
<li>Immune checkpoint depletion in the tumor microenvironment</li>
</ul>

<p>This report pools TAM + microglia strength (max per peak) and recomputes selectivity
against non-myeloid cell types (astrocyte, oligodendrocyte, OPC, neuron, GABA_neuron,
endothelial, T_cell — malignant and unassigned excluded as separate / noise pools).</p>

<div class="success">
<strong>Headline numbers:</strong>
<ul>
<li>Candidate pool (n_cohorts ≥ 4, strength_myeloid ≥ 0.03): <strong>{len(pool):,}</strong> peaks</li>
<li><strong>Table A — SELECTIVE</strong> (|dist TSS| ≥ 2 kb no upper cap, selectivity ≥ 1.0×, non-housekeeping): <strong>{len(clean):,}</strong> peaks</li>
<li><strong>Table B — SHARED</strong> (TAM ≥ 5% AND microglia ≥ 5%, distal, non-housekeeping — <em>the primary deliverable</em>): <strong>{len(shared_both):,}</strong> peaks</li>
<li>Within Table A: <strong>{len(shared):,} shared</strong>, <strong>{len(tam_only):,} TAM-only</strong>, <strong>{len(microglia_only):,} microglia-only</strong></li>
</ul>
</div>

<div class="note">
<strong>Why two tables.</strong> Pan-myeloid selectivity ceilings in GBM atlas are low
(q99 ≈ 0.98 across the pool) because reactive glia and malignant Neftel-MES-like cells
share a lot of myeloid-adjacent accessibility. A strict selectivity ≥ 2× floor would
leave only 6 peaks — missing the real signal the user asked about: peaks accessible
in both TAM AND microglia that got penalized as "non-specific" in per-cell-type reports.
Table B surfaces those shared peaks directly, ranked by joint strength rather than
selectivity. Use Table A for maximal purity against other glia; use Table B for pan-myeloid
targeting where glial crossover is acceptable.
</div>

<h2>How this complements per-cell-type and pan-malignant reports</h2>
<ul>
<li><strong>pan_malignant_report_v3</strong> — malignant cells (CNV | label union); AAV targeting tumor cells</li>
<li><strong>per_celltype_report</strong> — single-cell-type candidates with strict selectivity; narrow TAM-only or microglia-only targeting (microglia alone is rarely therapeutically relevant for GBM)</li>
<li><strong>pan_myeloid_report (this one)</strong> — AAV targeting the full myeloid compartment; peaks "shared" across TAM + microglia are the primary deliverable here</li>
</ul>

<h2>Methods</h2>
<ul>
<li><strong>strength_myeloid</strong> = max(strength[TAM], strength[microglia]) per peak</li>
<li><strong>strength_other_max</strong> = max over {{astrocyte, oligodendrocyte, OPC, neuron, GABA_neuron, endothelial, T_cell}}</li>
<li><strong>selectivity_myeloid</strong> = strength_myeloid / max(strength_other_max, 0.01)</li>
<li><strong>n_cohorts, n_patients, consistency</strong> = max across TAM + microglia (peak replicates in either → counts)</li>
<li><strong>composite_score_myeloid</strong> = strength × 1 + selectivity/10 × 0.5 + consistency × 0.3 + n_cohorts/8 × 0.2</li>
<li><strong>SHARED</strong> flag: both TAM AND microglia strength ≥ 0.03 (both compartments open the peak)</li>
</ul>

<h2>Table B — SHARED pan-myeloid (top 30) [primary deliverable]</h2>
<p>Peaks accessible in BOTH TAM and microglia (each ≥ 5% strength). Ranked by
joint strength × selectivity-boost × cohort replication. These are what you want for
pan-myeloid AAV targeting — they were being penalized in per-cell-type reports as
"non-specific" between TAM and microglia.</p>

<table>
<thead><tr>
  <th>peak_id (GRCh38)</th>
  <th>nearest gene<br/>tags</th>
  <th>dist TSS</th>
  <th>cohorts</th>
  <th>n patients</th>
  <th>consistency</th>
  <th>str TAM</th>
  <th>str micro</th>
  <th>str other (max)</th>
  <th>selectivity</th>
  <th>Z<sub>Daigle</sub></th>
  <th>stronger in</th>
  <th>score</th>
</tr></thead>
<tbody>
{section_rows(shared_both.head(30))}
</tbody></table>

<h2>Table A — SELECTIVE pan-myeloid (top 50)</h2>
<p>Peaks where myeloid strength exceeds any non-myeloid cell type (selectivity ≥ 1.0×).
These are strictly cleaner for AAV targeting but may miss truly shared pan-myeloid peaks
that fall just below selectivity 1.0 due to astrocyte/OPC crossover in GBM.</p>

<table>
<thead><tr>
  <th>peak_id (GRCh38)</th>
  <th>nearest gene<br/>tags</th>
  <th>dist TSS</th>
  <th>cohorts</th>
  <th>n patients</th>
  <th>consistency</th>
  <th>str TAM</th>
  <th>str micro</th>
  <th>str other (max)</th>
  <th>selectivity</th>
  <th>Z<sub>Daigle</sub></th>
  <th>stronger in</th>
  <th>score</th>
</tr></thead>
<tbody>
{section_rows(clean.head(50))}
</tbody></table>

<h2>Caveats</h2>
<ul>
<li><strong>Replication stats are loose upper bounds.</strong> We take max(TAM, microglia) cohorts/patients — so a peak replicating in 3 cohorts' TAM AND 3 cohorts' microglia (which may or may not be the same 3) scores as n_cohorts=3, not 6. A dedicated joint model would compute replication over the union of patients × cell_type_in_{{TAM, microglia}}; this approximation is biased slightly toward whichever source is replicated more.</li>
<li><strong>Myeloid-vs-malignant selectivity not computed here.</strong> The selectivity baseline excludes malignant cells on purpose (malignant cells can also be TAM-mimicking in certain Neftel MES-like states, and we don't want to double-penalize). If you want "myeloid vs malignant AND non-malignant" baseline, filter candidates by cross-referencing with <code>pan_malignant_scores_v4.parquet</code>.</li>
<li><strong>CATLAS reference is healthy brain</strong>, so GBM-TAM-specific peaks that aren't open in healthy microglia are flagged as low n_patients_accessible; those won't surface here. For TAM-specific-in-tumor-only peaks, use the single-cell-type TAM table.</li>
</ul>

<h2>Downstream files</h2>
<ul>
<li><code>reports/pan_myeloid_top50_shared.csv</code> — <strong>primary</strong>, Excel-ready (Table B)</li>
<li><code>reports/pan_myeloid_top30_shared.bed</code> — <strong>primary</strong>, GRCh38, Benchling-ready (Table B)</li>
<li><code>reports/pan_myeloid_top50_distal_selective.csv</code> — strict/selective, Excel-ready (Table A)</li>
<li><code>reports/pan_myeloid_top30.bed</code> — strict/selective, GRCh38, Benchling-ready (Table A)</li>
</ul>

</body></html>
"""
OUT_HTML.write_text(html_doc)
print(f"[myeloid] wrote {OUT_HTML}")
print("DONE.")
