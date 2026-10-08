"""
Per-cell-type AAV enhancer candidate report.

For each cell type in the atlas, produce a ranked list of distal, selective,
cross-cohort-replicated peaks suitable for AAV cell-type-targeting.

Cell types reported: OPC, neuron, GABA_neuron, TAM, microglia, astrocyte,
oligodendrocyte, endothelial, T_cell, malignant_unresolved.
(unassigned is excluded — it's a catch-all, not a biological class.)

Scoring approach (per cell type):
  - Start from candidate_scores.parquet (per peak × cell type with strength,
    selectivity, n_cohorts, consistency, posterior pass2 fields when available)
  - Join nearest-gene annotation + atlas peak coords from pan_malignant_scores
  - Clean filter: distal (|TSS| ≥ 2 kb and ≤ 100 kb), selectivity ≥ 1.5,
    non-housekeeping, non-TME-marker (for malignant), skip chr7 (for malignant)
  - Rank by composite_score (pass2 posterior_mean when pass2_was_run, else pass1)

Produces:
  reports/per_celltype_report.html         — one HTML with sections per cell type
  reports/per_celltype/<ct>_top50.csv      — one CSV per cell type (ranked top 50)
  reports/per_celltype/<ct>_top30.bed      — GRCh38 BED, 30 best per cell type
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
OUT_HTML = HERE / "reports/per_celltype_report.html"
OUT_DIR = HERE / "reports/per_celltype"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Cell types to emit, in display order (malignant last)
CELL_TYPES = [
    ("microglia",       "Microglia (homeostatic / resident)"),
    ("TAM",             "TAM (tumor-associated macrophage)"),
    ("astrocyte",       "Astrocyte"),
    ("oligodendrocyte", "Oligodendrocyte"),
    ("OPC",             "OPC"),
    ("neuron",          "Neuron (excitatory, incl. glutamatergic)"),
    ("GABA_neuron",     "GABAergic interneuron"),
    ("endothelial",     "Endothelial"),
    ("T_cell",          "T cell"),
    ("malignant_unresolved", "Malignant (label-based; use pan-malignant v3 for CNV union)"),
]

# Gene set tags (same as v3)
GLIOMA_DRIVERS = {"EGFR","PDGFRA","MET","FGFR1","FGFR3","NTRK1","NTRK2","NTRK3","TP53","PTEN","NF1","CDKN2A","CDKN2B","RB1","PIK3CA","PIK3R1","IDH1","IDH2","ATRX","TERT","MGMT","MYC","MYCN","BRAF","CIC","FUBP1","TCF12"}
GSC_LINEAGE = {"SOX2","SOX9","SOX10","SOX11","OLIG1","OLIG2","NES","PROM1","CD44","FABP7","GFAP","VIM","MSI1","ASCL1","HES1","HES5","NOTCH1","NOTCH2","DLL3","JAG1","NFIA","NFIB","DLX1","DLX2","DLX5","DLX6","POU3F2","POU3F3"}
NEFTEL_MES = {"CD44","CHI3L1","VIM","SERPINA3","ANXA1","TGFBI","CSTB","S100A11"}
NEFTEL_AC = {"GFAP","S100B","APOE","AQP4","ALDH1L1","SLC1A3","MLC1","FGFR3","SPARCL1"}
NEFTEL_OPC = {"OLIG1","OLIG2","PDGFRA","NFIA","APOD","CSPG4","SOX10"}
NEFTEL_NPC = {"DLL3","DLX1","DLX2","DLX5","DLX6","SOX11","NEUROD4","STMN2","TUBB3"}
TAM_MARKERS = {"TMEM119","P2RY12","CX3CR1","CD68","ITGAM","AIF1","CSF1R","TREM2","IBA1","MRC1","CD163","CD206"}
ASTRO_MARKERS = {"GFAP","AQP4","SLC1A2","SLC1A3","GJA1","SOX9","NDRG2","ALDH1L1","S100B"}
OLIGO_MARKERS = {"MBP","PLP1","MOG","MAG","CNP","OLIG1","OLIG2","SOX10"}
OPC_MARKERS = {"PDGFRA","CSPG4","SOX10","OLIG2","NG2","APOD"}
NEURON_MARKERS = {"RBFOX3","NEFL","NEFM","NEFH","SYN1","SNAP25","MAP2","TUBB3","STMN2","GRIN1","GRIA1"}
ENDO_MARKERS = {"PECAM1","VWF","CDH5","CLDN5","ENG","TEK","KDR"}

HOUSEKEEPING = {"ACTB","GAPDH","TUBB","B2M","RPLP0","RPL13A","HPRT1","PPIA","UBC",
                "VCP","MCM3","MCM4","MCM5","CCT1","CCT2","CCT3","CCT4","CCT5","CCT6A","CCT7","CCT8",
                "HSPA5","HSPA8","CALR","PDI","PPIB","RAN","RAB1A",
                "HMGCS1","SREBF1","FASN",
                "POLR2A","POLR2B","POLR2C",
                "RPS6","RPL7","EIF4A1","EIF4E"}

def gene_tags(gene: str, ct: str) -> list[str]:
    """Return GBM/lineage tags for a nearest-gene, highlighting what matters for THIS cell type."""
    tags = []
    if gene in HOUSEKEEPING: tags.append("⚠️ housekeeping")
    if ct in ("microglia", "TAM") and gene in TAM_MARKERS: tags.append(f"🟢 {ct} marker")
    if ct == "astrocyte" and gene in ASTRO_MARKERS: tags.append("🟢 astrocyte marker")
    if ct in ("oligodendrocyte",) and gene in OLIGO_MARKERS: tags.append("🟢 oligo marker")
    if ct == "OPC" and gene in OPC_MARKERS: tags.append("🟢 OPC marker")
    if ct in ("neuron","GABA_neuron") and gene in NEURON_MARKERS: tags.append(f"🟢 {ct} marker")
    if ct == "endothelial" and gene in ENDO_MARKERS: tags.append("🟢 endothelial marker")
    if gene in GLIOMA_DRIVERS: tags.append("glioma driver")
    if gene in GSC_LINEAGE: tags.append("GSC/lineage TF")
    if gene in NEFTEL_MES: tags.append("Neftel MES-like")
    if gene in NEFTEL_AC: tags.append("Neftel AC-like")
    if gene in NEFTEL_OPC: tags.append("Neftel OPC-like")
    if gene in NEFTEL_NPC: tags.append("Neftel NPC-like")
    return tags


print("[per-ct] loading per-cell-type scores")
ct_df = pl.read_parquet(CELLTYPE_SCORES).to_pandas()
print(f"[per-ct]   {len(ct_df):,} rows x {len(ct_df.columns)} cols")

print("[per-ct] loading pan scores for peak coords")
pan = pl.read_parquet(PAN_SCORES, columns=["peak_id","chrom","start","end"]).to_pandas()
pan.drop_duplicates("peak_id", inplace=True)
print(f"[per-ct]   {len(pan):,} atlas peaks")

ct_df = ct_df.merge(pan, on="peak_id", how="left")

# Nearest-gene (binary search over TSSes)
print("[per-ct] building TSS index")
gene_df = pd.read_csv(GENE_BED, sep="\t", header=None,
                     names=["Chromosome","Start","End","Name","GeneType","Strand"])
gene_df["tss"] = gene_df.apply(
    lambda r: r["Start"] if r["Strand"] == "+" else r["End"], axis=1
).astype(int)
tss_by_chr = {}
for ch, grp in gene_df.groupby("Chromosome"):
    g = grp.sort_values("tss").reset_index(drop=True)
    tss_by_chr[ch] = (g["tss"].values, g["Name"].values, g["Strand"].values)

def nearest_gene(chrom, center):
    if chrom not in tss_by_chr: return (None, None)
    tss_arr, names, _ = tss_by_chr[chrom]
    idx = bisect.bisect_left(tss_arr, center)
    cands = []
    if idx > 0:            cands.append(idx - 1)
    if idx < len(tss_arr): cands.append(idx)
    best = min(cands, key=lambda i: abs(center - tss_arr[i]))
    return (names[best], int(center - tss_arr[best]))


# --- Per cell type ----------------------------------------------------------
def esc(s): return html.escape(str(s)) if pd.notna(s) else ""


def process_ct(ct, label):
    print(f"\n[per-ct] === {ct} ({label}) ===")
    sub = ct_df[ct_df["cell_type"] == ct].copy()
    if len(sub) == 0:
        print("  no rows; skip")
        return None

    # Candidate pool: n_cohorts >= 4
    sub = sub[sub["n_cohorts"] >= 4].copy()
    if len(sub) == 0:
        print("  no rows after n_cohorts>=4; skip")
        return None

    # Annotate nearest gene
    annots = []
    for _, r in sub.iterrows():
        if pd.isna(r["chrom"]):
            annots.append((None, None)); continue
        center = (int(r["start"]) + int(r["end"])) // 2
        annots.append(nearest_gene(r["chrom"], center))
    sub[["nearest_gene","dist_to_tss_signed"]] = pd.DataFrame(annots, index=sub.index)
    sub["dist_to_tss_abs_kb"] = (sub["dist_to_tss_signed"].abs() / 1000).round(1)
    sub["gene_tags_list"] = [gene_tags(g if pd.notna(g) else "", ct) for g in sub["nearest_gene"]]
    sub["gene_tags"] = sub["gene_tags_list"].map(lambda L: ", ".join(L))

    # Final rank score: pass2 posterior_mean when available, else pass1 composite_score
    sub["final_score"] = sub["posterior_mean"].fillna(sub["composite_score"])

    # Filters
    # Raw selectivity has very different ceilings per cell type in a tumor atlas:
    # neuron (small, rare in tumor) hits >2× easily; astrocyte / microglia / TAM
    # confounded with reactive glia + malignant plasticity and often top out <1×.
    # Strategy: primary filter is distal + non-housekeeping + selectivity ≥ 1.0;
    # if that yields < 10 candidates, drop the selectivity floor and rely on
    # composite_score ranking (which already weights selectivity jointly).
    is_housekeeping = sub["gene_tags"].str.contains("housekeeping", na=False)
    sel_floor_strict = 2.0 if ct == "malignant_unresolved" else 1.0
    skip_chr7 = (ct == "malignant_unresolved")
    base_mask = (
        (sub["dist_to_tss_abs_kb"] >= 2.0)
        & (sub["dist_to_tss_abs_kb"] <= 100.0)
        & (~is_housekeeping)
        & ((~skip_chr7) | (sub["chrom"] != "chr7"))
    )
    strict = sub[base_mask & (sub["selectivity"] >= sel_floor_strict)].copy()
    if len(strict) >= 10 or ct == "malignant_unresolved":
        clean = strict
        filter_mode = f"selectivity ≥ {sel_floor_strict}×, ranked by composite/posterior"
        clean = clean.sort_values("final_score", ascending=False).reset_index(drop=True)
    else:
        # Fallback: selectivity floor emptied the table (tumor-abundant cell type
        # confounded with reactive glia / malignant plasticity). To avoid showing
        # broadly accessible peaks that rank well for every cell type via raw
        # composite_score, re-rank by selectivity-weighted score so the MOST
        # cell-type-biased peaks surface — even if their absolute selectivity is <1.
        clean = sub[base_mask].copy()
        clean["selectivity_adj_score"] = clean["selectivity"].clip(lower=0.01) * clean["final_score"]
        clean = clean.sort_values(
            ["selectivity_adj_score", "selectivity"], ascending=[False, False]
        ).reset_index(drop=True)
        filter_mode = (
            f"strict (sel ≥ {sel_floor_strict}×) yielded <10; fallback ranks by "
            "selectivity × composite — these peaks are MORE open in {ct} than any other "
            "but not by enough margin for strict floor"
        ).replace("{ct}", ct)
    print(f"  pool (≥4 cohorts): {len(sub):,}  clean ({filter_mode}): {len(clean):,}")

    # CSV + BED
    out_csv = OUT_DIR / f"{ct}_top50.csv"
    out_bed = OUT_DIR / f"{ct}_top30.bed"
    cols = [
        "peak_id","chrom","start","end",
        "nearest_gene","dist_to_tss_signed","gene_tags",
        "n_cohorts","n_patients","consistency",
        "strength","selectivity",
        "composite_score","posterior_mean","posterior_ci_lo","posterior_ci_hi",
        "pass2_was_run","final_score",
    ]
    clean[cols].head(50).to_csv(out_csv, index=False)
    print(f"  wrote {out_csv}")

    with open(out_bed, "w") as f:
        for _, r in clean.head(30).iterrows():
            name = f"{r['nearest_gene']}_{int(r['dist_to_tss_signed']/1000):+d}kb_{r['peak_id']}" if pd.notna(r["nearest_gene"]) else r["peak_id"]
            f.write(f"{r['chrom']}\t{int(r['start'])}\t{int(r['end'])}\t{name}\t{int(r['final_score']*1000)}\t.\n")
    print(f"  wrote {out_bed}")

    return {"cell_type": ct, "label": label, "n_pool": len(sub),
            "n_clean": len(clean), "filter_mode": filter_mode, "top": clean}


# --- Build HTML -------------------------------------------------------------
summaries = []
for ct, label in CELL_TYPES:
    r = process_ct(ct, label)
    if r is not None:
        summaries.append(r)


def tag_html_for(tags_str):
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
        pass2 = "✓" if r.get("pass2_was_run", False) else ""
        chr7_flag = ' <span style="color:#b91c1c;" title="chr7 CNV territory">⚠️</span>' if r["chrom"] == "chr7" else ""
        out.append(f"""
        <tr>
          <td><code>{esc(r['peak_id'])}</code>{chr7_flag}</td>
          <td><strong>{esc(r['nearest_gene'])}</strong><br/>{tag_html_for(r['gene_tags'])}</td>
          <td style="text-align:right;">{d_kb:+.1f} kb</td>
          <td style="text-align:center;">{int(r['n_cohorts'])}/8</td>
          <td style="text-align:center;">{int(r['n_patients'])}</td>
          <td style="text-align:center;">{r['consistency']*100:.0f}%</td>
          <td style="text-align:center;">{r['strength']*100:.1f}%</td>
          <td style="text-align:center;"><strong>{r['selectivity']:.2f}×</strong></td>
          <td style="text-align:center;">{r['final_score']:.3f}</td>
          <td style="text-align:center;">{pass2}</td>
        </tr>""")
    return "\n".join(out)


sections_html = ""
for s in summaries:
    ct = s["cell_type"]
    label = s["label"]
    pool_sz = s["n_pool"]
    clean_sz = s["n_clean"]
    top = s["top"]
    header_note = ""
    if ct == "malignant_unresolved":
        header_note = ("<div class='note'><strong>For malignant targeting prefer the "
                       "<a href='pan_malignant_report_v3.html'>pan-malignant v3 report</a></strong> "
                       "— it includes CNV-called malignant cells that this label-only pool misses "
                       "(~438K cells in gbm_space alone).</div>")
    top30_html = section_rows(top.head(30)) if clean_sz else ""
    sections_html += f"""

<h2 id="{ct}">{esc(label)}</h2>
{header_note}
<p>Candidate pool (≥4/8 cohorts): <strong>{pool_sz:,}</strong> peaks.
Clean pool ({esc(s['filter_mode'])}, distal 2-100 kb, non-housekeeping{'/non-chr7' if ct=='malignant_unresolved' else ''}): <strong>{clean_sz:,}</strong>.
Shortlist and BED: <code>reports/per_celltype/{ct}_top50.csv</code>, <code>reports/per_celltype/{ct}_top30.bed</code>.</p>
"""
    if clean_sz:
        sections_html += f"""
<table>
<thead><tr>
  <th>peak_id (GRCh38)</th>
  <th>nearest gene / tags</th>
  <th>dist to TSS</th>
  <th>cohorts</th>
  <th>n patients</th>
  <th>consist.</th>
  <th>strength</th>
  <th>selectivity</th>
  <th>final score</th>
  <th>pass2?</th>
</tr></thead>
<tbody>
{top30_html}
</tbody></table>
"""
    else:
        sections_html += "<p><em>No peaks passed clean filters for this cell type.</em></p>"

# TOC
toc = "<ul>" + "".join(f'<li><a href="#{s["cell_type"]}">{esc(s["label"])}</a> — {s["n_clean"]:,} clean candidates</li>'
                       for s in summaries) + "</ul>"

html_doc = f"""<!DOCTYPE html>
<html><head>
<title>Per-cell-type AAV enhancer candidates — GBM atlas</title>
<meta charset="utf-8"/>
<style>
  body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; max-width: 1500px; margin: 2em auto; padding: 0 2em; color: #111; line-height: 1.5; }}
  h1 {{ font-size: 1.9em; border-bottom: 2px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.4em; margin-top: 2.5em; color: #1e40af; border-bottom: 1px solid #e5e7eb; padding-bottom: 0.2em; }}
  code {{ font-size: 0.85em; background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.85em; margin: 1em 0; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 5px 8px; vertical-align: top; }}
  th {{ background: #f9fafb; text-align: left; font-weight: 600; }}
  tr:nth-child(even) {{ background: #fafbfc; }}
  tr:hover {{ background: #eff6ff; }}
  .note {{ background: #fef3c7; border-left: 4px solid #f59e0b; padding: 1em; margin: 1em 0; }}
  .success {{ background: #d1fae5; border-left: 4px solid #059669; padding: 1em; margin: 1em 0; }}
  nav {{ background: #f9fafb; padding: 1em 1.3em; border-left: 4px solid #0284c7; margin: 1em 0; }}
  nav ul {{ margin: 0; padding-left: 1.2em; }}
</style></head><body>

<h1>Per-cell-type AAV enhancer candidates — GBM atlas (2026-10-07)</h1>

<p><strong>Source atlas:</strong> 8 cohorts, 51 patients, 1,523,668 nuclei.
Each cell type scored independently against all other cell types, with CATLAS 544,735 peaks and
hierarchical Bayesian replication (pass1 = closed-form; pass2 = PyMC posterior for top 2000 per cell type).</p>

<p><strong>How to use this report:</strong> Pick a cell type you want to target with AAV. Open its
section below for the top-30 cleaned candidates (distal, selective, non-housekeeping). Pull the
<code>_top50.csv</code> or <code>_top30.bed</code> for full detail / cloning. For malignant-cell
targeting, use the <a href="pan_malignant_report_v3.html">pan-malignant v3 report</a> (adds
CNV-called malignant cells that marker-peak scoring misses).</p>

<nav><strong>Jump to cell type:</strong>{toc}</nav>

<div class="success">
<strong>Filter definitions:</strong>
<ul>
<li><strong>n_cohorts ≥ 4</strong> — peak must be scored in at least 4/8 cohorts for cross-cohort reproducibility.</li>
<li><strong>distal</strong> — nearest TSS is 2-100 kb away. Excludes promoter-proximal peaks (which drive broad expression) and gene deserts.</li>
<li><strong>selectivity ≥ 1.0×</strong> (2.0× for malignant) — strength in target cell type exceeds max strength in any other cell type. Soft floor — tumor-microenvironment cell types in snATAC-of-GBM have low selectivity ceilings (microglia tops ~5×, neuron/OPC ~2×) because of cell-type mixing. Ranking is by composite_score (and pass2 posterior when available), which already weights selectivity jointly with strength and replication.</li>
<li><strong>non-housekeeping</strong> — nearest gene is not in a standard housekeeping panel (ACTB, GAPDH, UBC, RPL/RPS, etc.).</li>
<li><strong>pass2</strong> — hierarchical Bayesian posterior was run on this peak (top 2000 per cell type by pass1 score); when available, posterior_mean is used for final ranking.</li>
</ul>
</div>

<h2>Per-cell-type sections</h2>
{sections_html}

<h2>Caveats</h2>
<ul>
<li><strong>Selectivity ceilings differ per cell type.</strong> Microglia and oligodendrocyte reach 10×+ selectivity; neuron, OPC, TAM typically max out ≤ 4×. Compare within cell type, not across.</li>
<li><strong>TAM vs microglia overlap.</strong> Both have overlapping marker panels (CSF1R, ITGAM) and this scoring can't fully resolve them in a mixed tumor — peaks that discriminate TAM-from-microglia specifically require a dedicated contrast (not done here).</li>
<li><strong>malignant_unresolved label-based pool is incomplete</strong> for CNV-fragmented cohorts. For malignant targeting use pan_malignant v3 (adds malignant_cnv==1 cells).</li>
<li><strong>Shared cohort QC caveats</strong> — see pan_malignant_report_v3.html for cohort-level notes (hra004942 low-cell, gse165037 sci-ATAC, etc.).</li>
</ul>

</body></html>
"""

OUT_HTML.write_text(html_doc)
print(f"\n[per-ct] wrote {OUT_HTML}")
print("DONE.")
