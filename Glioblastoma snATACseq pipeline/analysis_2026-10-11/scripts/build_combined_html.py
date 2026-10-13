"""Combined HTML report — all 1200 polished top-100 distal peaks in one scrollable,
sortable, filterable table.

Data dependency: expects reports/top100_distal_combined.parquet (produced by
scripts/build_polished_reports.py from the Zenodo-deposited data:
https://doi.org/10.5281/zenodo.23252324). Run build_reports.py +
build_polished_reports.py first.
"""
import sys
from pathlib import Path
import pandas as pd
import json

ROOT = Path(__file__).resolve().parent.parent
REPORTS = ROOT / "reports"
IN_PARQUET = REPORTS / "top100_distal_combined.parquet"
OUT_HTML = REPORTS / "top100_distal_combined.html"


# Colors for the report badge
REPORT_COLORS = {
    "pan_malignant":    "#e63946",
    "pan_myeloid":      "#f77f00",
    "microglia":        "#fcbf49",
    "TAM":              "#e67e22",
    "astrocyte":        "#2a9d8f",
    "oligodendrocyte":  "#06a77d",
    "OPC":              "#52b788",
    "neuron":           "#4361ee",
    "GABA_neuron":      "#7209b7",
    "T_cell":           "#8e44ad",
    "endothelial":      "#c0392b",
    "unassigned":       "#95a5a6",
}


# Column-header + filter-label hover tooltips. Each entry: equation (Unicode)
# then a plain-English description. Shown on hover via pure CSS (no JS).
TOOLTIPS = {
    # columns
    "report": (
        "Target cell-type (or pooled bucket) this row was scored against.\n"
        "pan_malignant = all malignant cells pooled across Neftel states + unresolved.\n"
        "pan_myeloid = TAM + microglia pooled.\n"
        "Everything else is a single-cell-type target vs. all other cell types as background."
    ),
    "rank": (
        "Rank within this report, 1 = highest z_arm_adj after filters\n"
        "(distance_category ≠ promoter, target_mean_cpm ≥ 1.0, n_target_patients ≥ 2)."
    ),
    "peak_id": (
        "hg38 coordinates of the union peak: chrN:start-end.\n"
        "The union peak set is 1,958,753 peaks across 8 GBM cohorts (74% de-novo vs CATLAS)."
    ),
    "gene": (
        "RefGene symbol of the nearest annotated TSS.\n"
        "Nearest-TSS is informational only — it does NOT imply this peak regulates that gene."
    ),
    "dist_to_tss": (
        "Signed distance from peak center to the nearest TSS, in bp.\n"
        "Negative = peak sits upstream of the TSS; positive = downstream."
    ),
    "category": (
        "Distance bucket:\n"
        "  promoter   ≤ 2 kb   (filtered out of this view)\n"
        "  proximal   2–10 kb\n"
        "  distal     10–100 kb\n"
        "  far-distal > 100 kb"
    ),
    "z_arm_adj": (
        "PRIMARY RANKING COLUMN.\n"
        "Equation:  z_arm_adj = (target_mean_cpm − μ_arm) / σ_arm\n"
        "where μ_arm, σ_arm are mean and SD of non-target CPM over ALL peaks on the same chromosome arm.\n\n"
        "Why arm-adjusted: GBM malignant cells have chr7 gain / chr10 loss, which uniformly inflates\n"
        "every peak on those arms against a genome-wide null. Scoring the null per arm removes that bias\n"
        "so remaining signal reflects cell-type biology, not CNV."
    ),
    "z_raw": (
        "Equation:  z_raw = (target_mean_cpm − μ_bg) / σ_bg\n"
        "where μ_bg, σ_bg are over ALL non-target columns genome-wide (no arm correction).\n\n"
        "Daigle 2024 Z ≥ 2 is the Allen Armamentarium specificity criterion.\n"
        "Shown alongside z_arm_adj to see how much CNV inflation was removed."
    ),
    "target_cpm": (
        "Mean CPM (counts per million) across target patients, patient-level pseudobulks.\n"
        "CPM floor of 1.0 applied in this polished view — kills peaks where both target and background\n"
        "are at near-zero signal (which otherwise produce artifactually huge Z scores)."
    ),
    "n_pat": (
        "Number of target-cell-type patients with nonzero signal at this peak.\n"
        "Primary replication axis. Filtered to ≥ 2 in this view."
    ),
    "n_coh": (
        "Number of distinct cohorts those target patients span.\n"
        "n_coh ≥ 2 means the peak replicates across independent studies, not just one dataset's batch."
    ),
    "cohort_bal": (
        "Equation:  cohort_balance = max_c(signal_c) / Σ_c(signal_c)\n"
        "Fraction of target signal from the single highest-contributing cohort.\n\n"
        "1.0 = single-cohort peak (batch-dominated, risky). Rows ≥ 0.9 shaded orange.\n"
        "< 0.5 = well-balanced across cohorts."
    ),
    "offtarg_max": (
        "max CPM across any non-target cell type (patient-averaged).\n"
        "Rows > 5 shaded pink — this peak also fires in some other cell type, poor AAV restriction."
    ),
    "offtarg_top_ct": (
        "Which non-target cell type contributes offtarg_max.\n"
        "Tells you which cell type an AAV using this enhancer would also hit."
    ),
    "aav_safety": (
        "Equation:  aav_safety_ratio = target_mean_cpm / offtarget_max_cpm\n"
        "Higher = cleaner cell-type restriction. An AAV enhancer with ratio > 10 is a strong candidate;\n"
        "< 3 means significant off-target activity. Clipped at a large value when offtarg_max is near zero."
    ),
    "n_offtarg>1": (
        "Count of non-target cell types with mean CPM > 1 at this peak.\n"
        "0 = clean (only the target cell type fires); higher = broader off-target activity."
    ),
    "cnv_class": (
        "Known-malignant-CNV status of the chromosome arm:\n"
        "  chr7_gain  — commonly amplified in GBM malignant cells\n"
        "  chr10_loss — commonly deleted in GBM malignant cells\n"
        "  none       — no characteristic arm-level CNV\n\n"
        "A chr7_gain peak with high z_raw but modest z_arm_adj is likely CNV-driven, not cell-type-specific."
    ),
    "arm": (
        "Chromosome arm (hg38 centromere split), e.g. 7p, 7q, 10p, 10q.\n"
        "Used as the stratum for arm-adjusted Z."
    ),
    "%pat_chr7_elev": (
        "For malignant targets: fraction of target patients where THIS peak falls on an arm with per-patient\n"
        "elevated signal (patient-specific CNV call).\n\n"
        "High % → this peak's malignant-specificity may be driven by CNV in those patients rather than by\n"
        "cell-type biology. Blank for non-malignant targets."
    ),
    # filter labels
    "f_report": "Restrict to one target cell type / pool.",
    "f_distcat": (
        "Restrict by distance bucket (promoter is pre-filtered out of this view).\n"
        "proximal (2–10kb), distal (10–100kb), far-distal (>100kb)."
    ),
    "f_search": "Substring match against peak_id or nearest-gene symbol (case-insensitive).",
    "f_minpat": "Hide rows where fewer than this many target patients have signal.",
    "f_maxbal": (
        "Hide rows where cohort_balance exceeds this.\n"
        "Set to e.g. 0.7 to see only cross-cohort-replicating peaks."
    ),
    "f_maxoff": (
        "Hide rows where offtarget_max_cpm exceeds this.\n"
        "Set to e.g. 2 to see only cell-type-clean peaks suitable for AAV restriction."
    ),
}


def tt_attr(key, override=None):
    """HTML-escape a tooltip for use in data-tip attribute."""
    t = override if override is not None else TOOLTIPS[key]
    return (t.replace('&', '&amp;')
             .replace('"', '&quot;')
             .replace('<', '&lt;')
             .replace('>', '&gt;'))


# Reports shown in the headline combined view (same order as the UI).
HEADLINE_REPORTS = [
    "pan_malignant", "pan_myeloid",
    "astrocyte", "endothelial", "GABA_neuron", "microglia",
    "neuron", "oligodendrocyte", "OPC", "T_cell", "TAM", "unassigned",
]

NEFTEL_LABELS = {"malignant_AC", "malignant_MES", "malignant_OPC",
                 "malignant_NPC", "malignant"}


def compute_report_denominators(matrix_path):
    """For each headline report, return (n_patients, n_cohorts) that could
    contribute to the target — i.e. the maximum possible denominator for
    n_target_patients and n_target_cohorts in that report."""
    df = pd.read_parquet(matrix_path, columns=None)
    cols = [c for c in df.columns if "__" in c]
    parsed = [c.split("__", 2) for c in cols]
    atlas_patients = {(c, p) for c, p, _ in parsed}
    atlas_cohorts = {c for c, _, _ in parsed}

    def label_matches(report, label):
        if report == "pan_malignant":
            return label in NEFTEL_LABELS
        if report == "pan_myeloid":
            return label in {"TAM", "microglia"}
        return label == report

    denoms = {}
    for r in HEADLINE_REPORTS:
        targ = [(c, p) for c, p, lab in parsed if label_matches(r, lab)]
        denoms[r] = (len({(c, p) for c, p in targ}), len({c for c, _ in targ}))
    return denoms, len(atlas_patients), len(atlas_cohorts)


def fmt_dist(d):
    if abs(d) >= 1000:
        return f"{d/1000:+.1f}kb"
    return f"{d:+d}"


def main():
    df = pd.read_parquet(IN_PARQUET)
    print(f"loaded {len(df)} rows across {df['report'].nunique()} reports")

    # Per-report denominators for n_pat / n_coh tooltips.
    denoms, atlas_pts, atlas_cohs = compute_report_denominators(
        ROOT / "data" / "matrix" / "enhancer_candidate_matrix.parquet"
    )
    # Format aligned table rows: "  <report>: patients=X, cohorts=Y"
    name_w = max(len(r) for r in HEADLINE_REPORTS)
    denom_rows = "\n".join(
        f"  {r:<{name_w}s}  up to {denoms[r][0]:>2d} patients, {denoms[r][1]} cohort{'s' if denoms[r][1] != 1 else ''}"
        for r in HEADLINE_REPORTS
    )
    denom_block = (
        f"Atlas totals (patient × cell_type pseudobulks with ≥50 cells):\n"
        f"  {atlas_pts} unique (cohort, patient) across {atlas_cohs} cohorts.\n\n"
        f"Max possible denominator per report (patients where the target\n"
        f"cell type was pseudobulked at all — the ceiling for n_target_patients\n"
        f"and n_target_cohorts in that report):\n\n"
        f"{denom_rows}"
    )
    n_pat_tip = TOOLTIPS["n_pat"] + "\n\n" + denom_block
    n_coh_tip = TOOLTIPS["n_coh"] + "\n\n" + denom_block

    import math
    # Precompute per-row HTML cells (one tuple per row, JSON-serialized for the client)
    rows = []
    for _, r in df.iterrows():
        pct_elev = r.get("patient_fraction_chr7_elevated", math.nan)
        pct_elev = -1 if (pct_elev is None or (isinstance(pct_elev, float) and math.isnan(pct_elev))) else float(pct_elev)
        rows.append([
            r["report"],
            int(r["rank"]),
            r["peak_id"],
            r["nearest_gene"],
            int(r["dist_to_tss_signed"]),
            r["distance_category"],
            float(r["daigle_z_arm_adjusted"]),
            float(r.get("daigle_z_raw", 0)),
            float(r["target_mean_cpm"]),
            int(r["n_target_patients"]),
            int(r["n_target_cohorts"]),
            float(r["cohort_balance"]),
            float(r["offtarget_max_cpm"]),
            r["offtarget_top_cell_type"],
            float(r["aav_safety_ratio"]),
            int(r["n_offtarget_ct_above_1"]),
            r.get("chromosome_cnv_class", "none"),
            r.get("chromosome_arm", ""),
            pct_elev,
        ])

    report_colors_json = json.dumps(REPORT_COLORS)
    rows_json = json.dumps(rows)
    all_reports = sorted(df["report"].unique().tolist())
    all_reports_json = json.dumps(all_reports)
    all_distcats = sorted(df["distance_category"].unique().tolist())
    all_distcats_json = json.dumps(all_distcats)

    tt_report = tt_attr("report")
    tt_rank = tt_attr("rank")
    tt_peak_id = tt_attr("peak_id")
    tt_gene = tt_attr("gene")
    tt_dist_to_tss = tt_attr("dist_to_tss")
    tt_category = tt_attr("category")
    tt_z_arm_adj = tt_attr("z_arm_adj")
    tt_z_raw = tt_attr("z_raw")
    tt_target_cpm = tt_attr("target_cpm")
    tt_n_pat = tt_attr("n_pat", override=n_pat_tip)
    tt_n_coh = tt_attr("n_coh", override=n_coh_tip)
    tt_cohort_bal = tt_attr("cohort_bal")
    tt_offtarg_max = tt_attr("offtarg_max")
    tt_offtarg_top_ct = tt_attr("offtarg_top_ct")
    tt_aav_safety = tt_attr("aav_safety")
    tt_n_offtarg_gt1 = tt_attr("n_offtarg>1")
    tt_cnv_class = tt_attr("cnv_class")
    tt_arm = tt_attr("arm")
    tt_pct_elev = tt_attr("%pat_chr7_elev")
    tt_f_report = tt_attr("f_report")
    tt_f_distcat = tt_attr("f_distcat")
    tt_f_search = tt_attr("f_search")
    tt_f_minpat = tt_attr("f_minpat")
    tt_f_maxbal = tt_attr("f_maxbal")
    tt_f_maxoff = tt_attr("f_maxoff")

    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>GBM Atlas v1 — Combined Top-100 Distal</title>
<style>
*, *::before, *::after {{ box-sizing: border-box; }}
body {{
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
  margin: 0; padding: 0; color: #222; background: #f8f9fa;
  font-size: 13px;
}}
.topbar {{
  position: sticky; top: 0; z-index: 100;
  box-shadow: 0 2px 4px rgba(0,0,0,0.1);
}}
header {{
  background: #212529; color: #f8f9fa; padding: 14px 24px;
}}
header h1 {{ margin: 0 0 4px 0; font-size: 18px; font-weight: 600; }}
header .sub {{ font-size: 11px; opacity: 0.75; }}
.controls {{
  background: white; padding: 10px 24px; display: flex; gap: 12px; align-items: center;
  flex-wrap: wrap; border-bottom: 1px solid #dee2e6;
}}
.controls label {{ font-size: 11px; font-weight: 600; color: #495057; text-transform: uppercase; letter-spacing: 0.3px; }}
.controls input, .controls select {{
  padding: 5px 10px; font-size: 12px; border: 1px solid #ced4da; border-radius: 4px;
  background: white; font-family: inherit;
}}
.controls input[type="search"] {{ min-width: 180px; }}
.controls .pill {{
  background: #e9ecef; padding: 4px 10px; border-radius: 12px; font-size: 11px;
  font-family: ui-monospace, monospace; color: #495057;
}}
table {{
  border-collapse: separate; border-spacing: 0; width: 100%; background: white; font-size: 11px;
}}
thead th {{
  background: #343a40; color: #f8f9fa; padding: 8px 10px; text-align: left;
  cursor: pointer; user-select: none; position: sticky; top: var(--topbar-h, 112px); z-index: 50;
  font-weight: 600; font-size: 10px; text-transform: uppercase; letter-spacing: 0.4px;
  white-space: nowrap;
}}
thead th:hover {{ background: #495057; }}
thead th::after {{ content: ' ⇅'; opacity: 0.3; }}
thead th.asc::after {{ content: ' ▲'; opacity: 1; }}
thead th.desc::after {{ content: ' ▼'; opacity: 1; }}
td {{
  padding: 4px 10px; border-bottom: 1px solid #e9ecef;
  font-family: ui-monospace, 'SF Mono', Menlo, monospace;
  white-space: nowrap;
}}
td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
tr:hover {{ background: #f1f3f5; }}
.badge {{
  display: inline-block; padding: 2px 7px; border-radius: 10px;
  font-size: 10px; font-weight: 600; color: white; font-family: inherit;
}}
.gene {{ color: #0066cc; font-weight: 600; }}
.cat-promoter {{ color: #adb5bd; }}
.cat-proximal {{ color: #868e96; }}
.cat-distal {{ color: #087f5b; font-weight: 600; }}
.cat-far-distal {{ color: #099268; font-weight: 600; }}
.batch-warn {{ background: #ffe8cc; }}
.leak-warn {{ background: #ffd3d3; }}
#count {{ font-weight: 600; }}
.footer {{ padding: 14px 24px; font-size: 11px; color: #6c757d; background: white; border-top: 1px solid #dee2e6; }}

/* Hover tooltips: pure CSS, no JS. Uses ::before so it never conflicts with the
   sort-indicator ::after on <th>, which keeps column widths rock-steady on hover.
   NOTE: deliberately NOT setting `position: relative` on [data-tip] globally,
   because that would clobber `position: sticky` on <th>. The sticky th already
   establishes a containing block for its absolute ::before. The .controls <label>
   gets its own relative rule below. */
.controls label[data-tip] {{ position: relative; }}
[data-tip]:hover::before {{
  content: attr(data-tip);
  position: absolute;
  top: 100%;
  left: 0;
  margin-top: 4px;
  padding: 8px 10px;
  background: #212529;
  color: #f8f9fa;
  border: 1px solid #495057;
  border-radius: 4px;
  font-family: ui-monospace, 'SF Mono', Menlo, monospace;
  font-size: 11px;
  font-weight: normal;
  text-transform: none;
  letter-spacing: normal;
  white-space: pre-wrap;
  width: max-content;
  max-width: 480px;
  z-index: 200;
  box-shadow: 0 4px 12px rgba(0,0,0,0.25);
  pointer-events: none;
  line-height: 1.5;
}}
.controls label[data-tip] {{ cursor: help; border-bottom: 1px dotted #adb5bd; }}
</style></head>
<body>

<div class="topbar" id="topbar">
<header>
  <h1>GBM Enhancer Atlas v2026-10-11 — Combined Top-100 Distal</h1>
  <div class="sub">1,958,753 union peaks across 8 cohorts · 12 reports · filtered to distance_category ≠ promoter, target_mean_cpm ≥ 1.0, n_target_patients ≥ 2 · ranked by Daigle Z</div>
</header>

<div class="controls">
  <label data-tip="{tt_f_report}">Report</label>
  <select id="filter-report"><option value="">all</option></select>

  <label data-tip="{tt_f_distcat}">Distance</label>
  <select id="filter-distcat"><option value="">all</option></select>

  <label data-tip="{tt_f_search}">Search</label>
  <input type="search" id="search" placeholder="gene / peak_id">

  <label data-tip="{tt_f_minpat}">Min n_patients</label>
  <input type="number" id="filter-minpat" min="0" step="1" value="0" style="width: 60px;">

  <label data-tip="{tt_f_maxbal}">Max cohort_balance</label>
  <input type="number" id="filter-maxbal" min="0" max="1" step="0.05" value="1.0" style="width: 60px;">

  <label data-tip="{tt_f_maxoff}">Max offtarget_cpm</label>
  <input type="number" id="filter-maxoff" min="0" step="0.5" value="999" style="width: 60px;">

  <span class="pill">showing <span id="count">1200</span> of 1200</span>
  <span class="pill" style="background:#fff3bf;color:#5c3d00">hover any column header or filter label for definition</span>
</div>
</div>

<table id="t">
  <thead>
    <tr>
      <th data-key="0" data-tip="{tt_report}">report</th>
      <th data-key="1" class="num" data-tip="{tt_rank}">rank</th>
      <th data-key="2" data-tip="{tt_peak_id}">peak_id</th>
      <th data-key="3" data-tip="{tt_gene}">gene</th>
      <th data-key="4" class="num" data-tip="{tt_dist_to_tss}">dist_to_tss</th>
      <th data-key="5" data-tip="{tt_category}">category</th>
      <th data-key="6" class="num" data-tip="{tt_z_arm_adj}">z_arm_adj</th>
      <th data-key="7" class="num" data-tip="{tt_z_raw}">z_raw</th>
      <th data-key="8" class="num" data-tip="{tt_target_cpm}">target_cpm</th>
      <th data-key="9" class="num" data-tip="{tt_n_pat}">n_pat</th>
      <th data-key="10" class="num" data-tip="{tt_n_coh}">n_coh</th>
      <th data-key="11" class="num" data-tip="{tt_cohort_bal}">cohort_bal</th>
      <th data-key="12" class="num" data-tip="{tt_offtarg_max}">offtarg_max</th>
      <th data-key="13" data-tip="{tt_offtarg_top_ct}">offtarg_top_ct</th>
      <th data-key="14" class="num" data-tip="{tt_aav_safety}">aav_safety</th>
      <th data-key="15" class="num" data-tip="{tt_n_offtarg_gt1}">n_offtarg&gt;1</th>
      <th data-key="16" data-tip="{tt_cnv_class}">cnv_class</th>
      <th data-key="17" data-tip="{tt_arm}">arm</th>
      <th data-key="18" class="num" data-tip="{tt_pct_elev}">%pat_chr7_elev</th>
    </tr>
  </thead>
  <tbody id="body"></tbody>
</table>

<div class="footer">
  Sort: click a column header to sort by it (click again for descending). Hover any header for its equation + plain-language definition. Filters stack.
  Batch-dominated rows (cohort_balance ≥ 0.9) are shaded orange. Likely-off-target rows (offtarget_max_cpm &gt; 5) are shaded pink.
  Full per-report HTML + annotated CSVs live in <code>reports/detail/</code>.
</div>

<script>
const COLORS = {report_colors_json};
const REPORTS = {all_reports_json};
const DISTCATS = {all_distcats_json};
const DATA = {rows_json};

function fmtDist(d) {{
  if (Math.abs(d) >= 1000) return (d/1000 >= 0 ? '+' : '') + (d/1000).toFixed(1) + 'kb';
  return (d >= 0 ? '+' : '') + d;
}}
function badge(report) {{
  const c = COLORS[report] || '#6c757d';
  return `<span class="badge" style="background:${{c}}">${{report}}</span>`;
}}
function catCss(cat) {{ return 'cat-' + cat.replace('_','-'); }}

function render(rows) {{
  const body = document.getElementById('body');
  const html = rows.map(r => {{
    const batch = r[11] >= 0.9 ? ' class="batch-warn"' : '';
    const leak = r[12] > 5 ? ' class="leak-warn"' : '';
    const cnvClass = r[16] || '';
    const cnvBadge = cnvClass === 'none' ? '' : `<span style="background:#ffe8cc;padding:1px 4px;border-radius:3px;font-size:9px">${{cnvClass}}</span>`;
    const pctElev = r[18] < 0 ? '' : (r[18] * 100).toFixed(0) + '%';
    return `<tr${{batch || leak}}>
      <td>${{badge(r[0])}}</td>
      <td class="num">${{r[1]}}</td>
      <td>${{r[2]}}</td>
      <td class="gene">${{r[3]}}</td>
      <td class="num">${{fmtDist(r[4])}}</td>
      <td class="${{catCss(r[5])}}">${{r[5]}}</td>
      <td class="num"><b>${{r[6].toFixed(2)}}</b></td>
      <td class="num">${{r[7].toFixed(1)}}</td>
      <td class="num">${{r[8].toFixed(2)}}</td>
      <td class="num">${{r[9]}}</td>
      <td class="num">${{r[10]}}</td>
      <td class="num"${{batch}}>${{r[11].toFixed(2)}}</td>
      <td class="num"${{leak}}>${{r[12].toFixed(2)}}</td>
      <td>${{r[13]}}</td>
      <td class="num">${{r[14].toFixed(1)}}</td>
      <td class="num">${{r[15]}}</td>
      <td>${{cnvBadge}}</td>
      <td>${{r[17]}}</td>
      <td class="num">${{pctElev}}</td>
    </tr>`;
  }}).join('');
  body.innerHTML = html;
  document.getElementById('count').textContent = rows.length;
}}

// Populate filter dropdowns
const reportSel = document.getElementById('filter-report');
REPORTS.forEach(r => {{
  const o = document.createElement('option');
  o.value = r; o.textContent = r;
  reportSel.appendChild(o);
}});
const distcatSel = document.getElementById('filter-distcat');
DISTCATS.forEach(c => {{
  const o = document.createElement('option');
  o.value = c; o.textContent = c;
  distcatSel.appendChild(o);
}});

let sortKey = 6, sortDir = -1; // default sort by z_arm_adj desc

function applyFiltersAndSort() {{
  const rep = reportSel.value;
  const cat = distcatSel.value;
  const q = document.getElementById('search').value.toLowerCase();
  const minPat = parseInt(document.getElementById('filter-minpat').value) || 0;
  const maxBal = parseFloat(document.getElementById('filter-maxbal').value) || 1.0;
  const maxOff = parseFloat(document.getElementById('filter-maxoff').value) || 999;
  let rows = DATA.filter(r =>
    (!rep || r[0] === rep) &&
    (!cat || r[5] === cat) &&
    (!q || r[2].toLowerCase().includes(q) || r[3].toLowerCase().includes(q)) &&
    r[9] >= minPat &&
    r[11] <= maxBal &&
    r[12] <= maxOff
  );
  rows.sort((a, b) => {{
    const av = a[sortKey], bv = b[sortKey];
    if (typeof av === 'number') return (av - bv) * sortDir;
    return String(av).localeCompare(String(bv)) * sortDir;
  }});
  render(rows);
}}

document.querySelectorAll('thead th').forEach((th, i) => {{
  th.addEventListener('click', () => {{
    document.querySelectorAll('thead th').forEach(x => x.classList.remove('asc','desc'));
    const key = parseInt(th.dataset.key);
    if (sortKey === key) sortDir *= -1;
    else {{ sortKey = key; sortDir = -1; }}
    th.classList.add(sortDir === 1 ? 'asc' : 'desc');
    applyFiltersAndSort();
  }});
}});

['filter-report','filter-distcat','search','filter-minpat','filter-maxbal','filter-maxoff'].forEach(id => {{
  document.getElementById(id).addEventListener('input', applyFiltersAndSort);
}});

applyFiltersAndSort();
// mark default sort indicator
document.querySelectorAll('thead th')[6].classList.add('desc');

// Keep the sticky thead flush under the (possibly multi-line) topbar.
function syncTopbarHeight() {{
  const h = document.getElementById('topbar').offsetHeight;
  document.documentElement.style.setProperty('--topbar-h', h + 'px');
}}
syncTopbarHeight();
window.addEventListener('resize', syncTopbarHeight);
new ResizeObserver(syncTopbarHeight).observe(document.getElementById('topbar'));
</script></body></html>
"""
    OUT_HTML.write_text(html)
    print(f"wrote {OUT_HTML}  ({OUT_HTML.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    sys.exit(main() or 0)
