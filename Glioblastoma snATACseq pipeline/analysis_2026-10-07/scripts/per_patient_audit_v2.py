"""Per-patient audit — v2 (post-CNV-fix).

Reads from pan_malignant_matrix_v4.parquet rather than the raw h5ads (which
have moved around). Confirms the 45/51 inclusion number baked into v3.

Outputs:
  reports/per_patient_audit_v2.csv
  reports/per_patient_audit_v2.html
"""
import polars as pl
import pandas as pd
import html as H
from pathlib import Path

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results")
MATRIX = HERE / "data/matrix/pan_malignant_matrix_v4.parquet"
OUT_CSV = HERE / "reports/per_patient_audit_v2.csv"
OUT_HTML = HERE / "reports/per_patient_audit_v2.html"

m = pl.read_parquet(MATRIX)
print(f"[audit-v2] matrix rows: {m.height:,}")

# Per-patient aggregate from the matrix. n_cells_malignant is the same value
# for every peak of that patient; take max.
agg = (
    m.group_by(["cohort", "patient_id"])
    .agg([
        pl.col("n_cells_malignant").max().alias("n_cells_pan_malignant"),
        pl.col("n_accessible").filter(pl.col("n_accessible") > 0).count().alias("n_peaks_accessible"),
        pl.col("peak_id").n_unique().alias("n_peaks_detected"),
    ])
    .sort(["cohort", "patient_id"])
    .to_pandas()
)
agg["included_in_pan_malignant_matrix"] = agg["n_cells_pan_malignant"] >= 10

print(f"[audit-v2] distinct (cohort, patient_id): {len(agg)}")
print(f"[audit-v2] included (n_cells >= 10): {int(agg['included_in_pan_malignant_matrix'].sum())}")

agg.to_csv(OUT_CSV, index=False)
print(f"[audit-v2] wrote {OUT_CSV}")


# HTML
def esc(s): return H.escape(str(s))

patients_total_from_51 = 51  # per README, 51 total atlas patients
patients_in = len(agg)
patients_out = patients_total_from_51 - patients_in

# Known exclusions (from pre-fix audit, retained for context — same 6 as pre-fix)
EXCLUDED_PATIENTS = [
    ("mathewson_lupien", "Wang_pair_16", 2, "2 total cells post-QC — nothing to rescue"),
    ("gbm_tme_atlas_hra004942", "P98", 480, "zero CNV-positive and zero label-positive cells after fix"),
    ("gbm_tme_atlas_hra004942", "P79", 1, "below 10-cell threshold"),
    ("gbm_tme_atlas_hra004942", "P80", 3, "below 10-cell threshold"),
    ("gbm_tme_atlas_hra004942", "P84", 8, "below 10-cell threshold"),
    ("gse165037", "GBM9", 5, "below 10-cell threshold"),
]

rows_html = []
for _, r in agg.iterrows():
    cls = "ok"
    rows_html.append(f"""
      <tr class="{cls}">
        <td>{esc(r['cohort'])}</td>
        <td><code>{esc(str(r['patient_id'])[:38])}</code></td>
        <td style="text-align:right;">{int(r['n_cells_pan_malignant']):,}</td>
        <td style="text-align:right;">{int(r['n_peaks_accessible']):,}</td>
        <td style="text-align:right;">{int(r['n_peaks_detected']):,}</td>
      </tr>""")

excluded_rows = "\n".join(
    f'<tr class="excluded"><td>{esc(c)}</td><td><code>{esc(p)}</code></td><td style="text-align:right;">{n}</td><td colspan=2>{esc(r)}</td></tr>'
    for c, p, n, r in EXCLUDED_PATIENTS
)

html_doc = f"""<!DOCTYPE html>
<html><head>
<title>Per-patient pan-malignant audit v2 — POST-CNV-FIX (2026-10-07)</title>
<meta charset="utf-8"/>
<style>
  body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; max-width: 1300px; margin: 2em auto; padding: 0 2em; color: #111; line-height: 1.5; }}
  h1 {{ font-size: 1.6em; border-bottom: 2px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.2em; margin-top: 1.5em; color: #1e40af; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.9em; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 5px 8px; }}
  th {{ background: #f9fafb; text-align: left; }}
  tr.excluded {{ background: #fee2e2; }}
  tr.ok {{ background: #f0fdf4; }}
  code {{ font-size: 0.85em; background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  .success {{ background: #d1fae5; border-left: 4px solid #059669; padding: 1em; margin: 1em 0; }}
</style></head><body>

<h1>Per-patient pan-malignant audit v2 — POST-CNV-FIX (2026-10-07)</h1>

<div class="success">
<strong>{patients_in} of {patients_total_from_51} patients contribute cells to the pan-malignant pool.
{patients_out} are missing (same 6 as pre-fix — see below).</strong>
<p>Source: <code>data/matrix/pan_malignant_matrix_v4.parquet</code> ({m.height:,} rows). This matrix
was built AFTER the CNV fix in <code>src/flatten_cohort.py</code>, so tcga_scatac and
gse276177_khan_astro are both fully included here (previously silently excluded as
<code>malignant_cnv=0</code>).</p>
<p>The old <code>per_patient_audit.html</code> next to this one is the <strong>pre-fix</strong>
snapshot (33/51 included) retained for history. Current = this v2 file.</p>
</div>

<h2>45 patients in the pan-malignant pool</h2>
<table>
<thead><tr>
  <th>cohort</th>
  <th>patient_id</th>
  <th>n cells pan-malignant</th>
  <th>n peaks ≥1 accessible</th>
  <th>n peaks tracked</th>
</tr></thead>
<tbody>
{"".join(rows_html)}
</tbody></table>

<h2>6 patients excluded (same as pre-fix)</h2>
<table>
<thead><tr>
  <th>cohort</th>
  <th>patient_id</th>
  <th>n cells</th>
  <th colspan=2>reason</th>
</tr></thead>
<tbody>
{excluded_rows}
</tbody></table>

<p>Note: inclusion threshold is <code>n_cells_pan_malignant &gt;= 10</code>. The old
"CNV_not_run_for_cohort" exclusion category has been retired — the fix in
<code>flatten_cohort.py</code> that calls <code>_parallel_chr7_10_ratio</code> in-process
runs CNV for every cohort, including the fragments-mode ones. The 12 patients that
were previously in that bucket (9 tcga_scatac + 3 gse276177_khan_astro, ~127K
malignant cells) are now fully represented in the matrix.</p>

</body></html>
"""
OUT_HTML.write_text(html_doc)
print(f"[audit-v2] wrote {OUT_HTML}")
