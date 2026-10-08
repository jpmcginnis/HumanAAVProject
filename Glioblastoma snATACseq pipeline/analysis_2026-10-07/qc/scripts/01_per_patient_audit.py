"""Per-patient audit of pan-malignant pool coverage.

Loops over all 8 cohort catlas_quantified.h5ads on Google Drive, computes
per-patient:
  - total_cells_post_qc
  - n_cells_malignant_cnv            (CNV-based call, if run for that cohort)
  - n_cells_malignant_label          (cell_type == "malignant_unresolved")
  - n_cells_pan_malignant            (CNV OR label — the inclusion rule)
  - frac_pan_malignant
  - exclusion_reason                 (free text if n_cells_pan_malignant == 0)

Outputs:
  reports/per_patient_audit.csv
  reports/per_patient_audit.html

Pre-QC counts are not recoverable (ingest-stage files were overwritten by QC
output). We report post-QC counts, which is what feeds the matrix anyway.
"""
from __future__ import annotations
import hdf5plugin  # noqa: register HDF5 filters
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

GDRIVE = Path("/Users/jpmcginnis1/Library/CloudStorage/GoogleDrive-jpmcginnis1@gmail.com/"
              "My Drive/AAV Gene Therapy/Enhancer ATACseq projects/October 2026 WT GBM analysis")
H5AD_ROOT = GDRIVE / "2026-10-6 dataset and results/data/processed_h5ads"
OUT_CSV = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results/reports/per_patient_audit.csv")
OUT_HTML = OUT_CSV.with_suffix(".html")

# Note: cohort folder "sundaram_gbm" on disk is really tcga_scatac
COHORT_DIR_TO_SLUG = {
    "sundaram_gbm": "tcga_scatac",
    "gse276177_khan_astro": "gse276177_khan_astro",
    "gbm_space": "gbm_space",
    "gbm_tme_atlas_hra004942": "gbm_tme_atlas_hra004942",
    "mathewson_lupien": "mathewson_lupien",
    "guilhamon": "guilhamon",
    "gse165037": "gse165037",
    "gse138794_guo": "gse138794_guo",
}
MIN_CELLS_FOR_PAN_MAL = 10

rows = []
for subdir in sorted(H5AD_ROOT.iterdir()):
    if not subdir.is_dir(): continue
    h5 = subdir / "catlas_quantified.h5ad"
    if not h5.exists():
        print(f"[audit] SKIP {subdir.name} — no catlas_quantified.h5ad"); continue
    cohort = COHORT_DIR_TO_SLUG.get(subdir.name, subdir.name)
    print(f"[audit] reading {subdir.name} ({cohort})")

    # backed='r' — don't load the giant X; just the obs
    a = ad.read_h5ad(h5, backed="r")
    obs = a.obs.copy()

    has_cnv_col = "malignant_cnv" in obs.columns
    if has_cnv_col:
        mal_cnv = obs["malignant_cnv"].astype(int).values
    else:
        mal_cnv = np.zeros(len(obs), dtype=int)
    mal_lab = obs["cell_type"].astype(str).values == "malignant_unresolved"
    pan_mal = (mal_cnv == 1) | mal_lab

    for pid, pid_df in obs.groupby("patient_id", observed=True):
        idx = pid_df.index
        n_total = len(idx)
        pos = obs.index.get_indexer(idx)
        n_cnv = int((mal_cnv[pos] == 1).sum())
        n_lab = int(mal_lab[pos].sum())
        n_pm = int(pan_mal[pos].sum())

        reason = ""
        if n_pm == 0:
            if not has_cnv_col or (mal_cnv == 0).all():
                reason = "CNV_not_run_for_cohort"
            elif n_cnv == 0 and n_lab == 0:
                reason = "zero_cnv_positive_zero_label"
        elif n_pm < MIN_CELLS_FOR_PAN_MAL:
            reason = f"pan_mal_cells={n_pm} < min_threshold={MIN_CELLS_FOR_PAN_MAL}"

        rows.append({
            "cohort": cohort,
            "patient_id": str(pid),
            "total_cells_post_qc": n_total,
            "n_cells_malignant_cnv": n_cnv,
            "n_cells_malignant_label": n_lab,
            "n_cells_pan_malignant": n_pm,
            "frac_pan_malignant": round(n_pm / max(n_total, 1), 4),
            "cnv_was_computed": has_cnv_col and (mal_cnv.max() > 0 if len(mal_cnv) else False),
            "included_in_pan_malignant_matrix": n_pm >= MIN_CELLS_FOR_PAN_MAL,
            "exclusion_reason": reason,
        })

    try: a.file.close()
    except Exception: pass

df = pd.DataFrame(rows).sort_values(["cohort", "patient_id"]).reset_index(drop=True)

# truncate ultra-long patient_ids (TCGA UUIDs) in display
df["patient_id_display"] = df["patient_id"].apply(lambda s: s if len(s) < 30 else s[:20] + "..." + s[-7:])

OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
df.drop(columns=["patient_id_display"]).to_csv(OUT_CSV, index=False)
print(f"\n[audit] wrote {OUT_CSV}")

# Build HTML
import html as H
def esc(s): return H.escape(str(s))
def row_cls(r):
    if r["n_cells_pan_malignant"] == 0: return "excluded"
    if not r["included_in_pan_malignant_matrix"]: return "warning"
    return "ok"

row_html = []
for _, r in df.iterrows():
    cls = row_cls(r)
    frac = f"{r['frac_pan_malignant']*100:.1f}%" if r["total_cells_post_qc"] else "—"
    cnv_badge = "✓" if r["cnv_was_computed"] else '<span style="color:#dc2626">CNV not run</span>'
    row_html.append(f"""
      <tr class="{cls}">
        <td>{esc(r['cohort'])}</td>
        <td><code title="{esc(r['patient_id'])}">{esc(r['patient_id_display'])}</code></td>
        <td style="text-align:right;">{r['total_cells_post_qc']:,}</td>
        <td style="text-align:right;">{r['n_cells_malignant_cnv']:,}</td>
        <td style="text-align:right;">{r['n_cells_malignant_label']:,}</td>
        <td style="text-align:right;"><strong>{r['n_cells_pan_malignant']:,}</strong></td>
        <td style="text-align:center;">{frac}</td>
        <td style="text-align:center;">{cnv_badge}</td>
        <td>{esc(r['exclusion_reason'])}</td>
      </tr>
    """)

total_patients = len(df)
patients_in = int(df["included_in_pan_malignant_matrix"].sum())
patients_out = total_patients - patients_in
reasons = df[~df["included_in_pan_malignant_matrix"]].groupby("exclusion_reason").size().sort_values(ascending=False)
reason_summary = "<ul>" + "".join(f"<li><strong>{int(n)} patients</strong>: <code>{esc(r)}</code></li>" for r, n in reasons.items()) + "</ul>"

html_doc = f"""<!DOCTYPE html>
<html><head>
<title>Per-patient pan-malignant audit — GBM atlas</title>
<meta charset="utf-8"/>
<style>
  body {{ font-family: -apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; max-width: 1300px; margin: 2em auto; padding: 0 2em; color: #111; line-height: 1.5; }}
  h1 {{ font-size: 1.6em; border-bottom: 2px solid #ddd; padding-bottom: 0.3em; }}
  h2 {{ font-size: 1.2em; margin-top: 1.5em; color: #1e40af; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.9em; }}
  th, td {{ border: 1px solid #e5e7eb; padding: 5px 8px; }}
  th {{ background: #f9fafb; text-align: left; }}
  tr.excluded {{ background: #fee2e2; }}
  tr.warning {{ background: #fef3c7; }}
  tr.ok {{ background: #f0fdf4; }}
  code {{ font-size: 0.85em; background: #f3f4f6; padding: 1px 4px; border-radius: 3px; }}
  .summary {{ background: #e0f2fe; border-left: 4px solid #0284c7; padding: 1em; margin: 1em 0; }}
  .note {{ background: #fef3c7; border-left: 4px solid #f59e0b; padding: 1em; margin: 1em 0; }}
</style></head><body>

<h1>Per-patient pan-malignant audit — GBM atlas 2026-10-07</h1>

<div class="summary">
<strong>{patients_in} of {total_patients} patients contribute cells to the pan-malignant pool.
{patients_out} are missing.</strong>
<p><strong>Reasons for exclusion:</strong></p>
{reason_summary}
</div>

<div class="note">
<strong>Known pipeline bug:</strong> for fragments-mode cohorts routed through
<code>src/flatten_cohort.py</code> (tcga_scatac and gse276177_khan_astro),
CNV calling was not run — <code>malignant_cnv</code> was set to a default of 0.
This silently drops those cohorts from any CNV-based malignant analysis.
Fix TODO: add CNV calling into the flatten path (currently only runs in
<code>src/cnv_malignant.py</code>, which only sees the standard
qc → cnv → labels chain).
</div>

<h2>Full per-patient table (sorted by cohort, patient_id)</h2>
<p>Row colors: 🟢 included, 🟡 included but low cell count, 🔴 excluded (0 cells in pan-malignant pool).</p>
<table>
<thead><tr>
  <th>cohort</th>
  <th>patient_id</th>
  <th>total cells (post-QC)</th>
  <th>n cells malignant_cnv=1</th>
  <th>n cells cell_type=malignant_unresolved</th>
  <th>n cells pan-malignant</th>
  <th>% pan-mal</th>
  <th>CNV run?</th>
  <th>exclusion reason</th>
</tr></thead>
<tbody>
{"".join(row_html)}
</tbody></table>

<h2>Fix path for the 12 CNV-missing patients</h2>
<p>Spin EC2 back up (<code>./infra/aws_restore.sh snap-04deb97b5e367bedf</code> → ~15 min), then for each of the two cohorts:</p>
<pre>
cd /data/projects/atacseq/claude/gbm-enhancer-atlas
source /data/miniconda3/etc/profile.d/conda.sh && conda activate atacseq

# Run CNV calling on the already-flattened h5ad in-place
python src/cnv_malignant.py \\
  --in /data/projects/atacseq/processed/sundaram_gbm/catlas_quantified.h5ad \\
  --out /data/projects/atacseq/processed/sundaram_gbm/catlas_quantified.h5ad.withcnv.h5ad
mv /data/projects/atacseq/processed/sundaram_gbm/catlas_quantified.h5ad.withcnv.h5ad \\
   /data/projects/atacseq/processed/sundaram_gbm/catlas_quantified.h5ad

# Same for khan
python src/cnv_malignant.py \\
  --in /data/projects/atacseq/processed/gse276177_khan_astro/catlas_quantified.h5ad \\
  --out /data/projects/atacseq/processed/gse276177_khan_astro/catlas_quantified.h5ad.withcnv.h5ad
mv ...

# Then re-run pan_malignant_matrix.py — should recover all 12 missing patients
python src/pan_malignant_matrix.py \\
  --inputs /data/projects/atacseq/processed/*/catlas_quantified.h5ad \\
  --matrix matrix/pan_malignant_matrix.parquet \\
  --scores matrix/pan_malignant_scores.parquet
</pre>

</body></html>
"""
OUT_HTML.write_text(html_doc)
print(f"[audit] wrote {OUT_HTML}")

print("\n=== SUMMARY ===")
print(f"total patients: {total_patients}")
print(f"included: {patients_in}")
print(f"excluded: {patients_out}")
print("\nreasons:")
for r, n in reasons.items():
    print(f"  {n:>3} patients: {r}")
