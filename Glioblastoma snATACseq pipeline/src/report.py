"""
Render HTML reports using Jinja templates and the Tufte aesthetic.
Each page corresponds to one Snakemake rule output.

Pages:
    index                — landing page with the top-10 clone list and verdict
    top_candidates       — full top-N table per cell type, sortable
    celltype <ct>        — per-cell-type deep dive
    surprises            — grouped by detector
    dataset_issues       — coverage audit and the "do we need our own data?" verdict
"""

from __future__ import annotations
import argparse, sys
from pathlib import Path


def log(msg): print(f"[report] {msg}", flush=True)


def tufte_css() -> str:
    # Minimal Tufte-inspired CSS, inlined to keep reports self-contained.
    return """
    body { font-family: Georgia, serif; max-width: 1100px; margin: 2rem auto;
           padding: 0 1.5rem; color: #111; line-height: 1.55; }
    h1 { font-size: 1.9rem; font-weight: 400; border-bottom: 1px solid #ccc; padding-bottom: 0.3rem; }
    h2 { font-size: 1.3rem; font-weight: 400; margin-top: 2.5rem; color: #333; }
    table { border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.92rem; }
    th { text-align: left; border-bottom: 1px solid #333; padding: 0.3rem 0.6rem; font-weight: 500; }
    td { padding: 0.25rem 0.6rem; border-bottom: 1px solid #eee; }
    tr:hover { background: #fafafa; }
    .verdict { padding: 1rem 1.3rem; border-left: 3px solid; margin: 1.5rem 0; }
    .verdict.pass { border-color: #2a7a2a; background: #f3faf3; }
    .verdict.fail { border-color: #b84a2a; background: #fcf5f3; }
    .caption { color: #666; font-style: italic; font-size: 0.88rem; margin-top: 0.3rem; }
    code { background: #f5f5f5; padding: 1px 4px; border-radius: 2px; font-size: 0.88em; }
    .footer { color: #888; font-size: 0.82rem; margin-top: 3rem; border-top: 1px solid #ccc; padding-top: 0.6rem; }
    """


def render_index(inputs, out):
    import polars as pl
    top = pl.read_parquet(inputs[0])
    issues = pl.read_parquet(inputs[1]) if len(inputs) > 1 else pl.DataFrame()

    clone_first = top.sort("composite_score", descending=True).head(10)
    rows = "\n".join(
        f"<tr><td>{i+1}</td><td><code>{r['peak_id']}</code></td><td>{r['cell_type']}</td>"
        f"<td>{r['composite_score']:.3f}</td><td>{r['strength']:.2f}</td>"
        f"<td>{r['consistency']:.2f}</td><td>{r['selectivity']:.1f}</td>"
        f"<td>{r.get('n_cohorts_replicating', r.get('n_cohorts', 0))}/{r.get('n_cohorts', 0)}</td></tr>"
        for i, r in enumerate(clone_first.iter_rows(named=True))
    )

    verdict_cls = "pass" if issues.filter(pl.col("severity") == "high").height == 0 else "fail"
    verdict_msg = ("Existing public data suffices. Write the R01 for validation."
                   if verdict_cls == "pass"
                   else "Public data has significant coverage gaps. Prospective n=24 Multiome cohort recommended.")

    html = f"""<!doctype html>
<html><head><meta charset="utf-8"><title>GBM Enhancer Atlas</title>
<style>{tufte_css()}</style></head><body>
<h1>GBM Enhancer Atlas</h1>
<p class="caption">Harmonized primary-tissue IDH-WT GBM snATAC + Multiome. Updated by the Snakemake pipeline on each run.</p>

<div class="verdict {verdict_cls}">
  <strong>Verdict:</strong> {verdict_msg}
</div>

<h2>Clone these first — top 10 overall</h2>
<table>
<thead><tr><th>#</th><th>peak</th><th>cell type</th><th>composite</th><th>strength</th>
<th>consistency</th><th>selectivity</th><th>cohorts</th></tr></thead>
<tbody>{rows}</tbody>
</table>
<p class="caption">Composite score is the weighted combination from <code>config/pipeline.yaml</code>.
Change the weights there and rerun <code>snakemake rank_top_n --forcerun</code>.</p>

<h2>Sections</h2>
<ul>
  <li><a href="top_candidates.html">Full top-N per cell type</a></li>
  <li><a href="surprises.html">Surprises — unexpected patterns worth knowing</a></li>
  <li><a href="dataset_issues.html">Dataset coverage audit — do we need our own data?</a></li>
</ul>

<p class="footer">Rendered by src/report.py — Tufte aesthetic, no chartjunk. Query the matrix directly from
<code>matrix/enhancer_candidate_matrix.parquet</code> for custom analyses.</p>
</body></html>"""
    out.write_text(html)
    log(f"Wrote index → {out}")


def render_top(inputs, out):
    import polars as pl
    top = pl.read_parquet(inputs[0])
    sections = []
    for ct in sorted(set(top["cell_type"])):
        sub = top.filter(pl.col("cell_type") == ct).sort("composite_score", descending=True)
        rows = "\n".join(
            f"<tr><td><code>{r['peak_id']}</code></td><td>{r['composite_score']:.3f}</td>"
            f"<td>{r['strength']:.2f}</td><td>{r['consistency']:.2f}</td>"
            f"<td>{r['selectivity']:.1f}</td><td>{r.get('n_cohorts_replicating', r.get('n_cohorts', 0))}</td></tr>"
            for r in sub.iter_rows(named=True)
        )
        sections.append(f"<h2>{ct}</h2><table>"
                        f"<thead><tr><th>peak</th><th>composite</th><th>strength</th>"
                        f"<th>consistency</th><th>selectivity</th><th>n cohorts replicating</th></tr></thead>"
                        f"<tbody>{rows}</tbody></table>")
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Top candidates per cell type</title><style>{tufte_css()}</style></head>
<body><h1>Top candidates per cell type</h1>
<p><a href="index.html">&larr; back to index</a></p>
{''.join(sections)}
</body></html>"""
    out.write_text(html)
    log(f"Wrote top candidates → {out}")


def render_surprises(inputs, out):
    import polars as pl
    sur = pl.read_parquet(inputs[0])
    sections = []
    for detector in sorted(set(sur["detector"])) if "detector" in sur.columns else []:
        sub = sur.filter(pl.col("detector") == detector).head(50)
        cols = [c for c in sub.columns if c != "detector"]
        head = "".join(f"<th>{c}</th>" for c in cols)
        rows = "\n".join(
            "<tr>" + "".join(f"<td>{r.get(c, '')}</td>" for c in cols) + "</tr>"
            for r in sub.iter_rows(named=True)
        )
        sections.append(f"<h2>{detector}</h2><table><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table>")
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Surprises</title><style>{tufte_css()}</style></head>
<body><h1>Surprises — patterns not explicitly asked for</h1>
<p><a href="index.html">&larr; back to index</a></p>
{''.join(sections) if sections else '<p>No surprises detected on this run.</p>'}
</body></html>"""
    out.write_text(html)
    log(f"Wrote surprises → {out}")


def render_dataset_issues(inputs, out):
    import polars as pl
    df = pl.read_parquet(inputs[0])
    high = df.filter(pl.col("severity") == "high") if df.height > 0 else df
    verdict_cls = "pass" if high.height == 0 else "fail"
    verdict_msg = ("Existing public data suffices — write for validation."
                   if verdict_cls == "pass"
                   else f"{high.height} high-severity gaps — prospective n=24 Multiome cohort recommended.")
    rows = "\n".join(
        f"<tr><td>{r['cell_type']}</td><td>{r['issue']}</td>"
        f"<td>{r['observed']}</td><td>{r['required']}</td>"
        f"<td>{r['severity']}</td></tr>"
        for r in df.iter_rows(named=True)
    ) if df.height > 0 else "<tr><td colspan='5'>No issues detected.</td></tr>"
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Dataset issues</title><style>{tufte_css()}</style></head>
<body><h1>Dataset coverage audit</h1>
<p><a href="index.html">&larr; back to index</a></p>
<div class="verdict {verdict_cls}"><strong>{verdict_msg}</strong></div>
<table><thead><tr><th>cell type</th><th>issue</th><th>observed</th><th>required</th><th>severity</th></tr></thead>
<tbody>{rows}</tbody></table>
</body></html>"""
    out.write_text(html)
    log(f"Wrote dataset issues → {out}")


def render_celltype(inputs, out, ct):
    import polars as pl
    df = pl.read_parquet(inputs[0])
    rows = "\n".join(
        f"<tr><td><code>{r['peak_id']}</code></td><td>{r['composite_score']:.3f}</td>"
        f"<td>{r['strength']:.2f}</td><td>{r['consistency']:.2f}</td>"
        f"<td>{r['selectivity']:.1f}</td></tr>"
        for r in df.iter_rows(named=True)
    )
    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>{ct} candidates</title><style>{tufte_css()}</style></head>
<body><h1>{ct} — top enhancer candidates</h1>
<p><a href="../index.html">&larr; back to index</a></p>
<table><thead><tr><th>peak</th><th>composite</th><th>strength</th><th>consistency</th><th>selectivity</th></tr></thead>
<tbody>{rows}</tbody></table>
</body></html>"""
    out.write_text(html)
    log(f"Wrote {ct} page → {out}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--page", required=True,
                    choices=["index", "top_candidates", "surprises", "dataset_issues", "celltype"])
    ap.add_argument("--inputs", nargs="+", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--celltype", default=None)
    args = ap.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.page == "index":
        render_index(args.inputs, args.out)
    elif args.page == "top_candidates":
        render_top(args.inputs, args.out)
    elif args.page == "surprises":
        render_surprises(args.inputs, args.out)
    elif args.page == "dataset_issues":
        render_dataset_issues(args.inputs, args.out)
    elif args.page == "celltype":
        render_celltype(args.inputs, args.out, args.celltype)
    return 0


if __name__ == "__main__":
    sys.exit(main())
