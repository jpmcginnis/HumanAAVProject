#!/usr/bin/env python3
"""
Step 3 of the pipeline: DE results → drift headline + figures.

Reads limma_voom_de.csv and intersects the significantly-DOWN gene set
with Gosselin's 477-gene conserved microglia signature to produce the
primary drift metric:

    fraction of the conserved signature dropping >2-fold at FDR<0.05

Alongside the headline, produces:
  * a per-gene table restricted to the 477 signature with per-gene
    log2FC, FDR, and DOWN/UP/NS call
  * a scatter of Gosselin's own per-gene logFC (from Table S6) vs ours,
    if the reference file is available (for concordance sanity)
  * a bar chart comparing the drift headline to Gosselin's 33 %

Usage:
  python 03_drift_summary.py --de output/limma_voom_de.csv \\
                              --signature references/gosselin_conserved_signature_477.csv \\
                              --outdir output/
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GOSSELIN_DRIFT_PCT = 0.33   # Gosselin 2017 published headline


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--de", required=True, help="limma_voom_de.csv")
    ap.add_argument("--signature", required=True,
                     help="CSV with column 'human_gene'")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--label", default="this study",
                     help="label for the y-axis of the comparison bar chart")
    ap.add_argument("--gosselin-per-gene", default=None,
                     help="optional: gosselin's own per-gene logFC CSV "
                          "for concordance scatter")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    de = pd.read_csv(args.de)
    sig_genes = pd.read_csv(args.signature)["human_gene"].tolist()

    de = de.set_index("gene")
    present = [g for g in sig_genes if g in de.index]
    sub = de.loc[present].copy()
    sub["call"] = np.select(
        [(sub["adj.P.Val"] < 0.05) & (sub["logFC"] < -1),
         (sub["adj.P.Val"] < 0.05) & (sub["logFC"] >  1)],
        ["DOWN", "UP"], default="NS"
    )
    sub.to_csv(outdir / "drift_per_gene_signature.csv")

    n_input = len(sig_genes)
    n_present = len(sub)
    n_down = int((sub["call"] == "DOWN").sum())
    n_up   = int((sub["call"] == "UP").sum())
    frac_down = n_down / n_present if n_present else float("nan")
    frac_up   = n_up   / n_present if n_present else float("nan")

    print("=== Drift headline ===")
    print(f"  signature: {n_input} genes ; {n_present} present in our DE test")
    print(f"  DROPPED >2-fold at FDR<.05: {n_down} ({frac_down:.1%})")
    print(f"  ROSE    >2-fold at FDR<.05: {n_up} ({frac_up:.1%})")
    print(f"  Gosselin 2017 reference (dissociated 7 d): {GOSSELIN_DRIFT_PCT:.0%} DOWN")

    with open(outdir / "drift_headline.txt", "w") as fh:
        fh.write(f"drift_headline_signature_input\t{n_input}\n")
        fh.write(f"drift_headline_signature_present\t{n_present}\n")
        fh.write(f"n_down_gt2fold_FDR05\t{n_down}\n")
        fh.write(f"n_up_gt2fold_FDR05\t{n_up}\n")
        fh.write(f"frac_down_gt2fold_FDR05\t{frac_down}\n")
        fh.write(f"frac_up_gt2fold_FDR05\t{frac_up}\n")
        fh.write(f"gosselin_reference_frac_down\t{GOSSELIN_DRIFT_PCT}\n")

    # Bar chart: our number vs Gosselin's 33 %
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    labels = [f"Gosselin 2017\n(dissociated 7 d)", args.label]
    vals = [GOSSELIN_DRIFT_PCT, frac_down]
    colors = ["#888", "#2ca02c"]
    bars = ax.bar(labels, vals, color=colors, edgecolor="black")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width()/2, v + 0.008, f"{v:.1%}",
                 ha="center", fontsize=11)
    ax.set_ylabel("Fraction of Gosselin's 477-gene conserved microglia signature\n"
                   "dropping > 2-fold in vitro  (FDR<0.05)")
    ax.set_ylim(0, max(0.5, max(vals) * 1.35))
    ax.set_title("Drift headline: this pipeline on user data vs Gosselin's published")
    plt.tight_layout()
    plt.savefig(outdir / "drift_headline_bar.png", dpi=300)
    plt.close()

    # Per-gene concordance scatter, if we have Gosselin's per-gene logFC
    if args.gosselin_per_gene and Path(args.gosselin_per_gene).exists():
        g = pd.read_csv(args.gosselin_per_gene).set_index("gene")
        col = [c for c in g.columns if "log2FC" in c][0]
        common = sorted(set(g.index) & set(de.index))
        gx = g.loc[common, col].values
        ox = de.loc[common, "logFC"].values
        r = float(np.corrcoef(gx, ox)[0, 1])
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(gx, ox, s=6, alpha=0.4, color="#444", rasterized=True)
        lo, hi = -8, 8
        ax.plot([lo, hi], [lo, hi], "--", color="grey", lw=0.6)
        ax.axhline(0, color="grey", lw=0.4); ax.axvline(0, color="grey", lw=0.4)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("Gosselin log2FC (7 d in-vitro − ex-vivo)")
        ax.set_ylabel(f"{args.label}: log2FC")
        ax.set_title(f"Per-gene concordance ({len(common)} genes)   Pearson r = {r:.2f}")
        plt.tight_layout()
        plt.savefig(outdir / "gosselin_concordance_scatter.png", dpi=300)
        plt.close()
        print(f"  wrote gosselin_concordance_scatter.png  (r = {r:.3f})")

    print(f"\nAll outputs in {outdir}")


if __name__ == "__main__":
    main()
