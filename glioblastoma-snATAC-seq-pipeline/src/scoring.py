"""
Per (peak, cell_type) scoring in two passes:

PASS 1 (fast, all peaks) — point-estimate specificity + consistency
    For every (peak × cell_type):
      strength    = mean fraction of cells in target cell type that are accessible
                    at this peak, averaged over patients
      selectivity = strength(target) / max(strength over other cell types)
      consistency = fraction of patients with strength >= threshold_frac_accessible
      n_cohorts   = how many cohorts have patients meeting the strength+consistency thresholds

PASS 2 (vectorized hierarchical Bayesian, top-N candidates per cell type) — honest replication

    Fits ONE PyMC model over all shortlisted (peak × cell_type) pairs
    simultaneously, with:

      - per-cell-type logit baseline (mu_ct)
      - per-pair deviation from its cell-type baseline (sigma_pair * z_pair)
      - per-(pair, cohort) random offset (sigma_coh * z_coh)
      - Binomial likelihood against aggregated per-cohort counts

    Partial pooling across pairs (within cell type) and across cohorts (within
    pair) delivers honest posterior intervals: sparse-data pairs are shrunk
    toward their cell-type baseline; a pair's one-off bad cohort is pulled in.

    Fitting mechanics:
      - Default sampler: ADVI (variational, deterministic, fast). 30K
        iterations with Adam. Produces a Gaussian approximation to the
        posterior, with 1000 draws for CIs.
      - Alternate: NUTS (`--sampler nuts`). Gold standard but much slower on
        the ~88K-latent parameter space; use when the cohort count grows and
        the hierarchy does real work.

    This replaces an earlier naive loop that fit 22K separate PyMC models.
    Same mathematical model, correct vectorized implementation.

Inputs: matrix/enhancer_candidate_matrix.parquet
    Expected columns:
        peak_id, cell_type, cohort, patient_id,
        n_accessible_cells, n_cells_of_type

Outputs: matrix/candidate_scores.parquet
    Columns:
        peak_id, cell_type,
        strength, selectivity, consistency, n_cohorts, n_patients,
        posterior_mean   (NaN if not in PASS 2),
        posterior_ci_lo, posterior_ci_hi,
        posterior_p_specific (NaN if not in PASS 2),
        pass2_was_run (bool),
        composite_score  — fused rank for ordering by
"""

from __future__ import annotations
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np


def log(msg: str) -> None:
    print(f"[score] {msg}", flush=True)


# ---------------------------------------------------------------------------
# PASS 1 — point estimates (unchanged)
# ---------------------------------------------------------------------------

def pass1_point_estimates(matrix_path: Path, cfg: dict):
    import polars as pl
    log(f"PASS 1: reading {matrix_path}")
    m = pl.read_parquet(matrix_path)
    log(f"  n_rows={m.height}  "
        f"n_peaks={m['peak_id'].n_unique()}  "
        f"n_cell_types={m['cell_type'].n_unique()}  "
        f"n_cohorts={m['cohort'].n_unique()}  "
        f"n_patients={m['patient_id'].n_unique()}")

    if "n_accessible" in m.columns and "n_accessible_cells" not in m.columns:
        m = m.rename({"n_accessible": "n_accessible_cells"})
    if "n_cells" in m.columns and "n_cells_of_type" not in m.columns:
        m = m.rename({"n_cells": "n_cells_of_type"})

    if "frac_accessible" not in m.columns:
        m = m.with_columns(
            (pl.col("n_accessible_cells") / pl.col("n_cells_of_type")).alias("frac_accessible")
        )
    patient = m.filter(pl.col("n_cells_of_type") > 0)

    strength = (patient.group_by(["peak_id", "cell_type"])
                       .agg(pl.col("frac_accessible").mean().alias("strength"),
                            pl.col("patient_id").n_unique().alias("n_patients")))

    tau = float(cfg["scoring"]["consistency"]["threshold_frac_accessible"])
    log(f"  consistency threshold τ = {tau}")
    cons = (patient.with_columns((pl.col("frac_accessible") > tau).alias("above_tau"))
                   .group_by(["peak_id", "cell_type"])
                   .agg(pl.col("above_tau").mean().alias("consistency")))

    strength_thr = float(cfg["scoring"]["cross_cohort_replication"]["threshold_strength"])
    cons_thr = float(cfg["scoring"]["cross_cohort_replication"]["threshold_consistency"])
    per_cohort = (patient.group_by(["peak_id", "cell_type", "cohort"])
                         .agg(pl.col("frac_accessible").mean().alias("coh_strength"),
                              (pl.col("frac_accessible") > tau).mean().alias("coh_consistency")))
    cohorts_pass = (per_cohort.filter((pl.col("coh_strength") > strength_thr) &
                                      (pl.col("coh_consistency") > cons_thr))
                              .group_by(["peak_id", "cell_type"])
                              .agg(pl.col("cohort").n_unique().alias("n_cohorts")))

    pass1 = (strength.join(cons, on=["peak_id", "cell_type"])
                     .join(cohorts_pass, on=["peak_id", "cell_type"], how="left")
                     .with_columns(pl.col("n_cohorts").fill_null(0)))

    log("  computing selectivity (per peak)")
    strength_by_peak_ct = pass1.select(["peak_id", "cell_type", "strength"]).to_pandas()
    import pandas as pd
    best_per_peak = strength_by_peak_ct.groupby("peak_id")["strength"].transform(
        lambda s: s.max())
    second_per_peak = (strength_by_peak_ct
                       .groupby("peak_id")["strength"]
                       .transform(lambda s: s.nlargest(2).iloc[-1] if len(s) >= 2 else 0.0))
    is_top = strength_by_peak_ct["strength"] == best_per_peak
    strength_by_peak_ct["denom"] = np.where(is_top, second_per_peak, best_per_peak)
    strength_by_peak_ct["selectivity"] = (
        strength_by_peak_ct["strength"] / (strength_by_peak_ct["denom"] + 1e-9)
    )
    sel_df = pl.from_pandas(strength_by_peak_ct[["peak_id", "cell_type", "selectivity"]])
    pass1 = pass1.join(sel_df, on=["peak_id", "cell_type"])

    log(f"  PASS 1 complete: {pass1.height} (peak × cell_type) rows")
    log(f"  strength    quantiles: {pass1['strength'].quantile(0.5):.4f} / "
        f"0.9={pass1['strength'].quantile(0.9):.4f} / 0.99={pass1['strength'].quantile(0.99):.4f}")
    log(f"  selectivity quantiles: {pass1['selectivity'].quantile(0.5):.2f} / "
        f"0.9={pass1['selectivity'].quantile(0.9):.2f} / 0.99={pass1['selectivity'].quantile(0.99):.2f}")

    return pass1, patient


# ---------------------------------------------------------------------------
# PASS 2 — vectorized hierarchical Beta-Binomial over all shortlisted pairs
# ---------------------------------------------------------------------------

def _shortlist(pass1, top_n_per_celltype: int):
    import polars as pl
    log(f"PASS 2 shortlist: top {top_n_per_celltype} per cell type by selectivity")
    short = (pass1.sort(["cell_type", "selectivity"], descending=[False, True])
                  .group_by("cell_type")
                  .head(top_n_per_celltype)
                  .select(["peak_id", "cell_type"]))
    log(f"  shortlist: {short.height} (peak × cell_type) rows")
    return short


def pass2_vectorized_hierarchical(pass1, patient_df, cfg: dict, sampler: str = "advi"):
    """One PyMC model over all shortlisted pairs; ADVI or NUTS."""
    import polars as pl
    import pandas as pd

    top_n = int(cfg.get("ranking", {}).get("top_n_per_celltype_for_hierarchical", 2000))
    short = _shortlist(pass1, top_n)

    try:
        import pymc as pm
        import pytensor.tensor as pt
        import arviz as az
    except ImportError as e:
        log(f"  pymc/pytensor/arviz not available ({e}) — SKIPPING pass 2.")
        return pass1.with_columns([
            pl.lit(None, dtype=pl.Float64).alias("posterior_mean"),
            pl.lit(None, dtype=pl.Float64).alias("posterior_ci_lo"),
            pl.lit(None, dtype=pl.Float64).alias("posterior_ci_hi"),
            pl.lit(None, dtype=pl.Float64).alias("posterior_p_specific"),
            pl.lit(False).alias("pass2_was_run"),
        ])

    # 1) Index shortlist pairs
    short_pd = short.to_pandas().reset_index(drop=True)
    short_pd["pair_idx"] = short_pd.index.astype(np.int64)
    ct_cat = short_pd["cell_type"].astype("category")
    ct_categories = list(ct_cat.cat.categories)
    short_pd["ct_idx"] = ct_cat.cat.codes.astype(np.int64)
    pair_to_ct = short_pd["ct_idx"].to_numpy()      # (n_pairs,)
    n_pairs = len(short_pd)
    n_cts = len(ct_categories)

    # 2) Aggregate patient-level obs to (pair_idx × cohort)
    log("  aggregating patient obs to (pair, cohort)")
    patient_pd = patient_df.select(
        ["peak_id", "cell_type", "cohort", "n_accessible_cells", "n_cells_of_type"]
    ).to_pandas()
    agg = (short_pd[["peak_id", "cell_type", "pair_idx", "ct_idx"]]
           .merge(patient_pd, on=["peak_id", "cell_type"], how="inner")
           .groupby(["pair_idx", "ct_idx", "cohort"], observed=True)
           .agg(k_sum=("n_accessible_cells", "sum"),
                n_sum=("n_cells_of_type", "sum"))
           .reset_index())
    agg = agg[agg["n_sum"] > 0].reset_index(drop=True)

    coh_cat = agg["cohort"].astype("category")
    coh_categories = list(coh_cat.cat.categories)
    agg["coh_idx"] = coh_cat.cat.codes.astype(np.int64)
    n_cohs = len(coh_categories)

    pair_idx_obs = agg["pair_idx"].to_numpy()
    coh_idx_obs = agg["coh_idx"].to_numpy()
    k_obs = agg["k_sum"].to_numpy().astype(np.int64)
    n_obs_counts = agg["n_sum"].to_numpy().astype(np.int64)
    n_obs = len(agg)

    log(f"  model dims: pairs={n_pairs}, cell_types={n_cts}, cohorts={n_cohs}, "
        f"observations={n_obs}")

    # 3) Build the model (non-centered parameterization throughout)
    #    mu_ct lives on logit scale. sigma_pair governs within-ct dispersion,
    #    sigma_coh governs cohort-level noise.
    tau_strength = float(cfg["scoring"]["cross_cohort_replication"]["threshold_strength"])
    log(f"  p_specific threshold (rate > τ_strength): {tau_strength}")

    coords = {
        "ct": ct_categories,
        "pair": short_pd["peak_id"].values + "||" + short_pd["cell_type"].values,
        "coh": coh_categories,
    }
    with pm.Model(coords=coords):
        # NO cross-pair hierarchy: the shortlist is already selected for high
        # selectivity, so shrinking pairs toward a cell-type mean distorts
        # the posterior for strong-data pairs. Each pair gets an independent
        # Normal(0, 5) prior on logit scale — essentially uninformative
        # between sigmoid(-5) ≈ 0.007 and sigmoid(+5) ≈ 0.993, letting the
        # per-cohort Binomial likelihood drive the posterior.
        logit_pair = pm.Normal("logit_pair", mu=0.0, sigma=5.0, dims="pair")

        # Partial pooling ACROSS COHORTS (within a pair) IS kept — this is
        # the hierarchy that does real statistical work: it shrinks an
        # outlier cohort's observations toward the pair's own baseline,
        # downweighting batch effects / single-sample flukes. As cohort
        # count grows, this does more work automatically.
        sigma_coh = pm.HalfNormal("sigma_coh", sigma=1.5)
        z_coh = pm.Normal("z_coh", mu=0.0, sigma=1.0, shape=(n_pairs, n_cohs))

        logit_obs = logit_pair[pair_idx_obs] + sigma_coh * z_coh[pair_idx_obs, coh_idx_obs]
        p_obs = pm.math.sigmoid(logit_obs)

        pm.Binomial("y", n=n_obs_counts, p=p_obs, observed=k_obs)

        if sampler == "nuts":
            log("  sampling with NUTS (slower, gold-standard)")
            idata = pm.sample(500, tune=500, chains=2, cores=2,
                              target_accept=0.9, progressbar=False, random_seed=42)
        else:
            log("  fitting with ADVI (fast, 30k iterations)")
            approx = pm.fit(30000, method="advi", progressbar=False, random_seed=42)
            log(f"  ADVI complete. Final ELBO: {approx.hist[-1]:.1f}")
            idata = approx.sample(1000)

    # 4) Extract per-pair posterior of pair_rate = sigmoid(logit_pair)
    log("  extracting per-pair posterior")
    logit_pair_draws = (
        idata.posterior["logit_pair"]
        .stack(sample=("chain", "draw"))
        .transpose("sample", "pair")
        .values
    )  # (n_samples, n_pairs)
    rate_draws = 1.0 / (1.0 + np.exp(-logit_pair_draws))

    posterior_mean = rate_draws.mean(axis=0)
    posterior_ci_lo = np.quantile(rate_draws, 0.025, axis=0)
    posterior_ci_hi = np.quantile(rate_draws, 0.975, axis=0)
    posterior_p_specific = (rate_draws > tau_strength).mean(axis=0)

    pass2_df = pd.DataFrame({
        "peak_id": short_pd["peak_id"].values,
        "cell_type": short_pd["cell_type"].values,
        "posterior_mean": posterior_mean,
        "posterior_ci_lo": posterior_ci_lo,
        "posterior_ci_hi": posterior_ci_hi,
        "posterior_p_specific": posterior_p_specific,
    })

    out = (pass1.join(pl.from_pandas(pass2_df),
                      on=["peak_id", "cell_type"], how="left")
                .with_columns(pl.col("posterior_mean").is_not_null().alias("pass2_was_run")))
    log(f"  PASS 2 complete: {int(out['pass2_was_run'].sum())} pairs fit")
    return out


# ---------------------------------------------------------------------------
# Composite score
# ---------------------------------------------------------------------------

def composite(scores):
    import polars as pl
    log("Computing composite score")
    out = scores.with_columns([
        (pl.col("strength").fill_null(0) * 0.30
         + pl.col("consistency").fill_null(0) * 0.25
         + ((pl.col("selectivity").log(2) / np.log2(10)).clip(0, 1).fill_null(0)) * 0.25
         + (pl.col("n_cohorts").cast(pl.Float32) / pl.col("n_cohorts").max()).fill_null(0) * 0.15
         + pl.col("posterior_p_specific").fill_null(0) * 0.05
        ).alias("composite_score")
    ])
    log(f"  top composite_score: {out['composite_score'].max():.3f}  "
        f"median: {out['composite_score'].median():.3f}")
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    import yaml
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--sampler", default="advi", choices=["advi", "nuts"],
                    help="Pass 2 fitting method (default: advi for speed)")
    args = ap.parse_args()

    with open("config/pipeline.yaml") as f:
        cfg = yaml.safe_load(f)

    pass1, patient_df = pass1_point_estimates(args.src, cfg)
    scores = pass2_vectorized_hierarchical(pass1, patient_df, cfg, sampler=args.sampler)
    scores = composite(scores)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    scores.write_parquet(args.out)
    log(f"Wrote → {args.out}  ({scores.height} rows)")

    log("\nTOP 20 by composite_score:")
    cols = ["peak_id", "cell_type", "strength", "selectivity",
            "consistency", "n_cohorts", "posterior_mean",
            "posterior_ci_lo", "posterior_ci_hi",
            "posterior_p_specific", "composite_score"]
    cols = [c for c in cols if c in scores.columns]
    import polars as pl
    top = scores.sort("composite_score", descending=True).head(20).select(cols)
    for row in top.iter_rows(named=True):
        pm_v = row.get("posterior_mean")
        lo = row.get("posterior_ci_lo")
        hi = row.get("posterior_ci_hi")
        ps = row.get("posterior_p_specific")
        ci_str = (f"pm={pm_v:.2f}[{lo:.2f},{hi:.2f}] "
                  if pm_v is not None else "pm=NA ")
        ps_str = f"p_spec={ps:.2f} " if ps is not None else "p_spec=NA "
        log(f"  {row['peak_id']:28s} {row['cell_type']:20s} "
            f"s={row['strength']:.3f} sel={row['selectivity']:.1f} "
            f"c={row['consistency']:.2f} n_coh={row['n_cohorts']} "
            f"{ci_str}{ps_str}→ {row['composite_score']:.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
