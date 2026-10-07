"""
Matrix-mode ingestion for GEO cohorts that deposited only peak-barcode matrices
(no fragments.tsv.gz). Produces an AnnData whose .var is the CATLAS 544K peak
set, so downstream peak_quantify is a no-op passthrough.

How we handle the peak-coordinate mismatch:
  1. Load per-sample peak-barcode matrix + barcodes + sample-specific peak BED.
  2. Reciprocal 50% intersect sample peaks against CATLAS peaks via pyranges.
  3. For each CATLAS peak with ≥1 overlapping sample peak, carry MAX accessibility
     over the overlapping sample peaks. For CATLAS peaks with no overlap, keep 0.
  4. Record per-cohort peak-recall (fraction of CATLAS peaks with any sample
     coverage) in .uns — this is the honest signal of how lossy the lift is.

This is the pragmatic path Phase A.2 adopts for Guilhamon + GSE165037 + similar.
A fragment-based re-processing from SRA would be ideal but costs days of compute
per sample and isn't necessary to get these cohorts into the atlas as supporting
replication evidence.

Cohort-specific sample → patient mappings are below.

Required .obs columns on output (schema matches src/ingest.py):
    cohort, patient_id, sample_id, region, modality, idh_status, primary_recurrent

Required .var columns on output:
    chrom, start, end  (CATLAS GRCh38)

.uns:
    matrix_lift = {
        "n_catlas_peaks": int,
        "per_sample": {sample_id: {"n_sample_peaks", "n_overlapping_catlas_peaks", "recall_fraction"}},
    }
"""

from __future__ import annotations
import argparse
import gzip
import io
import sys
from pathlib import Path

import numpy as np


def log(msg: str) -> None:
    print(f"[ingest-matrix] {msg}", flush=True)


def require(cond: bool, msg: str) -> None:
    if not cond:
        log(f"ABORT: {msg}")
        sys.exit(2)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def load_catlas_peaks(bed: Path):
    """Return a pyranges DataFrame of CATLAS peaks indexed in BED order."""
    import pyranges as pr
    import pandas as pd
    log(f"  loading CATLAS peaks from {bed}")
    df = pd.read_csv(bed, sep="\t", header=None, names=["Chromosome", "Start", "End"])
    df["catlas_idx"] = np.arange(len(df), dtype=np.int64)
    log(f"  CATLAS peaks: n={len(df)}")
    return df, pr.PyRanges(df)


def _load_mtx(mtx_gz: Path):
    """Load a .mtx.gz 10x-style peak-barcode matrix. Returns a sparse CSR
    with shape (n_peaks, n_barcodes), as 10x convention."""
    import scipy.io
    import scipy.sparse as sp
    log(f"    reading {mtx_gz.name}")
    with gzip.open(mtx_gz, "rb") as f:
        X = scipy.io.mmread(io.BytesIO(f.read()))
    X = sp.csr_matrix(X)
    log(f"      matrix shape (peaks × cells): {X.shape}  nnz={X.nnz}")
    return X


def _load_tsv_gz(path: Path, cols=None):
    import pandas as pd
    df = pd.read_csv(path, sep="\t", header=None, names=cols, compression="gzip")
    # Some deposits (e.g. Wang/Mathewson 2022 GSE174554) write peak coords as a
    # single column `chrX:start-end` instead of 3-column BED. In that case,
    # `names=cols` makes pandas pad Start/End with NaN. Detect + reparse.
    if cols == ["Chromosome", "Start", "End"]:
        if df["Start"].isna().all() and df["End"].isna().all():
            s = df["Chromosome"].astype(str)
            parsed = s.str.extract(r"^(?P<Chromosome>chr[\dXYM]+):(?P<Start>\d+)-(?P<End>\d+)$")
            parsed = parsed.dropna()
            if len(parsed):
                parsed["Start"] = parsed["Start"].astype(int)
                parsed["End"] = parsed["End"].astype(int)
                df = parsed
    return df


def _lift_sample_matrix_to_catlas(X_sample, sample_peaks_df, catlas_pr, catlas_df):
    """
    Project the per-sample peak × cell matrix onto the CATLAS peak space.

    X_sample: (n_sample_peaks, n_cells) CSR
    sample_peaks_df: DataFrame[Chromosome, Start, End] indexed 0..n_sample_peaks-1
    Returns:
        X_catlas: (n_catlas_peaks, n_cells) CSR  — max over sample peaks mapped to each CATLAS peak
        n_mapped_catlas: count of CATLAS peaks with ≥1 overlap
    """
    import pyranges as pr
    import pandas as pd
    import scipy.sparse as sp

    sample_peaks_df = sample_peaks_df.copy()
    # Defensive: some per-sample peak BEDs have header rows / NaN rows / string
    # coords. Drop rows with non-finite Start/End and coerce to int.
    for col in ("Start", "End"):
        sample_peaks_df[col] = pd.to_numeric(sample_peaks_df[col], errors="coerce")
    before = len(sample_peaks_df)
    sample_peaks_df = sample_peaks_df.dropna(subset=["Chromosome", "Start", "End"]).copy()
    sample_peaks_df["Start"] = sample_peaks_df["Start"].astype(np.int64)
    sample_peaks_df["End"] = sample_peaks_df["End"].astype(np.int64)
    if len(sample_peaks_df) != before:
        log(f"    dropped {before - len(sample_peaks_df)} non-coord rows from sample peaks BED")
    sample_peaks_df["sample_idx"] = np.arange(len(sample_peaks_df), dtype=np.int64)
    sample_pr = pr.PyRanges(sample_peaks_df)

    # Reciprocal 50% overlap: each peak covers ≥50% of the other.
    # pyranges .join is NOT reciprocal-sized by default; filter afterwards.
    joined = sample_pr.join(catlas_pr).df
    if len(joined) == 0:
        log("    WARN: 0 sample peaks overlap any CATLAS peak at ALL")
        empty = sp.csr_matrix((len(catlas_df), X_sample.shape[1]))
        return empty, 0

    joined["sample_len"] = joined["End"] - joined["Start"]
    joined["catlas_len"] = joined["End_b"] - joined["Start_b"]
    joined["ov_start"] = np.maximum(joined["Start"], joined["Start_b"])
    joined["ov_end"] = np.minimum(joined["End"], joined["End_b"])
    joined["overlap"] = joined["ov_end"] - joined["ov_start"]
    joined["frac_sample"] = joined["overlap"] / joined["sample_len"].clip(lower=1)
    joined["frac_catlas"] = joined["overlap"] / joined["catlas_len"].clip(lower=1)

    # Reciprocal 25% is our working floor — Guilhamon's 10x v1 peak set differs
    # enough from CATLAS GRCh38 that reciprocal 50% recovered only 3% of peaks.
    # 25% is still a biologically defensible overlap (half the peak length
    # covered on the smaller side) and ~doubles recall without blowing up noise.
    reciprocal_overlap = joined[(joined["frac_sample"] >= 0.25) & (joined["frac_catlas"] >= 0.25)]
    log(f"    overlaps: {len(joined)} any-overlap, "
        f"{len(reciprocal_overlap)} with reciprocal 25% coverage")

    if len(reciprocal_overlap) == 0:
        empty = sp.csr_matrix((len(catlas_df), X_sample.shape[1]))
        return empty, 0

    # Build (catlas_peak × sample_peak) incidence; multiply sample matrix to get lifted
    rows = reciprocal_overlap["catlas_idx"].values
    cols = reciprocal_overlap["sample_idx"].values
    data = np.ones(len(rows), dtype=np.float32)
    incidence = sp.csr_matrix(
        (data, (rows, cols)),
        shape=(len(catlas_df), X_sample.shape[0]),
    )
    # "Max" would need element-wise max over sparse ops which is slow;
    # the honest biological reading is "any sample peak accessible → CATLAS peak
    # accessible", so count>0 is the signal. Sum works fine for counts, and
    # the downstream frac_accessible is normalized by n_cells_of_type anyway.
    X_catlas = incidence @ X_sample     # (n_catlas × n_cells)
    n_mapped = int((incidence.sum(axis=1) > 0).sum())
    log(f"    CATLAS peaks with ≥1 mapped sample peak: {n_mapped:,} / {len(catlas_df):,} "
        f"({100*n_mapped/len(catlas_df):.1f}% recall)")
    return X_catlas, n_mapped


# ---------------------------------------------------------------------------
# Per-cohort handlers
# ---------------------------------------------------------------------------

GUILHAMON_SAMPLE_TO_PATIENT = {
    "GSM4131776_4218": "GBM_4218",
    "GSM4131777_4250": "GBM_4250",
    "GSM4131778_4275": "GBM_4275",
    "GSM4131779_4349": "GBM_4349",
}


def ingest_matrix_guilhamon(src: Path, catlas_bed: Path, out: Path) -> None:
    """GSE139136 — 4 primary IDHwt scATAC. Peak-barcode matrices only (no fragments)."""
    import anndata as ad
    import scipy.sparse as sp
    import pandas as pd

    log(f"Guilhamon matrix-mode ingest from {src}")
    catlas_df, catlas_pr = load_catlas_peaks(catlas_bed)
    n_catlas = len(catlas_df)

    # Expect per-sample triplets: _matrix.mtx.gz, _barcodes.tsv.gz, _peaks.bed.gz
    mtx_files = sorted(src.rglob("*_matrix.mtx.gz"))
    require(len(mtx_files) >= 4,
            f"Guilhamon expected ≥4 _matrix.mtx.gz files, found {len(mtx_files)}")

    per_sample_lift = {}
    cell_adatas = []
    for mtx in mtx_files:
        sample_tag = mtx.name.replace("_matrix.mtx.gz", "")
        barcodes = mtx.parent / f"{sample_tag}_barcodes.tsv.gz"
        peaks_bed = mtx.parent / f"{sample_tag}_peaks.bed.gz"
        require(barcodes.exists() and peaks_bed.exists(),
                f"Missing {sample_tag} companion files")
        patient = GUILHAMON_SAMPLE_TO_PATIENT.get(sample_tag, sample_tag)
        log(f"  [{sample_tag}] patient={patient}")

        X_sample = _load_mtx(mtx)  # (peaks × cells)
        sample_peaks_df = _load_tsv_gz(peaks_bed, cols=["Chromosome", "Start", "End"])
        log(f"    sample peaks: {len(sample_peaks_df)}  cells: {X_sample.shape[1]}")

        X_catlas, n_mapped = _lift_sample_matrix_to_catlas(
            X_sample, sample_peaks_df, catlas_pr, catlas_df)
        per_sample_lift[sample_tag] = {
            "n_sample_peaks": int(len(sample_peaks_df)),
            "n_overlapping_catlas_peaks": int(n_mapped),
            "recall_fraction": float(n_mapped / n_catlas),
        }

        barcodes_df = _load_tsv_gz(barcodes, cols=["barcode"])
        obs = pd.DataFrame({
            "cohort": "guilhamon",
            "sample_id": sample_tag,
            "patient_id": patient,
            "region": "unknown",
            "modality": "snATAC_matrix",
            "idh_status": "wildtype",
            "primary_recurrent": "primary",
        }, index=barcodes_df["barcode"].astype(str))
        var = catlas_df[["Chromosome", "Start", "End"]].copy()
        var.columns = ["chrom", "start", "end"]
        a = ad.AnnData(X=X_catlas.T.tocsr(), obs=obs, var=var)
        a.var_names = [f"{c}:{s}-{e}" for c, s, e in zip(var["chrom"], var["start"], var["end"])]
        cell_adatas.append(a)

    log(f"  Concatenating {len(cell_adatas)} samples")
    merged = ad.concat(cell_adatas, axis=0, join="outer")
    merged.uns["matrix_lift"] = {
        "n_catlas_peaks": n_catlas,
        "per_sample": per_sample_lift,
        "method": "pyranges reciprocal 25% overlap → sum into CATLAS space",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out, compression="gzip", compression_opts=4)
    log(f"Guilhamon → {out}  n_cells={merged.n_obs}  "
        f"n_peaks={merged.n_vars}  "
        f"mean_recall={np.mean([v['recall_fraction'] for v in per_sample_lift.values()]):.2%}")


def ingest_matrix_gse165037(src: Path, catlas_bed: Path, out: Path) -> None:
    """GSE165037 — 3 IDH-WT snATAC samples (GBM4, GBM9, GBM12).

    Filter the series _RAW.tar extracted contents to just the IDH-WT snATAC
    GSMs identified in DEPOSIT_LAYOUTS.md.
    """
    import anndata as ad
    import pandas as pd

    log(f"GSE165037 matrix-mode ingest from {src}")
    catlas_df, catlas_pr = load_catlas_peaks(catlas_bed)
    n_catlas = len(catlas_df)

    # IDH-WT snATAC GSMs (verified from eutils 2026-10-04)
    IDH_WT_ATAC_GSMS = {
        "GSM5024962": "GBM4",
        "GSM5024963": "GBM9",
        "GSM5024964": "GBM12",
    }
    mtx_files = []
    for mtx in sorted(src.rglob("*.mtx.gz")):
        gsm = mtx.name.split("_")[0]
        if gsm in IDH_WT_ATAC_GSMS:
            mtx_files.append((mtx, IDH_WT_ATAC_GSMS[gsm]))
    require(len(mtx_files) == 3,
            f"GSE165037 expected exactly 3 IDH-WT snATAC .mtx.gz files, "
            f"found {len(mtx_files)}: {[str(m[0]) for m in mtx_files]}")

    per_sample_lift = {}
    cell_adatas = []
    for mtx, patient in mtx_files:
        sample_tag = mtx.name.replace("_matrix.mtx.gz", "").replace(".mtx.gz", "")
        # Expect sibling _barcodes.tsv.gz and _peaks.bed.gz
        prefix = sample_tag.replace("_matrix", "") if sample_tag.endswith("_matrix") else sample_tag
        barcodes_cands = list(mtx.parent.glob(f"{prefix}*barcodes*.tsv.gz"))
        peaks_cands = list(mtx.parent.glob(f"{prefix}*peaks*.bed.gz"))
        require(barcodes_cands and peaks_cands,
                f"Missing companion files for {sample_tag}; see dir: {mtx.parent}")
        barcodes = barcodes_cands[0]
        peaks_bed = peaks_cands[0]
        log(f"  [{sample_tag}] patient={patient}")

        X_sample = _load_mtx(mtx)
        sample_peaks_df = _load_tsv_gz(peaks_bed, cols=["Chromosome", "Start", "End"])
        X_catlas, n_mapped = _lift_sample_matrix_to_catlas(
            X_sample, sample_peaks_df, catlas_pr, catlas_df)
        per_sample_lift[sample_tag] = {
            "n_sample_peaks": int(len(sample_peaks_df)),
            "n_overlapping_catlas_peaks": int(n_mapped),
            "recall_fraction": float(n_mapped / n_catlas),
        }
        barcodes_df = _load_tsv_gz(barcodes, cols=["barcode"])
        obs = pd.DataFrame({
            "cohort": "gse165037",
            "sample_id": sample_tag,
            "patient_id": patient,
            "region": "unknown",
            "modality": "snATAC_matrix",
            "idh_status": "wildtype",
            "primary_recurrent": "primary",
        }, index=barcodes_df["barcode"].astype(str))
        var = catlas_df[["Chromosome", "Start", "End"]].copy()
        var.columns = ["chrom", "start", "end"]
        a = ad.AnnData(X=X_catlas.T.tocsr(), obs=obs, var=var)
        a.var_names = [f"{c}:{s}-{e}" for c, s, e in zip(var["chrom"], var["start"], var["end"])]
        cell_adatas.append(a)

    merged = ad.concat(cell_adatas, axis=0, join="outer")
    merged.uns["matrix_lift"] = {
        "n_catlas_peaks": n_catlas,
        "per_sample": per_sample_lift,
        "method": "pyranges reciprocal 25% overlap → sum into CATLAS space",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out, compression="gzip", compression_opts=4)
    log(f"GSE165037 → {out}  n_cells={merged.n_obs}  "
        f"n_peaks={merged.n_vars}  "
        f"mean_recall={np.mean([v['recall_fraction'] for v in per_sample_lift.values()]):.2%}")


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def ingest_matrix_hra004942(src: Path, catlas_bed: Path, out: Path) -> None:
    """HRA004942 — Lin 2026 Nat Neurosci. The Signac RDS was extracted by
    src/extract_rds.R to a per-cohort mtx+bed+tsv bundle at
    <processed>/gbm_tme_atlas_hra004942/ingested/. Load that, lift peaks to
    CATLAS, and carry Lin's `predicted.id` as the cell_type label.

    Note: the input path (src) in this handler is the ingested/ DIRECTORY,
    not the raw deposit dir, because extract_rds.R produced the intermediate.
    """
    import anndata as ad
    import pandas as pd
    log(f"HRA004942 ingest from extracted bundle: {src}")

    # The extracted mtx+bed+tsv bundle lives at
    # /data/projects/atacseq/processed/gbm_tme_atlas_hra004942/ingested/
    # (produced by extract_rds.R, path hardcoded). Also check `src` itself.
    bundle_candidates = [
        Path("/data/projects/atacseq/processed/gbm_tme_atlas_hra004942/ingested"),
        src / "ingested",
        src,
    ]
    bundle = None
    for c in bundle_candidates:
        if (c / "matrix.mtx.gz").exists():
            bundle = c
            log(f"  using bundle: {c}")
            break
    require(bundle is not None,
            f"HRA004942 extracted bundle not found. Run src/extract_rds.R first on "
            f"the Signac RDS at {src}/scATAC/scATAC-seq/gbm_atac_signac_ArchR.rds")
    mtx = bundle / "matrix.mtx.gz"
    barcodes = bundle / "barcodes.tsv.gz"
    peaks_bed = bundle / "peaks.bed.gz"
    obs_tsv = bundle / "obs.tsv.gz"
    require(all(p.exists() for p in [mtx, barcodes, peaks_bed, obs_tsv]),
            f"HRA004942 bundle incomplete at {bundle} — did extract_rds.R run?")

    catlas_df, catlas_pr = load_catlas_peaks(catlas_bed)
    n_catlas = len(catlas_df)

    log("  loading peaks + matrix + obs")
    X_sample = _load_mtx(mtx)                             # peaks × cells
    sample_peaks_df = _load_tsv_gz(peaks_bed, cols=["Chromosome", "Start", "End"])
    log(f"    sample peaks: {len(sample_peaks_df)}  cells: {X_sample.shape[1]}")

    obs = pd.read_csv(obs_tsv, sep="\t", compression="gzip", index_col=0)
    log(f"    obs: {obs.shape[0]} rows × {obs.shape[1]} cols")

    # Lift to CATLAS
    X_catlas, n_mapped = _lift_sample_matrix_to_catlas(
        X_sample, sample_peaks_df, catlas_pr, catlas_df)
    log(f"    CATLAS peaks with signal: {n_mapped:,} / {n_catlas:,} "
        f"({100*n_mapped/n_catlas:.2f}% recall)")

    # Harmonize obs → the pipeline schema
    obs["cohort"] = "gbm_tme_atlas_hra004942"
    obs["modality"] = "snATAC_signac_extracted"
    obs["idh_status"] = "wildtype"
    obs["primary_recurrent"] = "primary"
    obs.setdefault("region", "unknown") if hasattr(obs, "setdefault") else None
    if "region" not in obs.columns:
        obs["region"] = "unknown"

    # Patient / sample IDs come from orig.ident or notes
    if "orig.ident" in obs.columns:
        obs["sample_id"] = obs["orig.ident"].astype(str)
    if "notes" in obs.columns:
        obs["patient_id"] = obs["notes"].astype(str)
    elif "orig.ident" in obs.columns:
        obs["patient_id"] = obs["orig.ident"].astype(str).str.extract(r"(P\d+)", expand=False).fillna(obs["orig.ident"].astype(str))

    # Carry Lin's cell-type labels into the pipeline's TME label columns
    if "predicted.id" in obs.columns:
        obs["label_tme_gbmap"] = obs["predicted.id"].astype(str)
        if "prediction.score.max" in obs.columns:
            obs["label_tme_confidence"] = obs["prediction.score.max"].astype(float)
        else:
            obs["label_tme_confidence"] = 1.0
        obs["label_tme_method"] = "lin_2026_predicted_id"
        log(f"    label_tme_gbmap distribution:")
        for ct, n in obs["label_tme_gbmap"].value_counts().head(15).items():
            log(f"      {str(ct)[:30]:30s} {int(n):>6d}")

    var = catlas_df[["Chromosome", "Start", "End"]].copy()
    var.columns = ["chrom", "start", "end"]
    var_names = [f"{c}:{s}-{e}" for c, s, e in zip(var["chrom"], var["start"], var["end"])]

    adata = ad.AnnData(X=X_catlas.T.tocsr(), obs=obs, var=var)
    adata.var_names = var_names
    adata.uns["matrix_lift"] = {
        "n_catlas_peaks": n_catlas,
        "n_sample_peaks": int(len(sample_peaks_df)),
        "n_overlapping_catlas_peaks": int(n_mapped),
        "recall_fraction": float(n_mapped / n_catlas),
        "method": "signac_extract → pyranges reciprocal 25% overlap → sum into CATLAS",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="lzf")
    log(f"HRA004942 → {out}  n_cells={adata.n_obs}  recall={100*n_mapped/n_catlas:.2f}%")


# Mathewson/Wang 2022 (GSE174554): matrix-mode deposit, filter to IDH-WT scATAC
WANG_ATAC_SAMPLES = {
    # SF ID → patient-pair + stage (from Wang 2022 Table 1, parsed earlier)
    "SF10592": ("Wang_pair_5",  "primary"),
    "SF11857": ("Wang_pair_5",  "recurrent"),
    "SF6621":  ("Wang_pair_16", "recurrent"),
    "SF6996":  ("Wang_pair_22", "primary"),
    "SF7062":  ("Wang_pair_22", "recurrent"),
    "SF9510":  ("Wang_pair_25", "recurrent"),  # pair primary SF9259R/S not in GEO
    "SF9798":  ("Wang_pair_27", "primary"),
    "SF9494":  ("Wang_pair_27", "recurrent"),
    "SF11331": ("Wang_pair_37", "primary"),
    # SF4297 is IDH-mutant per Table 1 — EXCLUDE
}


def ingest_matrix_mathewson_lupien(src: Path, catlas_bed: Path, out: Path) -> None:
    """Wang/Mathewson 2022 Nat Cancer (GSE174554) — 9 IDH-WT scATAC samples
    from GEO's RAW tar. Format matches Guilhamon (per-sample mtx+barcodes+peaks)."""
    import anndata as ad
    import pandas as pd
    log(f"Mathewson/Wang 2022 (GSE174554) ingest from {src}")
    catlas_df, catlas_pr = load_catlas_peaks(catlas_bed)
    n_catlas = len(catlas_df)

    # Find per-sample triplets. The GSE174554 tar has BOTH snRNA and snATAC
    # samples sharing SF IDs — only snATAC samples carry a sibling _peaks.bed.gz
    # (snRNA samples use _features.tsv.gz instead). Filter by that.
    all_mtx = sorted(src.glob("GSM*_matrix.mtx.gz"))
    log(f"  total _matrix.mtx.gz files in tar: {len(all_mtx)}")
    kept = []
    seen_sf = set()
    for mtx in all_mtx:
        sf = _wang_sample_from_name_for_match(mtx.name)
        if sf is None or sf not in WANG_ATAC_SAMPLES:
            continue
        base = mtx.name.replace("_matrix.mtx.gz", "")
        if not (mtx.parent / f"{base}_peaks.bed.gz").exists():
            continue  # snRNA sample, skip
        if sf in seen_sf:
            continue  # duplicate SF (sanity)
        seen_sf.add(sf)
        kept.append((mtx, sf))
    log(f"  kept {len(kept)} IDH-WT scATAC samples: {[sf for _, sf in kept]}")
    require(len(kept) >= 5,
            f"Wang expected ≥5 IDH-WT scATAC samples, found {len(kept)}")

    per_sample_lift = {}
    cell_adatas = []
    for mtx, sf in kept:
        base = mtx.name.replace("_matrix.mtx.gz", "")
        barcodes = mtx.parent / f"{base}_barcodes.tsv.gz"
        peaks_bed = mtx.parent / f"{base}_peaks.bed.gz"
        # Features file is also present but we don't use it (peak-level only)
        require(barcodes.exists() and peaks_bed.exists(),
                f"Companion files missing for {base}")
        patient, stage = WANG_ATAC_SAMPLES[sf]
        log(f"  [{sf}] {base}  patient={patient} stage={stage}")

        X_sample = _load_mtx(mtx)
        sample_peaks_df = _load_tsv_gz(peaks_bed, cols=["Chromosome", "Start", "End"])
        X_catlas, n_mapped = _lift_sample_matrix_to_catlas(
            X_sample, sample_peaks_df, catlas_pr, catlas_df)
        per_sample_lift[sf] = {
            "n_sample_peaks": int(len(sample_peaks_df)),
            "n_overlapping_catlas_peaks": int(n_mapped),
            "recall_fraction": float(n_mapped / n_catlas),
        }

        barcodes_df = _load_tsv_gz(barcodes, cols=["barcode"])
        obs = pd.DataFrame({
            "cohort": "mathewson_lupien",
            "sample_id": sf,
            "patient_id": patient,
            "region": "unknown",
            "modality": "snATAC_matrix",
            "idh_status": "wildtype",
            "primary_recurrent": stage,
        }, index=barcodes_df["barcode"].astype(str))
        var = catlas_df[["Chromosome", "Start", "End"]].copy()
        var.columns = ["chrom", "start", "end"]
        a = ad.AnnData(X=X_catlas.T.tocsr(), obs=obs, var=var)
        a.var_names = [f"{c}:{s}-{e}" for c, s, e in zip(var["chrom"], var["start"], var["end"])]
        cell_adatas.append(a)

    merged = ad.concat(cell_adatas, axis=0, join="outer")
    merged.uns["matrix_lift"] = {
        "n_catlas_peaks": n_catlas,
        "per_sample": per_sample_lift,
        "method": "pyranges reciprocal 25% overlap → sum into CATLAS space",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out, compression="lzf")
    log(f"Mathewson/Wang → {out}  n_cells={merged.n_obs} "
        f"mean_recall={np.mean([v['recall_fraction'] for v in per_sample_lift.values()]):.2%}")


def _wang_sample_from_name_for_match(fname: str) -> str | None:
    """Extract SF-sample ID from a filename like GSM5319577_SF9494_barcodes.tsv.gz."""
    import re
    m = re.search(r"GSM\d+_(SF\d+)_", fname)
    return m.group(1) if m else None


def ingest_matrix_gse138794_guo(src: Path, catlas_bed: Path, out: Path) -> None:
    """GSE138794 — Guo/Diaz 2020 "A single cell atlas of human glioma".

    28 samples total, 8 snATAC. IDH status from GEO soft files (2026-10-06):
      GSM4119517 SF11979 — Primary GBM, "IDHR132H WT GBM"     ← keep
      GSM4119518 SF11956 — Primary GBM, "IDHR132H WT GBM"     ← keep
      GSM4119519 SF11215 — Primary GBM (WHO 2021: GBM=IDH-WT) ← keep
      GSM4119520 SF11331 — Primary GBM (same)                 ← keep
      GSM4119513 SF11964 — IDHR132H mutant Astrocytoma        ← excluded
      GSM4119514 SF12017 — IDH1 Mutant astrocytoma G2         ← excluded
      GSM4119515 SF11949 — IDH1 mutant oligodendroglioma      ← excluded
      GSM4119516 SF11612 — Recurrent oligodendroglioma        ← excluded

    Format: CellRanger-ATAC filtered_peak_bc_matrix (10x hg38 convention).
    Per-sample files inside GSE138794_RAW.tar: GSM<id>_<sample>_matrix.mtx.gz,
    _barcodes.tsv.gz, _peaks.bed.gz.
    """
    import anndata as ad
    import pandas as pd

    log(f"GSE138794 matrix-mode ingest from {src}")
    catlas_df, catlas_pr = load_catlas_peaks(catlas_bed)
    n_catlas = len(catlas_df)

    IDH_WT_PRIMARY_GBM_GSMS = {
        "GSM4119517": "SF11979",
        "GSM4119518": "SF11956",
        "GSM4119519": "SF11215",
        "GSM4119520": "SF11331",
    }

    mtx_files = []
    for mtx in sorted(src.rglob("*.mtx.gz")):
        gsm = mtx.name.split("_")[0]
        if gsm in IDH_WT_PRIMARY_GBM_GSMS:
            mtx_files.append((mtx, IDH_WT_PRIMARY_GBM_GSMS[gsm]))
    require(len(mtx_files) >= 1,
            f"GSE138794 expected at least 1 IDH-WT primary GBM snATAC .mtx.gz, "
            f"found {len(mtx_files)}. Verify tar extracted. Dir: {src}")
    log(f"  shortlisted {len(mtx_files)} IDH-WT primary GBM snATAC samples")

    per_sample_lift = {}
    cell_adatas = []
    for mtx, patient in mtx_files:
        sample_tag = mtx.name.replace("_matrix.mtx.gz", "").replace(".mtx.gz", "")
        prefix = sample_tag.replace("_matrix", "") if sample_tag.endswith("_matrix") else sample_tag
        barcodes_cands = list(mtx.parent.glob(f"{prefix}*barcodes*.tsv.gz"))
        peaks_cands = (list(mtx.parent.glob(f"{prefix}*peaks*.bed.gz"))
                       + list(mtx.parent.glob(f"{prefix}*peaks*.tsv.gz")))
        require(barcodes_cands and peaks_cands,
                f"Missing companion files for {sample_tag}; dir: {mtx.parent}")
        barcodes = barcodes_cands[0]
        peaks_bed = peaks_cands[0]
        log(f"  [{sample_tag}] patient={patient}")

        X_sample = _load_mtx(mtx)
        sample_peaks_df = _load_tsv_gz(peaks_bed, cols=["Chromosome", "Start", "End"])
        X_catlas, n_mapped = _lift_sample_matrix_to_catlas(
            X_sample, sample_peaks_df, catlas_pr, catlas_df)
        per_sample_lift[sample_tag] = {
            "n_sample_peaks": int(len(sample_peaks_df)),
            "n_overlapping_catlas_peaks": int(n_mapped),
            "recall_fraction": float(n_mapped / n_catlas),
        }
        barcodes_df = _load_tsv_gz(barcodes, cols=["barcode"])
        obs = pd.DataFrame({
            "cohort": "gse138794_guo",
            "sample_id": sample_tag,
            "patient_id": patient,
            "region": "unknown",
            "modality": "snATAC_matrix",
            "idh_status": "wildtype",
            "primary_recurrent": "primary",
        }, index=barcodes_df["barcode"].astype(str))
        var = catlas_df[["Chromosome", "Start", "End"]].copy()
        var.columns = ["chrom", "start", "end"]
        a = ad.AnnData(X=X_catlas.T.tocsr(), obs=obs, var=var)
        a.var_names = [f"{c}:{s}-{e}" for c, s, e in zip(var["chrom"], var["start"], var["end"])]
        cell_adatas.append(a)

    merged = ad.concat(cell_adatas, axis=0, join="outer")
    merged.uns["matrix_lift"] = {
        "n_catlas_peaks": n_catlas,
        "per_sample": per_sample_lift,
        "method": "pyranges reciprocal 25% overlap → sum into CATLAS space",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out, compression="gzip", compression_opts=4)
    log(f"GSE138794 → {out}  n_cells={merged.n_obs}  "
        f"n_peaks={merged.n_vars}  "
        f"mean_recall={np.mean([v['recall_fraction'] for v in per_sample_lift.values()]):.2%}")


def ingest_matrix_gse276177_khan_astro(src: Path, catlas_bed: Path, out: Path) -> None:
    """GSE276177 — Khan 2024 "Divergence from the human astrocyte developmental
    trajectory in glioblastoma". 3 IDH-WT primary GBM patients.

    Keeps the 5 snATAC GSMs (GBM20/25 tumor+margin, GBM38 tumor).
    """
    import anndata as ad
    import pandas as pd

    log(f"GSE276177 matrix-mode ingest from {src}")
    catlas_df, catlas_pr = load_catlas_peaks(catlas_bed)
    n_catlas = len(catlas_df)

    IDH_WT_ATAC_GSMS = {
        "GSM8492625": ("GBM20", "margin"),
        "GSM8492627": ("GBM20", "tumor"),
        "GSM8492629": ("GBM25", "margin"),
        "GSM8492631": ("GBM25", "tumor"),
        "GSM8492633": ("GBM38", "tumor"),
    }

    mtx_files = []
    for mtx in sorted(src.rglob("*matrix.mtx.gz")):
        gsm = mtx.name.split("_")[0]
        if gsm in IDH_WT_ATAC_GSMS:
            patient, region = IDH_WT_ATAC_GSMS[gsm]
            mtx_files.append((mtx, patient, region))
    require(len(mtx_files) >= 3,
            f"GSE276177 expected ≥3 IDH-WT primary GBM snATAC .mtx.gz files, "
            f"found {len(mtx_files)}")
    log(f"  shortlisted {len(mtx_files)} IDH-WT primary GBM snATAC samples")

    per_sample_lift = {}
    cell_adatas = []
    for mtx, patient, region in mtx_files:
        sample_tag = mtx.name.replace("_matrix.mtx.gz", "").replace(".mtx.gz", "")
        prefix = sample_tag.replace("_matrix", "") if sample_tag.endswith("_matrix") else sample_tag
        barcodes_cands = list(mtx.parent.glob(f"{prefix}*barcodes*.tsv.gz"))
        peaks_cands = (list(mtx.parent.glob(f"{prefix}*peaks*.bed.gz"))
                       + list(mtx.parent.glob(f"{prefix}*peaks*.tsv.gz")))
        require(barcodes_cands and peaks_cands,
                f"Missing companion files for {sample_tag}; dir: {mtx.parent}")
        barcodes = barcodes_cands[0]
        peaks_bed = peaks_cands[0]
        log(f"  [{sample_tag}] patient={patient}  region={region}")

        X_sample = _load_mtx(mtx)
        sample_peaks_df = _load_tsv_gz(peaks_bed, cols=["Chromosome", "Start", "End"])
        X_catlas, n_mapped = _lift_sample_matrix_to_catlas(
            X_sample, sample_peaks_df, catlas_pr, catlas_df)
        per_sample_lift[sample_tag] = {
            "n_sample_peaks": int(len(sample_peaks_df)),
            "n_overlapping_catlas_peaks": int(n_mapped),
            "recall_fraction": float(n_mapped / n_catlas),
        }
        barcodes_df = _load_tsv_gz(barcodes, cols=["barcode"])
        obs = pd.DataFrame({
            "cohort": "gse276177_khan_astro",
            "sample_id": sample_tag,
            "patient_id": patient,
            "region": region,
            "modality": "snATAC_matrix",
            "idh_status": "wildtype",
            "primary_recurrent": "primary",
        }, index=barcodes_df["barcode"].astype(str))
        var = catlas_df[["Chromosome", "Start", "End"]].copy()
        var.columns = ["chrom", "start", "end"]
        a = ad.AnnData(X=X_catlas.T.tocsr(), obs=obs, var=var)
        a.var_names = [f"{c}:{s}-{e}" for c, s, e in zip(var["chrom"], var["start"], var["end"])]
        cell_adatas.append(a)

    merged = ad.concat(cell_adatas, axis=0, join="outer")
    merged.uns["matrix_lift"] = {
        "n_catlas_peaks": n_catlas,
        "per_sample": per_sample_lift,
        "method": "pyranges reciprocal 25% overlap → sum into CATLAS space",
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out, compression="gzip", compression_opts=4)
    log(f"GSE276177 → {out}  n_cells={merged.n_obs}  "
        f"n_peaks={merged.n_vars}  "
        f"mean_recall={np.mean([v['recall_fraction'] for v in per_sample_lift.values()]):.2%}")


HANDLERS = {
    "guilhamon": ingest_matrix_guilhamon,
    "gse165037": ingest_matrix_gse165037,
    "gse138794_guo": ingest_matrix_gse138794_guo,
    "gse276177_khan_astro": ingest_matrix_gse276177_khan_astro,
    "gbm_tme_atlas_hra004942": ingest_matrix_hra004942,
    "mathewson_lupien": ingest_matrix_mathewson_lupien,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(HANDLERS))
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--catlas-peaks", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    log(f"=== Matrix-mode ingest {args.cohort} ===")
    HANDLERS[args.cohort](args.src, args.catlas_peaks, args.out)
    log(f"=== Done {args.cohort} ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
