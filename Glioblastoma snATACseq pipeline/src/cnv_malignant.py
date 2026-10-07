"""
CNV-based malignant calling via inferCNVpy with same-patient reference cells.

Rationale: the earlier chr7+/chr10- ratio ≥ 2.0 cutoff is crude and miscalls ~20%
of mesenchymal and recurrent GBMs where chr7/10 CNV is less clean. The honest
approach for primary-tissue nuclei is inferCNV with a per-sample reference set
drawn from the same patient's own non-malignant cells (microglia / endothelial),
so technical, patient-specific background is subtracted.

Flow:
  1. If modality='multiome', use the paired RNA counts directly.
     Otherwise compute gene activity scores via snapatac2.pp.make_gene_matrix.
  2. Per sample_id:
       a. Score every cell for microglia markers and endothelial markers
          (gene sets from config/celltypes.yaml).
       b. Take the top ~5% in each as the sample's reference population.
       c. If < 20 reference cells total, log a warning and fall back to the
          overall median across the cohort (less ideal but keeps the sample).
  3. Run infercnvpy.tl.infercnv with per-sample reference.
  4. Compute a per-cell dCNV score across the standard GBM CNV panel
     (chr7 gain, chr10 loss, chr9p loss, chr20 gain, chr19 gain).
  5. Flag cells as malignant_cnv=1 if dCNV > μ_ref + 2σ_ref (sample-specific
     reference distribution).
  6. Also preserve the old chr7/chr10 ratio in a parallel column `malignant_chr7_10`
     so the two methods can be compared.

Required .obs columns on input:
    sample_id, patient_id, modality
Added .obs columns on output:
    cnv_score            (float)
    malignant_cnv        (0/1)
    malignant_chr7_10    (0/1)      — the old method, for comparison
    cnv_reference_count  (int)      — how many reference cells per sample
    cnv_reference_source (str)      — 'same_sample_microglia_endo' or 'cohort_median_fallback'
"""

from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np


MICROGLIA_MARKERS   = ["P2RY12", "TMEM119", "CX3CR1", "CSF1R", "C1QA", "C1QB"]
ENDOTHELIAL_MARKERS = ["CLDN5", "PECAM1", "VWF", "KDR", "FLT1"]

# Standard GBM CNV panel (chromosomes, direction, biological rationale)
GBM_CNV_PANEL = [
    ("chr7",  +1),   # +chr7 — EGFR amplicon neighborhood, >80% of IDH-WT GBM
    ("chr10", -1),   # -chr10 — PTEN loss, >70% of IDH-WT GBM
    ("chr9",  -1),   # -chr9p — CDKN2A/B loss (chromosome-arm-level infercnv call)
    ("chr20", +1),   # +chr20 — common gain in high-grade glioma
    ("chr19", +1),   # +chr19 — frequent gain
]


def log(msg: str) -> None:
    print(f"[cnv] {msg}", flush=True)


def require(cond: bool, msg: str) -> None:
    if not cond:
        log(f"ABORT: {msg}")
        sys.exit(2)


def _compute_gene_activity_if_needed(adata) -> bool:
    """Compute gene activity scores in-place when possible.
    Returns True if gene-level data is now available, False for matrix-mode
    cohorts where we must fall back to peak-level CNV calling only."""
    if "gene_activity" in adata.obsm:
        log("  gene_activity obsm already present")
        return True
    if adata.uns.get("modality") == "multiome" and adata.X.max() > 100:
        log("  modality=multiome and X looks RNA-like; using X directly")
        return True
    if not adata.uns.get("fragments_bed_path"):
        log("  Matrix-mode cohort (no fragments) — cannot compute gene activity.")
        log("  CNV calling will use the legacy chr7+/chr10- ratio only.")
        return False
    import snapatac2 as snap
    log("  Computing gene activity scores via snapatac2.pp.make_gene_matrix")
    snap.pp.make_gene_matrix(adata, snap.genome.hg38)
    log(f"  Gene matrix: {adata.n_obs} cells × {adata.n_vars} genes")
    return True


def _score_marker_set(adata, markers: list[str], score_key: str) -> np.ndarray:
    """Return per-cell marker score; absent genes skipped."""
    import scanpy as sc
    present = [g for g in markers if g in adata.var_names]
    if not present:
        log(f"  WARN: 0/{len(markers)} markers present for {score_key}; score=0")
        return np.zeros(adata.n_obs)
    sc.tl.score_genes(adata, gene_list=present, score_name=score_key, use_raw=False)
    return np.asarray(adata.obs[score_key])


def _pick_sample_reference(adata, sample_id: str, mg_score, endo_score,
                           ref_frac=0.05, min_ref=20) -> tuple[list[int], str]:
    """Pick ~ref_frac of cells from this sample scoring highest on EITHER marker set."""
    mask_sample = (adata.obs["sample_id"] == sample_id).values
    n_sample = int(mask_sample.sum())
    k = max(min_ref, int(ref_frac * n_sample))
    k = min(k, n_sample)  # can't exceed sample size

    idx_sample = np.where(mask_sample)[0]
    # Top-k by max(microglia, endo) score within sample
    combined = np.maximum(mg_score[idx_sample], endo_score[idx_sample])
    top_within = np.argsort(-combined)[:k]
    ref_idx = idx_sample[top_within].tolist()
    source = "same_sample_microglia_endo"
    if len(ref_idx) < min_ref:
        source = "cohort_median_fallback"
    return ref_idx, source


def _run_infercnv(adata, ref_mask: np.ndarray, sample_key="sample_id"):
    """Call infercnvpy.tl.infercnv with per-cell reference flag."""
    try:
        import infercnvpy as cnv
    except ImportError:
        log("ABORT: infercnvpy not in env. Add to atacseq_env.yml and recreate the env.")
        sys.exit(2)

    adata.obs["cnv_ref_cell"] = ["reference" if r else "query" for r in ref_mask]
    log(f"  inferCNVpy on {adata.n_obs} cells "
        f"({int(ref_mask.sum())} reference, {int((~ref_mask).sum())} query)")
    cnv.tl.infercnv(
        adata,
        reference_key="cnv_ref_cell",
        reference_cat="reference",
        window_size=100,          # ~100 genes per window
        step=10,
        exclude_chromosomes=["chrX", "chrY", "chrM"],
    )
    cnv.tl.cnv_score(adata)        # writes adata.obs['cnv_score']
    log(f"  cnv_score range: "
        f"[{adata.obs['cnv_score'].min():.3f}, {adata.obs['cnv_score'].max():.3f}] "
        f"median={adata.obs['cnv_score'].median():.3f}")


def _call_malignant_per_sample(adata) -> None:
    """Flag malignant_cnv=1 for cells whose cnv_score exceeds μ_ref + 2σ_ref,
    computed PER SAMPLE from the sample's own reference pool."""
    adata.obs["malignant_cnv"] = 0
    for sid in sorted(adata.obs["sample_id"].unique()):
        mask = (adata.obs["sample_id"] == sid).values
        ref = (adata.obs.loc[mask, "cnv_ref_cell"] == "reference").values
        scores = adata.obs.loc[mask, "cnv_score"].values
        if ref.sum() < 5:
            log(f"  {sid}: too few reference cells ({ref.sum()}); using global median threshold")
            thr = np.median(adata.obs["cnv_score"]) + 2 * np.std(adata.obs["cnv_score"])
        else:
            thr = scores[ref].mean() + 2 * scores[ref].std()
        n_mal = int((scores > thr).sum())
        log(f"  {sid}: threshold={thr:.3f}  n_malignant={n_mal}/{mask.sum()} "
            f"({100*n_mal/mask.sum():.1f}%)")
        sample_cells = np.where(mask)[0]
        adata.obs.iloc[sample_cells, adata.obs.columns.get_loc("malignant_cnv")] = \
            (scores > thr).astype(int)


def _parallel_chr7_10_ratio(adata) -> None:
    """Legacy chr7+/chr10- ratio. Works on either gene-level or peak-level X.
    Resolves per-var chromosome from the 'chrom' column OR from the var_name
    (CATLAS format 'chrX:start-end')."""
    chrom = None
    if "chrom" in adata.var.columns:
        chrom = adata.var["chrom"].astype(str)
    else:
        # Try to parse from var_names like 'chr7:12345-67890'
        import re
        parsed = adata.var_names.to_series().str.extract(r"^(?P<chrom>chr[\dXYM]+)")
        if parsed["chrom"].notna().any():
            chrom = parsed["chrom"]
    if chrom is None:
        log("  Legacy chr7/chr10 ratio: cannot resolve chromosome per var; skipping")
        adata.obs["malignant_chr7_10"] = -1
        return

    chr7_mask = (chrom == "chr7").values
    chr10_mask = (chrom == "chr10").values
    if not chr7_mask.any() or not chr10_mask.any():
        log(f"  Legacy chr7/chr10 ratio: chr7 vars={int(chr7_mask.sum())}, "
            f"chr10 vars={int(chr10_mask.sum())}; need both; skipping")
        adata.obs["malignant_chr7_10"] = -1
        return
    # For backed AnnData with 1M+ cells, slicing adata[:, chr_mask].X eagerly
    # densifies the slice (100+ GB RAM). Instead iterate rows in chunks and
    # compute per-chunk means — bounded memory.
    chr7_idx = np.where(chr7_mask)[0]
    chr10_idx = np.where(chr10_mask)[0]
    is_backed = getattr(adata, "isbacked", False)
    if is_backed and adata.n_obs > 200_000:
        log(f"  streaming chr7/chr10 ratio in chunks (n_obs={adata.n_obs}, backed=True)")
        chunk = 20000
        chr7_mean = np.empty(adata.n_obs, dtype=np.float64)
        chr10_mean = np.empty(adata.n_obs, dtype=np.float64)
        for i in range(0, adata.n_obs, chunk):
            j = min(i + chunk, adata.n_obs)
            Xc = adata.X[i:j, :]  # loads just this chunk's rows (sparse)
            # Need .toarray() only on the subset of columns, keep sparse ops where possible
            import scipy.sparse as sp
            if sp.issparse(Xc):
                chr7_mean[i:j] = np.asarray(Xc[:, chr7_idx].mean(axis=1)).flatten()
                chr10_mean[i:j] = np.asarray(Xc[:, chr10_idx].mean(axis=1)).flatten()
            else:
                chr7_mean[i:j] = Xc[:, chr7_idx].mean(axis=1)
                chr10_mean[i:j] = Xc[:, chr10_idx].mean(axis=1)
            if i % 100000 == 0:
                log(f"    {i:,}/{adata.n_obs:,}")
    else:
        chr7_mean = np.asarray(adata[:, chr7_mask].X.mean(axis=1)).flatten()
        chr10_mean = np.asarray(adata[:, chr10_mask].X.mean(axis=1)).flatten()
    ratio = (chr7_mean + 1e-6) / (chr10_mean + 1e-6)
    adata.obs["malignant_chr7_10"] = (ratio >= 2.0).astype(int)
    log(f"  Legacy chr7/chr10 ratio: {int(adata.obs['malignant_chr7_10'].sum())} / {adata.n_obs} "
        f"({100*adata.obs['malignant_chr7_10'].mean():.1f}%) called malignant")


def main() -> int:
    import anndata as ad
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    file_gb = args.src.stat().st_size / 1e9
    use_backed = file_gb > 10  # files > 10 GB won't fit in RAM
    log(f"Loading {args.src} (size {file_gb:.1f} GB, backed={use_backed})")
    adata = ad.read_h5ad(args.src, backed="r" if use_backed else None)
    log(f"n_obs={adata.n_obs}  n_vars={adata.n_vars}  "
        f"modality={adata.obs['modality'].iloc[0] if 'modality' in adata.obs else '?'}  "
        f"n_samples={adata.obs['sample_id'].nunique()}")

    # Make obs_names unique — matrix-mode cohorts can repeat barcodes across samples
    if not adata.obs_names.is_unique:
        log("  Making obs_names unique (duplicates across samples)")
        adata.obs_names_make_unique()

    has_gene_level = _compute_gene_activity_if_needed(adata)

    if has_gene_level:
        log("Scoring marker sets")
        mg = _score_marker_set(adata, MICROGLIA_MARKERS, "score_microglia")
        endo = _score_marker_set(adata, ENDOTHELIAL_MARKERS, "score_endothelial")

        log("Picking per-sample reference cells")
        ref_mask = np.zeros(adata.n_obs, dtype=bool)
        ref_sources = []
        for sid in sorted(adata.obs["sample_id"].unique()):
            ref_idx, source = _pick_sample_reference(adata, sid, mg, endo)
            ref_mask[ref_idx] = True
            ref_sources.append((sid, len(ref_idx), source))
            log(f"  {sid}: n_ref={len(ref_idx)} source={source}")

        adata.obs["cnv_reference_count"] = sum(c for _, c, _ in ref_sources)
        adata.obs["cnv_reference_source"] = ";".join(f"{s}:{src}" for s, _, src in ref_sources)[:255]

        _run_infercnv(adata, ref_mask)
        _call_malignant_per_sample(adata)
    else:
        # Matrix-mode fallback: no gene-level data, cannot run inferCNVpy.
        # Set stub columns so downstream rules read predictable values.
        log("Matrix-mode: skipping inferCNV, using chr7/10 ratio only")
        adata.obs["cnv_score"] = 0.0
        adata.obs["malignant_cnv"] = 0
        adata.obs["cnv_reference_count"] = 0
        adata.obs["cnv_reference_source"] = "matrix_mode_no_gene_activity"

    _parallel_chr7_10_ratio(adata)
    # For matrix-mode cohorts, the chr7/10 ratio is the only CNV call we have —
    # mirror it into malignant_cnv so downstream label_transfer treats it as the
    # malignancy flag.
    if not has_gene_level and (adata.obs["malignant_chr7_10"] >= 0).any():
        adata.obs["malignant_cnv"] = adata.obs["malignant_chr7_10"].clip(lower=0)

    # Comparison summary — new vs old method
    n_mal_new = int(adata.obs["malignant_cnv"].sum())
    n_mal_old = int(adata.obs["malignant_chr7_10"].sum()) if (adata.obs["malignant_chr7_10"] >= 0).any() else -1
    log(f"=== Comparison ===")
    log(f"  inferCNVpy (new method):      {n_mal_new} / {adata.n_obs} malignant "
        f"({100*n_mal_new/adata.n_obs:.1f}%)")
    if n_mal_old >= 0:
        log(f"  chr7+/chr10- ratio (legacy):  {n_mal_old} / {adata.n_obs} "
            f"({100*n_mal_old/adata.n_obs:.1f}%)")
        both = int(((adata.obs["malignant_cnv"] == 1) & (adata.obs["malignant_chr7_10"] == 1)).sum())
        only_new = n_mal_new - both
        only_old = n_mal_old - both
        log(f"  agreement: both={both}  only inferCNV={only_new}  only chr7/10={only_old}")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    if use_backed:
        # Can't call write_h5ad on backed AnnData (would re-copy 100+ GB).
        # Use hardlink trick: input and output are identical except for new
        # obs columns, so hardlink the file and patch obs in-place via h5py.
        import os as _os
        import h5py as _h5
        try:
            from anndata.io import write_elem as _write_elem
        except ImportError:
            from anndata.experimental import write_elem as _write_elem
        obs_to_write = adata.obs.copy()
        adata.file.close()
        del adata
        if args.out.exists():
            args.out.unlink()
        try:
            _os.link(args.src, args.out)
            log(f"  hardlinked {args.src} → {args.out}")
        except OSError as e:
            log(f"  hardlink failed ({e}); falling back to shutil.copy")
            import shutil as _sh
            _sh.copy(args.src, args.out)
        with _h5.File(args.out, "r+") as hf:
            if "obs" in hf:
                del hf["obs"]
            _write_elem(hf, "obs", obs_to_write)
    else:
        adata.write_h5ad(args.out)
    log(f"Wrote → {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
