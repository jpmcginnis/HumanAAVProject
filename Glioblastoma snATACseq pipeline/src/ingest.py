"""
Standardize each cohort's deposit format to a shared AnnData schema.

Required obs columns on the output AnnData:
    - cohort            (str)   — e.g. "gbm_space"
    - patient_id        (str)   — standardized across cohorts
    - sample_id         (str)   — unique per sequencing library
    - region            (str)   — if multi-region sampling; else "unknown"
    - modality          (str)   — "snATAC" or "multiome"
    - idh_status        (str)   — "wildtype"
    - primary_recurrent (str)   — "primary", "recurrent", or "unknown"

Required var columns:
    - chrom, start, end — in GRCh38

Required layers:
    - counts (int32)
    - fragments_bed_path (uns)  — pointer to the per-cell fragments file for peak requantification

Design notes:
  - Fails loudly if the raw deposit doesn't match DEPOSIT_LAYOUTS.md.
  - Verbose per JP's conventions: prints n_cells, n_samples, per-sample row.
"""

from __future__ import annotations
import argparse
import gzip
import shutil
import sys
from pathlib import Path

import numpy as np


def log(msg: str) -> None:
    print(f"[ingest] {msg}", flush=True)


def require(cond: bool, msg: str) -> None:
    if not cond:
        log(f"ABORT: {msg}")
        sys.exit(2)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _import_fragments_set(frag_paths: list[Path], cohort: str,
                          out: Path, obs_setters: dict) -> None:
    """Import a list of fragments.tsv(.gz) via snapatac2 and write an AnnDataSet.

    obs_setters maps column name → callable(sample_id) → value.
    """
    import snapatac2 as snap
    import anndata as ad
    log(f"  Importing {len(frag_paths)} fragments files (disk-backed; no in-memory concat)")
    out.parent.mkdir(parents=True, exist_ok=True)

    # Write each imported fragments file to its OWN disk-backed .h5ad so we
    # never hold all samples in RAM at once. Then combine via snap.AnnDataSet
    # (path-based list) which just references them lazily.
    per_sample_paths = []
    for i, frag in enumerate(frag_paths, 1):
        sample_id = frag.name.replace(".fragments.tsv.gz", "").replace(".fragments.tsv", "")
        log(f"    [{i}/{len(frag_paths)}] {sample_id}")
        _import_fn = getattr(snap.pp, "import_fragments", None) or snap.pp.import_data
        per_sample_h5ad = out.parent / f"{sample_id}.h5ad"
        a = _import_fn(
            fragment_file=str(frag),
            chrom_sizes=snap.genome.hg38,
            sorted_by_barcode=False,
            file=str(per_sample_h5ad),   # disk-backed — avoids in-memory concat OOM
        )
        # snapatac2's backed AnnData uses polars; scalar assignment fails
        # ("length 0 vs height N"). Use full-length arrays instead.
        n = a.n_obs
        a.obs["cohort"] = np.full(n, cohort, dtype=object)
        a.obs["sample_id"] = np.full(n, sample_id, dtype=object)
        for col, fn in obs_setters.items():
            a.obs[col] = np.full(n, fn(sample_id), dtype=object)
        # Flush changes to disk and release the handle before next sample.
        try:
            a.close()
        except Exception:
            pass
        per_sample_paths.append((sample_id, per_sample_h5ad))
        log(f"      wrote {per_sample_h5ad}")

    # Merge the per-sample backed h5ads into one output. snap.AnnDataSet
    # accepts list[(name, Path)] in 2.9+.
    log(f"  Combining {len(per_sample_paths)} per-sample h5ads via snap.AnnDataSet")
    _ = snap.AnnDataSet(
        adatas=[(sid, str(p)) for sid, p in per_sample_paths],
        filename=str(out),
    )
    log(f"  → {out} (n_samples={len(per_sample_paths)})")


def _parse_mathewson_primrec(sample_id: str) -> str:
    low = sample_id.lower()
    if "recur" in low or re.search(r"_r(?:[0-9]|_|$)", low):
        return "recurrent"
    return "primary"


# ---------------------------------------------------------------------------
# Cohort ingestors
# ---------------------------------------------------------------------------

import re  # kept late to make imports cheap in the no-op dispatch path


def ingest_guilhamon(src: Path, out: Path) -> None:
    """GSE139136 — 4 primary IDHwt scATAC. Fragments inside GSE139136_RAW.tar
    (extracted by fetch.py). One patient per sample."""
    log(f"Guilhamon ingest from {src}")
    frags = sorted(src.rglob("*fragments.tsv.gz"))
    require(len(frags) >= 4,
            f"Expected ≥4 fragments files for Guilhamon, found {len(frags)}. "
            f"Check that GSE139136_RAW.tar was extracted in-place.")
    _import_fragments_set(
        frags, cohort="guilhamon", out=out,
        obs_setters={
            "patient_id": lambda s: s,          # 1:1 patient:sample
            "region": lambda s: "unknown",
            "modality": lambda s: "snATAC",
            "idh_status": lambda s: "wildtype",
            "primary_recurrent": lambda s: "primary",
        })


def ingest_gse276177_khan_fragments(src: Path, out: Path) -> None:
    """GSE276177 — Khan 2024 'Divergence from the human astrocyte developmental
    trajectory in glioblastoma'. Fragment-mode (NOT CellRanger matrix output as
    my earlier matrix-mode handler incorrectly assumed).

    Each GSM is deposited as a *_atac_fragments.tsv.gz + .tbi.gz tabix index.
    3 primary IDH-WT GBM patients, 5 snATAC samples:
      GSM8492625 GBM20_margin, GSM8492627 GBM20_tumor,
      GSM8492629 GBM25_margin, GSM8492631 GBM25_tumor,
      GSM8492633 GBM38_tumor.
    """
    log(f"Khan 2024 GSE276177 fragment-mode ingest from {src}")

    IDH_WT_ATAC_GSMS = {
        "GSM8492625": ("GBM20", "margin"),
        "GSM8492627": ("GBM20", "tumor"),
        "GSM8492629": ("GBM25", "margin"),
        "GSM8492631": ("GBM25", "tumor"),
        "GSM8492633": ("GBM38", "tumor"),
    }

    # Match our 5 IDH-WT GSMs' fragment files (ignore RNA h5 and other GSMs)
    frags = []
    for gsm, (patient, region) in IDH_WT_ATAC_GSMS.items():
        candidates = sorted(src.rglob(f"{gsm}_*atac_fragments.tsv.gz"))
        if candidates:
            frags.append((candidates[0], patient, region))
    require(len(frags) >= 3,
            f"GSE276177 expected ≥3 IDH-WT GBM fragments.tsv.gz, found {len(frags)}")
    log(f"  found {len(frags)} IDH-WT GBM fragment files")

    def patient_from_path(sample_id: str) -> str:
        for gsm, (pt, _) in IDH_WT_ATAC_GSMS.items():
            if gsm in sample_id:
                return pt
        return sample_id

    def region_from_path(sample_id: str) -> str:
        for gsm, (_, rg) in IDH_WT_ATAC_GSMS.items():
            if gsm in sample_id:
                return rg
        return "unknown"

    _import_fragments_set(
        [f for f, _, _ in frags], cohort="gse276177_khan_astro", out=out,
        obs_setters={
            "patient_id": patient_from_path,
            "region": region_from_path,
            "modality": lambda s: "snATAC",
            "idh_status": lambda s: "wildtype",
            "primary_recurrent": lambda s: "primary",
        })


def ingest_tcga_scatac(src: Path, out: Path) -> None:
    """TCGA scATAC-seq GBM subset — Terekhanova et al. 2024 Science
    (PMID 39236169), "Epigenetic regulation during cancer transitions
    across 11 tumour types".

    Open-access distribution via Zenodo (``full_cancer_scatacseq.tar.gz``,
    which bundles processed fragments for 11 cancer types). The GBM subset
    is 9 samples named ``scATAC_GBMx_<UUID>_X<N>_S<N>_B1_T1.fragments.tsv.gz``.
    All are stated as primary IDH-wildtype GBM in the paper's Data Availability
    section; we accept that on faith (no per-sample IDH table in the deposit).

    5 of 9 fragments are directly extractable from the Zenodo tar. The other 4
    require dbGaP phs000178 controlled access for the raw BAMs; running them
    yourself gives you the full 9-sample cohort.

    Each UUID segment is treated as its own patient_id. Fragment format is
    10x-standard (fragments.tsv.gz + .tbi); ``_import_fragments_set`` writes
    one per-sample snapatac2 h5ad each, then ``src/flatten_cohort.py`` turns
    the per-sample files into a flat CATLAS-space cohort h5ad.

    Previously slugged ``sundaram_gbm`` — the "Sundaram" release is actually
    the same TCGA scATAC data, so the cohort is now called ``tcga_scatac`` to
    reflect the true primary source.
    """
    log(f"TCGA scATAC GBM ingest from {src}")
    frags = sorted(f for f in src.rglob("scATAC_GBMx_*.fragments.tsv.gz") if f.is_file())
    require(len(frags) >= 2,
            f"Expected >=2 TCGA GBM fragments files, found {len(frags)}. "
            f"Verify the pan-cancer tar was extracted (look under "
            f"Cancer_scATACseq_data/).")
    log(f"  found {len(frags)} GBM scATAC samples")

    def patient_from_path(sample_id: str) -> str:
        # scATAC_GBMx_<UUID>_X<N>_S<N>... — use the UUID segment as patient_id
        m = re.search(r"scATAC_GBMx_([A-F0-9_]{36})", sample_id)
        return m.group(1) if m else sample_id

    _import_fragments_set(
        frags, cohort="tcga_scatac", out=out,
        obs_setters={
            "patient_id": patient_from_path,
            "region": lambda s: "unknown",
            "modality": lambda s: "snATAC",
            "idh_status": lambda s: "wildtype",
            "primary_recurrent": lambda s: "primary",
        })


def ingest_raviram(src: Path, out: Path) -> None:
    """GSE165037 — 3 primary IDHwt snATAC. Fragments are in per-GSM subdirs
    (one dir per GSM, each with a fragments.tsv.gz)."""
    log(f"Raviram ingest from {src}")
    frags = sorted(src.rglob("*fragments.tsv.gz"))
    require(len(frags) == 3,
            f"Expected exactly 3 IDH-WT snATAC fragments files for Raviram, "
            f"found {len(frags)}: {[str(f) for f in frags]}")
    # Parse GBM4 / GBM9 / GBM12 patient tag from the parent directory name
    def patient(sample_id: str) -> str:
        m = re.search(r"GBM\d+", sample_id)
        return m.group(0) if m else sample_id

    _import_fragments_set(
        frags, cohort="raviram", out=out,
        obs_setters={
            "patient_id": patient,
            "region": lambda s: "unknown",
            "modality": lambda s: "snATAC",
            "idh_status": lambda s: "wildtype",
            "primary_recurrent": lambda s: "primary",
        })


def ingest_tcga_scatac(src: Path, out: Path) -> None:
    """TCGA pan-cancer scATAC (Terekhanova 2024). Prefer the pre-integrated
    h5ad in S3; fall back to ArchR .arrow files if only those are present."""
    import anndata as ad
    log(f"TCGA scATAC ingest from {src}")
    pre = src / "integrated.h5ad"
    if pre.exists():
        log(f"  Using pre-integrated {pre.name} ({pre.stat().st_size/1e9:.1f} GB)")
        a = ad.read_h5ad(pre)
        a.obs["cohort"] = "tcga_scatac"
        a.obs["modality"] = "snATAC"
        a.obs["idh_status"] = "wildtype"
        a.obs["primary_recurrent"] = "primary"
        if "region" not in a.obs.columns:
            a.obs["region"] = "unknown"
        # patient_id/sample_id should already be set by the TCGA integration step
        require("patient_id" in a.obs.columns,
                "TCGA integrated.h5ad missing patient_id in .obs")
        log(f"  n_obs={a.n_obs}  n_patients={a.obs['patient_id'].nunique()}")
        out.parent.mkdir(parents=True, exist_ok=True)
        a.write_h5ad(out)
        log(f"TCGA → {out}")
        return
    log("ABORT: no integrated.h5ad found; add an ArchR → AnnData converter rule.")
    sys.exit(2)


# Wang/Mathewson 2022 Nat Cancer (GSE174554): verified IDH + stage per sample
# from the paper's Supplementary Table 1 (parsed 2026-10-05).
# SF4297 is IDH-mutant — DROP. SF9259R/S is in the paper but not in the GEO deposit.
WANG_SAMPLE_METADATA = {
    "SF10592": {"idh": "wildtype", "stage": "primary",   "pair": 5},
    "SF11857": {"idh": "wildtype", "stage": "recurrent", "pair": 5},
    "SF6621":  {"idh": "wildtype", "stage": "recurrent", "pair": 16},
    "SF6996":  {"idh": "wildtype", "stage": "primary",   "pair": 22},
    "SF7062":  {"idh": "wildtype", "stage": "recurrent", "pair": 22},
    "SF9510":  {"idh": "wildtype", "stage": "recurrent", "pair": 25},  # pair primary SF9259R/S not in GEO
    "SF9798":  {"idh": "wildtype", "stage": "primary",   "pair": 27},
    "SF9494":  {"idh": "wildtype", "stage": "recurrent", "pair": 27},
    "SF11331": {"idh": "wildtype", "stage": "primary",   "pair": 37},
    "SF4297":  {"idh": "mutant",   "stage": "primary",   "pair": None, "DROP": True},
}


def _wang_sample_from_name(fname: str) -> str | None:
    """Extract SF-sample ID from a fragments/matrix filename."""
    m = re.search(r"(SF\d+)", fname)
    return m.group(1) if m else None


def ingest_mathewson_lupien(src: Path, out: Path) -> None:
    """GSE174554 (Wang, Jung, Babikir 2022 Nat Cancer) — 10 deposited scATAC samples.
    IDH + stage come from the paper's Supplementary Table 1 (baked in above).
    Drops SF4297 (IDH-mutant)."""
    log(f"Wang/GSE174554 ingest from {src}")
    frags = sorted(src.rglob("*fragments.tsv.gz"))
    log(f"  Found {len(frags)} fragments files total")

    kept = []
    dropped = []
    for f in frags:
        sid = _wang_sample_from_name(f.name)
        if sid is None:
            dropped.append((f, "no SF id in filename"))
            continue
        meta = WANG_SAMPLE_METADATA.get(sid)
        if meta is None:
            dropped.append((f, f"unknown sample {sid} (not in Table 1)"))
            continue
        if meta.get("DROP"):
            dropped.append((f, f"{sid} IDH-mutant"))
            continue
        kept.append((f, sid, meta))
    log(f"  Keeping {len(kept)} IDH-WT scATAC fragments:")
    for f, sid, meta in kept:
        log(f"    {sid}  pair={meta['pair']}  stage={meta['stage']}  path={f.name}")
    for f, reason in dropped:
        log(f"    DROP: {f.name}  ({reason})")
    require(len(kept) >= 5,
            f"Wang expected ≥5 IDH-WT scATAC fragments, found {len(kept)}")

    def patient_for(sid: str) -> str:
        p = WANG_SAMPLE_METADATA[sid]["pair"]
        return f"Wang_pair_{p}" if p else sid

    def stage_for(sid: str) -> str:
        return WANG_SAMPLE_METADATA[sid]["stage"]

    # Use the per-sample resolution of our shared importer
    frag_paths = [f for f, _, _ in kept]
    name_to_sid = {f.name.replace(".fragments.tsv.gz", ""): sid for f, sid, _ in kept}
    _import_fragments_set(
        frag_paths, cohort="mathewson_lupien", out=out,
        obs_setters={
            "patient_id": lambda s: patient_for(name_to_sid.get(s, "")),
            "region":     lambda s: "unknown",
            "modality":   lambda s: "snATAC",
            "idh_status": lambda s: "wildtype",
            "primary_recurrent": lambda s: stage_for(name_to_sid.get(s, "")),
        })


def _gunzip_if_needed(src_path: Path) -> Path:
    """Return a path to an uncompressed .h5ad. Decompresses .h5ad.gz lazily
    (large files). Crash-safe: writes to <name>.part, fsyncs, then atomically
    renames to final name — avoids leaving a truncated file after a crash."""
    import os
    if src_path.suffix != ".gz":
        return src_path
    plain = src_path.with_suffix("")  # strips .gz
    if plain.exists() and plain.stat().st_size > 1e8:
        # Validate size matches .gz's uncompressed length to catch truncated files
        with gzip.open(src_path, "rb") as g:
            # Reading the ISIZE field is cheap but gzip only stores it modulo 2^32,
            # which fails for >4 GB streams. Instead do a quick HDF5 header check.
            pass
        # Try to open as HDF5 to detect truncation
        try:
            import h5py
            with h5py.File(plain, "r") as hf:
                _ = list(hf.keys())
            log(f"  Reusing existing valid {plain.name} ({plain.stat().st_size/1e9:.1f} GB)")
            return plain
        except Exception as e:
            log(f"  Existing {plain.name} is corrupt/truncated ({e}); re-gunzipping")
            plain.unlink()

    part = plain.with_suffix(plain.suffix + ".part")
    if part.exists():
        part.unlink()
    log(f"  Gunzipping {src_path.name} → {plain.name}  "
        f"({src_path.stat().st_size/1e9:.1f} GB input)")
    with gzip.open(src_path, "rb") as f_in, open(part, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out, length=16 * 1024 * 1024)
        f_out.flush()
        os.fsync(f_out.fileno())
    os.replace(part, plain)
    log(f"  Gunzip complete: {plain.name} ({plain.stat().st_size/1e9:.1f} GB)")
    return plain


def ingest_gbm_space(src: Path, out: Path) -> None:
    """E-MTAB-17183 — De Jong 2025 GBM-Space. We read the annotated ATAC h5ad
    directly (DEPOSIT_LAYOUTS.md confirms the file lives at
    GBM_space_ATAC_filtered_peaks.h5ad.gz)."""
    import anndata as ad
    import pandas as pd
    log(f"GBM-Space Multiome ingest from {src}")

    atac_gz = src / "GBM_space_ATAC_filtered_peaks.h5ad.gz"
    require(atac_gz.exists(),
            f"GBM_space_ATAC_filtered_peaks.h5ad.gz missing in {src} — "
            f"did fetch_gbm_space run?")
    atac_h5ad = _gunzip_if_needed(atac_gz)
    file_gb = atac_h5ad.stat().st_size / 1e9
    log(f"  Reading annotated ATAC h5ad ({file_gb:.1f} GB) — backed='r' for huge files")
    # KEY FIX from a prior hang: the 127 GB file caused anndata.read_h5ad to
    # cycle on reads. Use backed='r' here so X isn't materialized in memory;
    # obs/var still load into pandas DataFrames.
    a = ad.read_h5ad(atac_h5ad, backed="r")
    log(f"  n_obs={a.n_obs}  n_vars={a.n_vars}  X is backed (not in RAM)")

    # SDRF carries hash-filename → donor/region. Parse it so we can audit the
    # patient column against the deposit's own metadata table.
    sdrf = src / "E-MTAB-17183.sdrf.txt"
    sdrf_map: dict = {}
    if sdrf.exists():
        df = pd.read_csv(sdrf, sep="\t")
        log(f"  Loaded SDRF: {len(df)} rows, {df.shape[1]} cols")
        donor_cols = [c for c in df.columns if "individual" in c.lower() or "donor" in c.lower()]
        if donor_cols:
            log(f"  SDRF donor-ID column: {donor_cols[0]}")
            sdrf_map["donor_col"] = donor_cols[0]
    a.uns["sdrf_donor_col"] = sdrf_map.get("donor_col", "")

    # Harmonize to the shared schema
    a.obs["cohort"] = "gbm_space"
    a.obs["modality"] = "multiome"
    a.obs["idh_status"] = "wildtype"
    a.obs["primary_recurrent"] = "primary"
    for src_col, dst_col in [
        ("donor", "patient_id"),
        ("donor_id", "patient_id"),     # GBM-Space uses this
        ("individual", "patient_id"),
        ("tumour_region", "region"),
        ("region", "region"),
        ("site_id", "region"),          # GBM-Space uses this
        ("sample", "sample_id"),
        ("Sample", "sample_id"),        # GBM-Space uses this (capitalized)
        ("library", "sample_id"),
    ]:
        if src_col in a.obs.columns:
            if dst_col not in a.obs.columns or a.obs[dst_col].isna().all():
                a.obs[dst_col] = a.obs[src_col].astype(str)
    if "region" not in a.obs.columns:
        a.obs["region"] = "unknown"
    require("patient_id" in a.obs.columns,
            f"patient_id could not be set — inspect .obs columns: {list(a.obs.columns)[:30]}")
    log(f"  n_patients={a.obs['patient_id'].nunique()}  "
        f"n_regions={a.obs.get('region', pd.Series(['unknown']*a.n_obs)).nunique()}")

    out.parent.mkdir(parents=True, exist_ok=True)
    # CRITICAL: the earlier `.to_memory()` call blew up to 100+ GB RAM (OOM'd
    # the box). A plain shutil.copy to the output path would need another
    # ~127 GB of disk headroom (we don't have it).
    # Fix: HARDLINK the backed h5ad into the output location (zero disk cost
    # — same inode, two names), then patch obs/var in-place via h5py. The
    # hardlink works because input/output are on the same filesystem. After
    # the patch, the "output" is a valid standalone h5ad; the input raw file
    # can be deleted later to reclaim nothing (same inode) but the output
    # keeps the data intact.
    import os as _os
    import h5py as _h5
    try:
        from anndata.io import write_elem as _write_elem
    except ImportError:
        from anndata.experimental import write_elem as _write_elem

    try:
        a.file.close()
    except Exception:
        pass

    log(f"  hardlinking backed h5ad → output (zero-cost, same inode)")
    if out.exists():
        out.unlink()
    try:
        _os.link(atac_h5ad, out)
    except OSError as e:
        log(f"  hardlink failed ({e}); falling back to shutil.copy ({file_gb:.0f} GB)")
        import shutil as _sh
        _sh.copy(atac_h5ad, out)

    log("  patching obs/var in the hardlinked file in-place via h5py")
    with _h5.File(out, "r+") as hf:
        if "obs" in hf:
            del hf["obs"]
        _write_elem(hf, "obs", a.obs)
        if "var" in hf:
            del hf["var"]
        _write_elem(hf, "var", a.var)
        if a.uns:
            if "uns" in hf:
                del hf["uns"]
            _write_elem(hf, "uns", dict(a.uns))
    log(f"GBM-Space → {out}  ({out.stat().st_size/1e9:.1f} GB on disk)")


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def ingest_rds_cohort_stub(src: Path, out: Path) -> None:
    """Stub for cohorts whose deposited data is R-serialized (ArchR, Signac, SnapATAC).

    Detects .rds files and emits a clear error directing to the R extraction
    workstream — DO NOT silently skip these cohorts. They need to be processed
    through src/extract_rds.R first (which converts to h5ad), then re-run
    through the normal ingest path.
    """
    rds_files = sorted(src.rglob("*.rds")) + sorted(src.rglob("*.rds.gz"))
    log(f"R-serialized data detected: {len(rds_files)} .rds files")
    for p in rds_files[:5]:
        log(f"  {p}  ({p.stat().st_size/1e9:.1f} GB)")
    log("ABORT: this cohort needs R-side extraction before Python ingest.")
    log("  Run: Rscript src/extract_rds.R --rds <input.rds> --out <output.h5ad>")
    log("  Then re-run this rule with the resulting h5ad as --in.")
    sys.exit(2)


HANDLERS = {
    "guilhamon": ingest_guilhamon,
    "raviram": ingest_raviram,
    "gse165037": ingest_rds_cohort_stub,          # R SnapATAC .rds.gz
    "tcga_scatac": ingest_tcga_scatac,            # TCGA scATAC GBM (ex-"sundaram_gbm")
    "mathewson_lupien": ingest_mathewson_lupien,
    "gbm_space": ingest_gbm_space,
    "gbm_tme_atlas_hra004942": ingest_rds_cohort_stub,   # Signac/ArchR .rds
    "gse276177_khan_astro": ingest_gse276177_khan_fragments,  # fragment-mode (not matrix)
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(HANDLERS))
    ap.add_argument("--in", dest="src", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    log(f"=== Ingesting {args.cohort} ===")
    HANDLERS[args.cohort](args.src, args.out)
    log(f"=== Ingest complete: {args.cohort} ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
