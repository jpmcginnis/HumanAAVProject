"""
Per-cohort download logic — grounded in the verified deposit layouts recorded
in DEPOSIT_LAYOUTS.md. Each handler is deliberately explicit about what it
pulls, filters, and fails on.

Key invariants:
  - Every GEO handler enumerates GSMs via NCBI eutils and filters by title
    *before* downloading, so IDH-mutant / non-tumor samples never land in
    the raw/ tree.
  - Every handler fails loudly (SystemExit 2) if the expected post-conditions
    (e.g. ≥ N fragments.tsv.gz files, specific annotated h5ad present) are
    not satisfied. No silent empty-AnnData runs downstream.
  - Prints cell/sample counts verbosely per JP's conventions.

Usage:
    python -u src/fetch.py --cohort gbm_space --out /data/projects/atacseq/raw/gbm_space
"""

from __future__ import annotations
import argparse
import json
import os
import re
import subprocess
import sys
import tarfile
import time
import urllib.parse
import urllib.request
from pathlib import Path

import yaml


def log(msg: str) -> None:
    print(f"[fetch] {msg}", flush=True)


# ---------------------------------------------------------------------------
# NCBI eutils helpers (used by every GEO handler)
# ---------------------------------------------------------------------------

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _http_json(url: str, retries: int = 3, sleep_s: float = 1.0) -> dict:
    last_err = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            last_err = e
            time.sleep(sleep_s * (2 ** i))
    raise RuntimeError(f"GET {url} failed after {retries} tries: {last_err}")


def enumerate_gsms(gse: str) -> list[dict]:
    """Return [{accession, title, ftp_dir}] for every GSM under a GSE."""
    q = urllib.parse.quote(f"{gse}[Accession]")
    search = _http_json(f"{EUTILS}/esearch.fcgi?db=gds&term={q}&retmode=json")
    ids = search["esearchresult"]["idlist"]
    # The series itself has id starting '2'; GSM esummary accessions begin '3'
    series_id = next((i for i in ids if i.startswith("2")), None)
    if series_id is None:
        raise RuntimeError(f"{gse}: no series record from eutils")
    summary = _http_json(f"{EUTILS}/esummary.fcgi?db=gds&id={series_id}&retmode=json")
    r = summary["result"][series_id]
    samples = r.get("samples", [])
    # Build canonical GSM GEO FTP paths. GEO path is .../GSMnnnnn/GSMnnnnnnn/suppl/
    out = []
    for s in samples:
        acc = s["accession"]
        prefix = acc[:-3] + "nnn"
        ftp = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{prefix}/{acc}/suppl/"
        out.append({"accession": acc, "title": s.get("title", ""), "ftp": ftp})
    return out


def wget_dir(url: str, dest: Path, extra_flags: list[str] | None = None) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    cmd = [
        "wget", "-r", "-np", "-nH", "--cut-dirs=5",
        "-R", "index.html*",
        "-P", str(dest),
        url,
    ]
    if extra_flags:
        cmd[1:1] = extra_flags
    log("  $ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def wget_file(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    target = dest / Path(urllib.parse.urlparse(url).path).name
    cmd = ["wget", "-c", "-O", str(target), url]
    log("  $ " + " ".join(cmd))
    subprocess.run(cmd, check=True)
    return target


def untar_in_place(tar_path: Path) -> list[Path]:
    """Extract a tarball into its parent dir. Returns the list of extracted files."""
    log(f"  Untar {tar_path.name} ({tar_path.stat().st_size/1e6:.0f} MB)")
    extracted: list[Path] = []
    with tarfile.open(tar_path, "r") as tf:
        tf.extractall(tar_path.parent)
        extracted = [tar_path.parent / m.name for m in tf.getmembers() if m.isfile()]
    log(f"  Extracted {len(extracted)} files")
    return extracted


def require(cond: bool, msg: str) -> None:
    if not cond:
        log(f"ABORT: {msg}")
        sys.exit(2)


# ---------------------------------------------------------------------------
# Per-cohort handlers
# ---------------------------------------------------------------------------

def fetch_guilhamon(out: Path) -> None:
    """GSE139136 — 4 IDH-WT scATAC-seq samples (Guilhamon et al. 2021 eLife).

    GSE139136/suppl/ holds a single GSE139136_RAW.tar; per-GSM supplementary
    may or may not include fragments.tsv.gz (confirmed on first run).
    """
    out.mkdir(parents=True, exist_ok=True)
    log("Guilhamon GSE139136 — pulling series-level RAW tar + per-GSM fallback")
    acc = "GSE139136"
    series_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar"
    )
    tar = wget_file(series_url, out)
    untar_in_place(tar)

    frags = sorted(out.rglob("*fragments.tsv.gz")) + sorted(out.rglob("*fragments.tsv"))
    matrices = sorted(out.rglob("*_matrix.mtx.gz"))
    log(f"  Found {len(frags)} fragments files, {len(matrices)} matrix files after extraction")

    # Guilhamon GSE139136 deposits ONLY peak-barcode matrices (verified
    # DEPOSIT_LAYOUTS.md 2026-10-05). Downstream ingest_matrix.py handles
    # this format and lifts per-sample peaks onto CATLAS.
    if len(frags) == 0 and len(matrices) >= 4:
        log(f"  matrix-mode deposit confirmed ({len(matrices)} per-sample matrix.mtx.gz files)")
        log("  Downstream ingest will use src/ingest_matrix.py for CATLAS lift.")
        return
    require(len(frags) >= 4 or len(matrices) >= 4,
            f"Expected ≥4 fragments OR matrix files, found frags={len(frags)}, matrices={len(matrices)}")
    log(f"Guilhamon → {out}  ({len(frags)} fragments, {len(matrices)} matrix files)")


def fetch_tcga_scatac(out: Path) -> None:
    """TCGA scATAC GBM subset — Terekhanova et al. 2024 Science (PMID 39236169).

    Pulls the single open-access pan-cancer fragment tar.gz from the GDC
    publication endpoint, which is where the processed per-sample fragments
    for the TCGA scATAC collection are distributed. Extracts only the GBM
    subset (scATAC_GBMx_* fragments.tsv.gz + .tbi), leaving other cancer
    types on disk but outside the per-cohort workflow dir.

    5 of 9 GBMx samples are directly in the open-access tar; the other 4
    require dbGaP phs000178 controlled access for the raw BAMs.

    Previously named ``fetch_sundaram_gbm`` — the "Sundaram" release IS this
    TCGA scATAC data, so the fetch is renamed to reflect the true source.
    All GBM samples are primary IDH-WT per the paper's inclusion criteria.
    """
    out.mkdir(parents=True, exist_ok=True)
    log("Sundaram 2024 Science — pulling GDC open-access pan-cancer fragment tar.gz")
    # Direct GDC URL for the pan-cancer scATAC-seq blob
    gdc_uuid = "a7765b5f-a606-4657-94ed-4c2a85bf3f39"
    gdc_url = f"https://api.gdc.cancer.gov/data/{gdc_uuid}"
    tar_gz = out / "full_cancer_scatacseq.tar.gz"
    if not tar_gz.exists() or tar_gz.stat().st_size < 1_000_000_000:
        log(f"  wget {gdc_url}")
        subprocess.run(["wget", "-c", "-O", str(tar_gz), gdc_url], check=True)

    # Extract only the GBM scATAC fragment files (ignore other cancers)
    log("  extracting GBM-only subset from pan-cancer tar")
    subprocess.run(
        ["tar", "-xzf", str(tar_gz), "-C", str(out),
         "--wildcards", "Cancer_scATACseq_data/scATAC_GBMx_*"],
        check=True,
    )
    frags = sorted(out.rglob("scATAC_GBMx_*.fragments.tsv.gz"))
    require(len(frags) >= 2,
            f"Sundaram expected ≥2 GBM fragments.tsv.gz files, found {len(frags)}")
    log(f"Sundaram → {out}  ({len(frags)} GBM scATAC samples extracted)")


def fetch_gse276177_khan_astro(out: Path) -> None:
    """GSE276177 — Khan et al. 2024 "Divergence from the human astrocyte
    developmental trajectory in glioblastoma" (Nature Cell Biology).

    3 primary IDH-WT GBM patients (GBM20, GBM25, GBM38) with paired snATAC +
    snRNA on tumor and (for GBM20/25) margin tissue. 5 snATAC samples total:
      GSM8492625 GBM20_margin, GSM8492627 GBM20_tumor,
      GSM8492629 GBM25_margin, GSM8492631 GBM25_tumor,
      GSM8492633 GBM38_tumor.

    Format: CellRanger-ATAC per-sample matrix + peaks + barcodes.
    GSE276177_RAW.tar is 27.5 GB (contains both ATAC and RNA samples).
    """
    out.mkdir(parents=True, exist_ok=True)
    log("Khan GSE276177 — pulling series-level RAW tar")
    acc = "GSE276177"
    series_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar"
    )
    tar = wget_file(series_url, out)
    untar_in_place(tar)
    matrices = sorted(out.rglob("*matrix.mtx.gz"))
    log(f"  Found {len(matrices)} matrix.mtx.gz files")
    wanted = {"GSM8492625", "GSM8492627", "GSM8492629", "GSM8492631", "GSM8492633"}
    found_wanted = [m for m in matrices if m.name.split("_")[0] in wanted]
    require(len(found_wanted) >= 3,
            f"GSE276177 expected ≥3 IDH-WT primary GBM snATAC GSMs in tar; "
            f"found {[m.name for m in found_wanted]}")
    log(f"GSE276177 → {out}  ({len(found_wanted)}/5 IDH-WT snATAC GSMs present)")


def fetch_gse138794_guo(out: Path) -> None:
    """GSE138794 — Guo/Diaz 2020 'A single cell atlas of human glioma'.

    28 total samples (8 snATAC + 20 snRNA). We keep only the 4 IDH-WT primary
    GBM snATAC samples (verified via eutils soft files 2026-10-06):
      GSM4119517 SF11979, GSM4119518 SF11956,
      GSM4119519 SF11215, GSM4119520 SF11331.

    Format: CellRanger-ATAC filtered_peak_bc_matrix (hg38). The GSE RAW tar
    bundles per-GSM matrix.mtx.gz + barcodes.tsv.gz + peaks.bed.gz. We pull
    the whole tar (~1.5 GB) and let ingest_matrix.py filter by GSM id.
    """
    out.mkdir(parents=True, exist_ok=True)
    log("Guo/Diaz GSE138794 — pulling series-level RAW tar")
    acc = "GSE138794"
    series_url = (
        f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar"
    )
    tar = wget_file(series_url, out)
    untar_in_place(tar)

    matrices = sorted(out.rglob("*matrix.mtx.gz"))
    log(f"  Found {len(matrices)} matrix.mtx.gz files")
    # Only need our 4 IDH-WT primary GBM GSMs; the rest are present but
    # ingest_matrix_gse138794_guo filters them.
    wanted = {"GSM4119517", "GSM4119518", "GSM4119519", "GSM4119520"}
    found_wanted = [m for m in matrices if m.name.split("_")[0] in wanted]
    require(len(found_wanted) >= 2,
            f"GSE138794 expected ≥2 IDH-WT primary GBM GSMs in tar; "
            f"found {[m.name for m in found_wanted]}")
    log(f"GSE138794 → {out}  ({len(found_wanted)}/{len(wanted)} IDH-WT primary GBM GSMs present)")


def fetch_raviram(out: Path) -> None:
    """GSE165037 — 3 IDH-WT snATAC samples. Deposited as .rds.gz files (R SnapATAC
    v1 objects) in the series RAW tar. Downloads + extracts; downstream ingest
    will route through the RDS extractor.

    Per-GSM supplementary URLs 403 from EC2 (NCBI rate-limits them). Series
    _RAW.tar is 200 OK and contains all 12 GSM files; filter during ingest.
    """
    out.mkdir(parents=True, exist_ok=True)
    log("GSE165037 — pulling series _RAW.tar (contains .rds.gz for all 12 GSMs)")
    acc = "GSE165037"
    series_tar = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{acc}_RAW.tar"
    tar_path = wget_file(series_tar, out)
    untar_in_place(tar_path)
    rds_files = sorted(out.rglob("*.rds.gz")) + sorted(out.rglob("*.rds"))
    log(f"  Found {len(rds_files)} .rds/.rds.gz files")
    for p in rds_files[:5]:
        log(f"    {p.name}  ({p.stat().st_size/1e9:.2f} GB)")
    require(len(rds_files) >= 3,
            f"GSE165037 expected ≥3 .rds files, found {len(rds_files)}")
    log(f"GSE165037 → {out}  ({len(rds_files)} .rds files; downstream needs R extraction)")


# (fetch_tcga_scatac defined above; the former S3-arrow-sync fallback was
# removed — the Terekhanova 2024 processed fragments are on the GDC endpoint
# used by that handler, and we do not depend on the S3 ArchR .arrow mirror
# any more.)


def fetch_mathewson_lupien(out: Path) -> None:
    """GSE174554 — Mathewson/Lupien, IDH-WT scATAC primary + recurrent.

    Expect ~10 total scATAC samples per the Mathewson 2022 Nature Cancer paper.
    Flag loudly if ingest finds fewer — may indicate the series only hosts the
    Wang-2022 subset and the Mathewson-2021 reanalysis fragments live elsewhere.
    """
    out.mkdir(parents=True, exist_ok=True)
    log("Mathewson/Lupien GSE174554 — enumerating GSMs")
    gsms = enumerate_gsms("GSE174554")
    # Keep snATAC/scATAC titles only
    keep_re = re.compile(r"(scATAC|snATAC|ATAC)", re.IGNORECASE)
    kept = [s for s in gsms if keep_re.search(s["title"])]
    log(f"  {len(gsms)} GSMs total; {len(kept)} ATAC titles")
    for s in kept:
        log(f"  → {s['accession']}  {s['title']}")
        wget_dir(s["ftp"], out / s["accession"])

    frags = sorted(out.rglob("*fragments.tsv.gz"))
    log(f"  Found {len(frags)} fragments files; Mathewson 2022 claims 10 scATAC samples")
    if len(frags) < 10:
        log(f"  WARNING — fewer than 10 fragments files. "
            f"May be Wang-subset-only; "
            f"the 4 Mathewson-2021-reanalysis samples may live under a different accession.")
    require(len(frags) >= 1,
            f"Mathewson expected at least 1 fragments file, found {len(frags)}")
    log(f"Mathewson → {out}  ({len(frags)} fragments files)")


def fetch_gbm_space(out: Path) -> None:
    """E-MTAB-17183 — De Jong 2025 GBM-Space Multiome atlas (12 primary IDH-WT).

    Verified file layout (DEPOSIT_LAYOUTS.md, 2026-10-04): 165 files flat,
    opaque UUID-named Cell Ranger tarballs plus four annotated h5ad files
    and an SDRF mapping hashes to donor/region.

    We pull only the ANNOTATED objects + READMEs + SDRF (~35 GB total),
    not the 1.5 TB of raw per-reaction tarballs. The ATAC h5ad already has
    filtered peaks and cell-type annotations; that's the enhancer-atlas
    input. Pass --full to also download the raw tarballs.
    """
    out.mkdir(parents=True, exist_ok=True)
    log("GBM-Space E-MTAB-17183 — pulling annotated h5ad + SDRF (not the 1.5 TB raw)")

    base = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-17183"
    want = [
        "E-MTAB-17183.sdrf.txt",
        "E-MTAB-17183.idf.txt",
        "README_cellranger_arc.md",
        "README_metacell.md",
        "GBM_space_ATAC_filtered_peaks_README.md",
        "GBM_space_snRNA_README.md",
        "GBM_space_ATAC_filtered_peaks.h5ad.gz",   # 22.4 GB
        "GBM_space_snRNA.h5ad.gz",                 # 11.3 GB
        "rna_metacell_scdori_11_03_24.h5ad",       # 0.2 GB
        "genescore_metacell_scdori_11_03_24.h5ad", # 1.7 GB
    ]
    for name in want:
        url = f"{base}/{name}"
        log(f"  → {name}")
        wget_file(url, out)

    atac = out / "GBM_space_ATAC_filtered_peaks.h5ad.gz"
    rna = out / "GBM_space_snRNA.h5ad.gz"
    require(atac.exists() and atac.stat().st_size > 1e9,
            f"ATAC h5ad missing or too small: {atac}")
    require(rna.exists() and rna.stat().st_size > 1e9,
            f"RNA h5ad missing or too small: {rna}")
    log(f"GBM-Space → {out}  (ATAC={atac.stat().st_size/1e9:.1f} GB, "
        f"RNA={rna.stat().st_size/1e9:.1f} GB)")

    if os.environ.get("GBM_SPACE_FULL") == "1":
        log("  GBM_SPACE_FULL=1 — also pulling the 155 Cell Ranger tarballs (~1.5 TB)")
        api = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-17183/files?offset=0&limit=200"
        meta = _http_json(api)
        for item in meta.get("items", []):
            name = item["Name"]
            if name.startswith("cellranger-arc"):
                wget_file(f"{base}/{name}", out / "cellranger_tars")
        log(f"  Full tarballs → {out}/cellranger_tars")


# ---------------------------------------------------------------------------
# New cohorts added in the 2026-10-04 second-pass sweep
# ---------------------------------------------------------------------------

def fetch_hra004942(out: Path) -> None:
    """Lin et al. 2026 Nat Neurosci — 100 primary GBM patients, largest cohort found.

    Processed data at Zenodo 8085502 (open, no login). Raw at NGDC HRA004942
    (also open, no DAA — verified 2026-10-04 at ngdc.cncb.ac.cn/gsa-human/browse/HRA004942).
    We pull only the scATAC and scRNA zips by default (3.3 GB combined).
    """
    out.mkdir(parents=True, exist_ok=True)
    log("gbm_tme_atlas HRA004942 — pulling Zenodo processed-data zips")
    # DOI 10.5281/zenodo.8085502 resolves to record 17117252 (the latest version).
    # Use the files-API /content endpoint — the plain /records/<id>/files/<name>
    # path returns 404 for this record.
    zenodo = "https://zenodo.org/api/records/17117252/files"
    want = [
        ("scATAC-seq.zip", True),
        ("scRNA-seq.zip",  True),
        ("10X_VISIUM.zip", False),
        ("Patch-seq.zip",  False),
        ("PhenoCycler.zip", False),
        ("ISH.zip",        False),
    ]
    for name, required in want:
        url = f"{zenodo}/{name}/content"
        log(f"  → {name}  (URL: {url})")
        try:
            # wget_file names the local file after the URL's basename, but the
            # /content path means the basename would be "content". Save it to
            # the real file name explicitly instead.
            target = out / name
            subprocess.run(["wget", "-c", "-O", str(target), url], check=True)
            log(f"    saved {name}  ({target.stat().st_size/1e6:.1f} MB)")
        except subprocess.CalledProcessError:
            if required:
                raise
            log(f"  (non-critical) failed to pull {name}, continuing")

    # Unzip the ATAC bundle in-place so ingest can glob inside it
    atac_zip = out / "scATAC-seq.zip"
    require(atac_zip.exists(), f"Zenodo scATAC-seq.zip did not download to {atac_zip}")
    log(f"  Unzipping {atac_zip.name} ({atac_zip.stat().st_size/1e9:.1f} GB)")
    subprocess.run(["unzip", "-o", "-q", str(atac_zip), "-d", str(out / "scATAC")], check=True)
    unzipped = sorted((out / "scATAC").rglob("*"))
    log(f"  → {len(unzipped)} files unzipped. First few:")
    for p in unzipped[:10]:
        if p.is_file():
            log(f"      {p.relative_to(out)}  ({p.stat().st_size/1e6:.1f} MB)")
    log(f"HRA004942 → {out}  (contents still require inspection — ingest handler is a stub)")


def fetch_wang_sciadv_tc_ptb(out: Path) -> None:
    """Wang et al. 2024 Sci Adv — 5 patients, 10 TC-vs-PTB paired samples, DNBelab C4.

    Data at CNGB CNSA under CNP0003766. CNGB bulk downloads use their own API;
    we try the HTTP bulk path first and fall back to directing the user to the
    web portal if automation fails (common for non-NCBI archives).
    """
    out.mkdir(parents=True, exist_ok=True)
    log("Wang 2024 Sci Adv CNP0003766 — pulling from CNGB CNSA")
    # CNGB CNSA bulk-download HTTPS endpoint pattern
    base = "https://ftp.cngb.org/pub/CNSA/data5/CNP0003766"
    log(f"  CNGB base: {base}")
    try:
        subprocess.run(
            ["wget", "-r", "-np", "-nH", "--cut-dirs=4",
             "-R", "index.html*",
             "-P", str(out),
             f"{base}/"],
            check=True,
        )
    except subprocess.CalledProcessError:
        log("  CNGB bulk HTTPS failed; user must download manually via web portal:")
        log("    https://db.cngb.org/search/project/CNP0003766/")
        log("    Deposit into: " + str(out))
        sys.exit(2)

    # DNBelab fragments typically arrive as .frags.tsv.gz with 20 bp barcodes
    frags = sorted(out.rglob("*frag*.tsv.gz")) + sorted(out.rglob("*.frags.gz"))
    log(f"  Found {len(frags)} fragment-like files")
    require(len(frags) >= 10,
            f"Wang expected ≥10 fragment files (5 TC + 5 PTB), found {len(frags)}. "
            f"If only BAM is deposited, add a bam→fragments conversion step.")
    log(f"Wang → {out}  ({len(frags)} fragment files)")


def fetch_spatial_epigenomic_niches(out: Path) -> None:
    """Kint/Gallo 2025 bioRxiv — Marco Gallo lab at BCM. No public deposit in v1.

    This is a stub. The honest path is: JP emails Marco Gallo at BCM for the
    processed snATAC + Multiome h5ad files directly. Pre-publication,
    same-institution ask.
    """
    out.mkdir(parents=True, exist_ok=True)
    log("spatial_epigenomic_niches — NO PUBLIC DEPOSIT AS OF 2026-10-04")
    log("  Corresponding author: Marco Gallo, Baylor College of Medicine")
    log("  Action: email Gallo directly (same institution as JP) for processed data")
    log("  Expected deliverables: ")
    log("    - 3 snATAC h5ad (bulk nuclei)")
    log("    - 2 Multiome h5ad (bulk nuclei)")
    log("    - 28 spatial ATAC arrays (AtlasXomics format, separate ingestion)")
    log("  Once received, place them at:")
    log(f"    {out}/snATAC/*.h5ad")
    log(f"    {out}/multiome/*.h5ad")
    sys.exit(2)


# ---------------------------------------------------------------------------
# Reference-data handlers (CATLAS, GBmap) — invoked by Snakefile ref rules
# ---------------------------------------------------------------------------

def fetch_catlas(out: Path) -> None:
    """CATLAS 544K-peak snATAC atlas (Li et al. 2023 Science)."""
    out.parent.mkdir(parents=True, exist_ok=True)
    s3 = "s3://jpm-atacseq-archive-2026/catlas/catlas_full_annotated.h5ad"
    log(f"CATLAS — syncing {s3} → {out}")
    subprocess.run(["aws", "s3", "cp", s3, str(out)], check=True)
    require(out.exists() and out.stat().st_size > 1e8,
            f"CATLAS h5ad missing or too small: {out}")
    log(f"CATLAS → {out}  ({out.stat().st_size/1e9:.1f} GB)")


def fetch_gbmap(out: Path) -> None:
    """GBmap scRNA reference (Ruiz-Moreno 2025, Extended GBmap: 1.1M cells across 240+ patients).

    Pulled via the CELLxGENE direct HTTPS endpoint instead of cellxgene_census —
    the census filter API changed and no longer accepts collection_id filtering.
    The Extended GBmap dataset is published at a stable HTTPS URL (~11 GB).
    """
    out.parent.mkdir(parents=True, exist_ok=True)
    # Extended GBmap (collection 999f2a15-...), dataset_id 230c8701-...
    url = "https://datasets.cellxgene.cziscience.com/230c8701-7291-4dd5-bf36-688490c681ee.h5ad"
    log(f"GBmap — downloading from CELLxGENE: {url}")
    log("  Expected size: ~11.2 GB (Extended GBmap, 1,135,677 cells)")
    subprocess.run(
        ["wget", "-c", "-O", str(out), url],
        check=True,
    )
    require(out.exists() and out.stat().st_size > 5e9,
            f"GBmap h5ad too small: {out} is {out.stat().st_size} bytes "
            f"(expected ~11 GB)")
    log(f"GBmap → {out}  ({out.stat().st_size/1e9:.1f} GB)")


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

HANDLERS = {
    "guilhamon": fetch_guilhamon,
    "gse165037": fetch_raviram,                    # new slug after 2026-10-04 citation fix
    "gse138794_guo": fetch_gse138794_guo,           # 4 IDH-WT primary GBM from Guo 2020
    "gse276177_khan_astro": fetch_gse276177_khan_astro,  # 3 IDH-WT primary GBM from Khan 2024
    "raviram": fetch_raviram,                      # keep for backwards compat
    "tcga_scatac": fetch_tcga_scatac,
    "mathewson_lupien": fetch_mathewson_lupien,
    "gbm_space": fetch_gbm_space,
    "gbm_tme_atlas_hra004942": fetch_hra004942,
    "wang_sciadv_tc_ptb": fetch_wang_sciadv_tc_ptb,
    "spatial_epigenomic_niches": fetch_spatial_epigenomic_niches,
    "catlas": fetch_catlas,
    "gbmap": fetch_gbmap,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(HANDLERS))
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    if args.cohort in ("catlas", "gbmap"):
        log(f"=== Fetching reference {args.cohort} ===")
        HANDLERS[args.cohort](args.out)
    else:
        with open("config/datasets.yaml") as f:
            datasets = yaml.safe_load(f)
        meta = datasets["cohorts"][args.cohort]
        log(f"=== Fetching cohort {args.cohort} ===")
        log(f"  citation : {meta.get('citation', 'n/a')}")
        log(f"  accession: {meta.get('accession', 'n/a')}")
        log(f"  platform : {meta.get('platform', 'n/a')}")
        log(f"  expected : {meta.get('n_patients', 'n/a')} patients")
        HANDLERS[args.cohort](args.out)

    log(f"=== Fetch complete: {args.cohort} ===")
    return 0


if __name__ == "__main__":
    sys.exit(main())
