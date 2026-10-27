# Per-cohort ingest reference

33 cohorts contributed to the NB atlas (713 atlas-eligible donors after case exclusion). Each cohort shipped at a different layer of the processing stack (raw FASTQ, 10x-matrix, ArchR RDS, Seurat RDS, SnapATAC2 h5ad, CELLxGENE h5ad, custom .mtx). This document lists every cohort with the ingest route and script. See [`../src/ingest/`](../src/ingest/) for the ingest scripts.

Every cohort pseudobulk lands as `<cohort>_direct_nnz1k.parquet` with columns:
`peak_num, major_class, donor_id, n_cells_access_direct, n_cells_total_direct, fraction_open_direct, study_id`.

QC gate universal to all cohorts: **cells with nnz ≥ 1000 against the 990K union peak set are retained**, before pseudobulking.

## Table

| cohort | n_donors | tissue | deposit layer | script |
|---|---|---|---|---|
| `hmba-bg` | 204 | basal ganglia (NAcc, caudate, putamen, GP, STN, SNr, SNc) | SnapATAC2 h5ad per donor on Allen S3 | [`hmba_bg_process.py`](../src/ingest/hmba_bg_process.py) + [`hmba_bg_continuation.py`](../src/ingest/hmba_bg_continuation.py) |
| `hmba-bg-retry` | 8 | same (donors that FAILed first download pass) | SnapATAC2 h5ad local dir (post-retry-dl) | [`hmba_bg_retry_process.py`](../src/ingest/hmba_bg_retry_process.py) |
| `menon_corces_pd` | 101 | multi-region PD atlas (cerebellum, cingulate, MTG, putamen, nigra) | CELLxGENE h5ad (open-access) | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `seaad` | 86 | MTG + DLPFC + Caudate (SEA-AD) | per-donor h5ad on `s3://sea-ad-single-cell-profiling/` | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `4dn_hippo` | 37 | hippocampus | permissive barcode-called h5ad on 4DN portal (median nnz 3–4; ~98% dropped by nnz≥1k) | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `GSE277313` | 78 | — | GEO h5ad | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `GSE336954` | 63 | central amygdala, AUD atlas | GEO 10x matrix | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) + per-cohort shim |
| `GSE322785` | 55+ | — | GEO permissive barcode h5ad (same median-nnz caveat as 4DN) | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `GSE174367` | ~24 | Morabito/Swarup AD cortex | GEO 10x matrix | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `cxg1` (Braun/Linnarsson) | 8 | first-trimester fetal brain (whole-brain regions) | CELLxGENE h5ad | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `cxg-zhu-roussos` | 12 | developing cerebral cortex | CELLxGENE h5ad via GSE204682 | R-side RDS extraction + `pseudobulk_from_h5ad.py` |
| `cxg-takeuchi-als` | 17 (8 HC + 9 ALS; brain + spinal cord as separate donor_ids) | motor cortex + lumbar spinal cord | Seurat `*_so_integrated_LinkePeaks.rds` on figshare | [`pseudobulk_takeuchi_als_rds.R`](../src/ingest/pseudobulk_takeuchi_als_rds.R) |
| `cxg-lattke-ds` | 30 (15 CON + 15 DS fetal) | fetal frontal cortex | CELLxGENE h5ad | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `geo-GSE246443` | 39 PASS controls (MDD cases excluded) | postmortem dlPFC (MDD study) | 31.5 GB uncompressed `.mtx` with 2.07 B nnz entries, sorted by cell_idx | [`gse246443_stream_mtx.py`](../src/ingest/gse246443_stream_mtx.py) |
| `GSE168408` | — | Herring PFC multiome (fetal→adult) | GEO | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `GSE193688` / `GSE207334` / `GSE212630` / `GSE214979` / `GSE259298` / `GSE266305` (iHS controls) / `GSE268807_margin` / `GSE272082` / `GSE275224` / `GSE288939` / `GSE296943_C1` / `GSE296943_C2` / `GSE305146` / `GSE325789` / `GSE330130` / `GSE331288` / `GSE338484` | 1–40 each | mixed | GEO h5ad or 10x matrix | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |
| `veiga` | 5 | adolescent cortex | CELLxGENE h5ad | [`pseudobulk_from_h5ad.py`](../src/ingest/pseudobulk_from_h5ad.py) |

Per-donor data URLs + PMID/DOI for every row above live in the metadata CSV at [data.mcginnislab.org/normal-brain-atlas/metadata.csv](https://data.mcginnislab.org/normal-brain-atlas/metadata.csv).

## Standard ingest flow (per cohort)

1. **Download** the deposit (raw FASTQ, 10x matrix, h5ad, mtx, RDS, SnapATAC2 h5ad).
2. **Load** into a per-cell matrix (fragment counts per cell per peak), using the cohort's native peak set.
3. **Overlap** the native peaks → the 990K-peak atlas union set via [`pyranges`](https://biocore-ntnu.github.io/pyranges/) (longest-overlap-wins). Each native peak maps to zero, one, or multiple union peaks.
4. **Binarize** cell × union_peak to accessibility.
5. **nnz filter**: keep cells with ≥ 1000 union peaks accessible.
6. **Label**: donor_id (from obs/metadata), major_class (cell-type-ontology-mapped from obs; see `src/atlas/nb_label_transfer.py` where author labels are missing).
7. **Pseudobulk per (donor, major_class)**: for each (peak, donor, class) group, record `n_cells_access`, `n_cells_total`, `fraction_open` where `n_cells_access ≥ 3` (otherwise dropped).
8. **Write** `<cohort>_direct_nnz1k.parquet`. Delete raw after confirming write.

## Non-standard ingests (why specific scripts exist)

- **HMBA-BG** (`hmba_bg_process.py` + `hmba_bg_continuation.py`): Allen's SnapATAC2 h5ads store fragments in `obsm/fragment_paired` as `(position, length)` pairs in a cells × total_genome_bp sparse CSR. The ingest decodes to `(chrom, start, end)`, overlaps with union peaks, binarizes, nnz-filters, then aggregates. 188 donors × ~3-5 GB each, so `hmba_bg_continuation.py` is the stateful batch runner that checkpoints after every donor and resumes on restart. `hmba_bg_retry_process.py` is the parallel re-attempt for the 8 donors whose first `aws s3 cp` failed.
- **GSE246443** (`gse246443_stream_mtx.py`): 31.5 GB uncompressed `.mtx` (200,071 cells × 825,492 peaks, 2.07 B nnz). Can't fit on disk, so we stream via `tar -xzO ... | python` cell-by-cell (the mtx is sorted by column 2 = cell_idx). ~50 min.
- **Takeuchi ALS** (`pseudobulk_takeuchi_als_rds.R`): figshare deposit is Seurat `*_so_integrated_LinkePeaks.rds` (5.6 GB brain + 3.2 GB spinal). Loaded in R, overlapped with union peaks, pseudobulked, written as TSV then converted to parquet.

## Cohort cases excluded

- **GSE266305 FCD IIIa × 4** → `EXCLUDE_epilepsy` (Series title explicitly "epileptogenic cortex").
- **GSE288939 epilepsy MEA × 1** → `EXCLUDE_case`.
- All ALS / PD / Down syndrome / MDD disease arms → `EXCLUDE_case`.
- `GSE325789` → `EXCLUDE_mixed` (patient-type ambiguous).

All exclusions are in `atlas_usability_tier` in the metadata CSV.

## Blocked / deferred (not ingested)

See `audit_2026-10-25/public_brain_snATAC_catalog.md` in the atlas-build directory for the full list of deposits examined. Blocked:
- SEA-AD raw multiome FASTQ, Clarence postnatal, Developing Human Brain snATAC — dbGaP / DUC.
- GSE214637 (bigwig only), GSE216270 (gene-activity only), GSE283715 (pooled no demux), GSE338996 (raw-only).
- CATlas raw FASTQ (GSE244618, 4.56 TB on NeMO) and Pančíková/Aerts PD (Zenodo 20077637, ASAP CRN FASTQ) — open-access but EC2-only (would need alignment pipeline).
