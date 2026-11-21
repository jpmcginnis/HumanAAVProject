# Verified deposit layouts (checked 2026-10-04)

Reality-checked against NCBI eutils and EBI BioStudies API before writing
`src/fetch.py`. If a deposit changes shape, re-verify this file before
re-running fetchers.

---

## GSE139136 — Guilhamon

At `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139136/suppl/`:

- `GSE139136_RAW.tar` (161 MB) — this contains the per-GSM supplementary files
- `filelist.txt`

**Implication:** the scaffold's `wget -r` finds exactly these two files
and never produces `.fragments.tsv.gz`. Must untar first, then glob
inside the extracted contents.

---

## GSE165037 — Raviram (actually Guilhamon lab IDHwt+IDHmt mix)

eutils returns 12 samples, of which only 3 are IDH-WT snATAC:

| GSM        | Sample                     | Keep?               |
|------------|----------------------------|---------------------|
| GSM5024962 | GBM4_IDHWT_snATAC          | yes                 |
| GSM5024963 | GBM9_IDHWT_snATAC          | yes                 |
| GSM5024964 | GBM12_IDHWT_snATAC         | yes                 |
| GSM5024965 | GBM1_IDHMT_snATAC          | no — IDH-mutant     |
| GSM5024966 | GBM11_IDHMT_snATAC         | no — IDH-mutant     |
| GSM5024967 | Non_tumor_brain_snATAC     | no — non-tumor      |
| GSM5024968 | GBM4_IDHWT_snRNA           | pair (for Multiome) |
| GSM5024969 | GBM9_IDHWT_snRNA           | pair                |
| GSM5024970 | GBM12_IDHWT_snRNA          | pair                |
| GSM5024971 | GBM1_IDHMT_snRNA           | no                  |
| GSM5024972 | GBM11_IDHMT_snRNA          | no                  |
| GSM5024973 | Non_tumor_brain_snRNA      | no                  |

**Discrepancies with `config/datasets.yaml`:**
- The config calls this "Raviram et al. 2023 PNAS" but the eutils title
  matches a Guilhamon-style IDHwt-vs-IDHmt study. Need to verify the
  citation is correct; the config may be mislabeled.
- The config claims "RL3 multi-region case (4 sections)" — no such
  sample in this accession.

**Implication:** fetcher must filter by GSM title (keep only `GBM*_IDHWT_snATAC`)
before pulling.

---

## E-MTAB-17183 — GBM-Space (De Jong 2025)

Verified via `https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-17183/files`.
165 files, flat (no per-sample subdirectories):

| What                                              | Count | Size   |
|---------------------------------------------------|-------|--------|
| `cellranger-arc201_count_<hash>.gz`               | 155   | ~10 GB each (~1.5 TB total) |
| `GBM_space_ATAC_filtered_peaks.h5ad.gz`           | 1     | 22.4 GB |
| `GBM_space_snRNA.h5ad.gz`                         | 1     | 11.3 GB |
| `rna_metacell_scdori_11_03_24.h5ad`               | 1     | 0.2 GB  |
| `genescore_metacell_scdori_11_03_24.h5ad`         | 1     | 1.7 GB  |
| `GBM_space_ATAC_filtered_peaks_README.md`         | 1     | small   |
| `GBM_space_snRNA_README.md`                       | 1     | small   |
| `README_cellranger_arc.md`                        | 1     | small   |
| `README_metacell.md`                              | 1     | small   |
| `E-MTAB-17183.sdrf.txt`                           | 1     | small — maps hash → donor/region |
| `E-MTAB-17183.idf.txt`                            | 1     | small   |

**Implication:**
- The scaffold's `*annotated*multiome*.h5ad` glob matches nothing. The actual
  annotated ATAC object is `GBM_space_ATAC_filtered_peaks.h5ad.gz`.
- The scaffold's `rglob("outs")` fallback finds nothing — files are flat with
  opaque UUID hashes, not Cell Ranger-style directory trees.
- We only need **22 GB** (the ATAC h5ad) to build the enhancer atlas, not the
  full 1.5 TB of raw tarballs. The scRNA h5ad and metacells add ~13 GB more
  if we want multiome enhancer-to-gene linkage.
- Hash-filename-to-donor mapping requires parsing `E-MTAB-17183.sdrf.txt`.

---

## Still to verify

- **GSE174554 (Mathewson/Lupien)** — need to call eutils and confirm which GSMs
  are IDH-WT snATAC primary vs recurrent, how many fragments.tsv.gz files are
  deposited vs how many need SRA re-processing from FASTQ.
- **TCGA-ATAC-Seq-2024** — the S3 archive at
  `s3://jpm-atacseq-archive-2026/gbm_pipeline_complete/tcga_scatac/` is assumed
  to contain `integrated.h5ad`. Confirm on first spin-up.
- **CATLAS (GSE244618)** — assumed to be at
  `s3://jpm-atacseq-archive-2026/catlas/catlas_full_annotated.h5ad`. Confirm.
- **GBmap** — must be pulled via cellxgene_census, not wget. Need to verify
  which census version has the Ruiz-Moreno 2025 GBmap collection and what the
  donor-level filter looks like.

---

## Verified ingest strategies (what the rewritten fetchers/ingestors do)

- **Guilhamon:** untar `GSE139136_RAW.tar` → look for fragments inside;
  if the tar only has matrices (no fragments), fail loudly rather than
  silently producing an empty AnnData.
- **Raviram:** enumerate GSMs via eutils, keep only those whose title matches
  `GBM.*IDHWT_snATAC`, download per-GSM supplementary.
- **Mathewson/Lupien:** same enumerate-via-eutils pattern.
- **GBM-Space:** pull only `GBM_space_ATAC_filtered_peaks.h5ad.gz` +
  `GBM_space_snRNA.h5ad.gz` + the three READMEs + the SDRF — not the 155
  tarballs. Read the ATAC h5ad directly; use SDRF to validate donor labels.
- **TCGA:** S3 sync (unchanged); verify `integrated.h5ad` exists post-sync.
- **CATLAS/GBmap:** two new producer rules in the Snakefile.
