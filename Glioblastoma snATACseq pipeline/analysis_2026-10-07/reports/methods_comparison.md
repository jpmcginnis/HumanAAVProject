# snATAC-seq / Multiome methods comparison — 8-cohort GBM enhancer atlas

**Compiled 2026-10-07.** Pulled from paper XML/HTML (Europe PMC, bioRxiv JATS, PMC) and GEO SOFT metadata; three cohorts' primary papers (Wang 2019 *Cancer Discovery*, Terekhanova 2024 *Science*, Lin/Qu Kun 2026 *Nat Neurosci*) are behind paywalls that WebFetch could not clear — their rows below are filled from GEO SOFT, the TCGA-ATAC-Seq-2024 data portal, and the NGDC project record. Manual pull URLs are noted in those rows and in the Appendix.

## 1. Summary table

| Cohort (slug) | Accession / primary paper | Cases, inclusion, IRB | Tissue handling | Nuclei isolation | Library prep (platform / chip / loading) | Sequencer & reads | Processing / peak calling |
|---|---|---|---|---|---|---|---|
| **guilhamon** | GSE139136 — Guilhamon P … Lupien M, *eLife* 2021, PMID 33427645 | 4 primary GBM, all IDH-wt (G4218 M64, G4250 M73, G4275 F52, G4349 M62). IRB: SickKids Toronto + Univ Calgary + HREBA Alberta. Fresh from OR. | Fresh fragments 0.3-0.7 cm³ blunt-dissected and slow-frozen in NeuroCult/BSA/DMSO freezing media on CoolCell at -80 °C. Cryopreserved → thawed, PBS wash to remove DMSO. | **Detergent/mechanical, no enzyme.** ATAC resuspension buffer (Tris-HCl pH7.4, NaCl, MgCl₂, 0.1% NP-40, 0.1% Tween-20, 0.01% Digitonin, 1% BSA/PBS); wide-bore P1000 + vortex + 10 min on ice; wash → ATAC-Tween wash buffer; 40 µm cell-strainer FACS tube. | **10x Chromium Single Cell ATAC v1.0** (first-gen scATAC kit). Loading target not stated. | **Illumina NextSeq 500**, 50 bp paired-end. | `cellranger-atac mkfastq` v1.0.0; `cellranger-atac count` v1.1.0; **hg19** reference v1.1.0. Downstream: chromVAR + **Signac v1.4.1**. CNV: CONICSmat on Signac gene-activity matrix. Peak calling: via Signac default. |
| **mathewson_lupien** | GSE174554 — Wang L … Diaz AA, *Nature Cancer* 2022, PMID 36539501 (Mathewson/Lupien+Diaz) | Longitudinal GBM specimens (paired primary/recurrent) — see paper Table 1 / "Wang_pair_*" IDs. IRB: UCSF, de-identified, Declaration of Helsinki. All de-identified. | Fresh-frozen and FFPE from UCSF Neurosurgery Tissue Bank. snATAC used **frozen only**. | **Mechanical Dounce** in Sigma lysis buffer, 40 µm strainer, pellet, wash with 10x nuclei wash buffer, **sucrose density gradient**, resuspend in 10x tagmentation buffer. | **10x Genomics scATAC (standard)** — "10x Genomics platform as per manufacturer's protocol" (methods do not state v1 vs v1.1 explicitly; CellRanger ATAC v1.1.0 downstream suggests v1 / early v1.1 chemistry). snRNA loaded ~15,000 nuclei/capture. | **Illumina NovaSeq**, 10x-recommended parameters (paired-end). | **CellRanger ATAC v1.1.0** → snapATAC pkg (r3fang/SnapATAC). QC: fragments > 1000 and promoter ratio > 0.2. Clustering: Seurat v3 SNN on snapATAC gene-body accessibility. chromVAR v1.6.0 TF deviations; `findDAR`, `runMACSForAll`, `runHomer`; deepTools v3.4.0. CNV: CONICSmat v1.0 on snapATAC gene activity. |
| **gbm_tme_atlas_hra004942** | HRA004942 (NGDC/GSA-Human) — Lin, Chen, Li, Chen, Fang … Qu Kun, *Nature Neuroscience* 2026, DOI 10.1038/s41593-026-02265-5 | **Not accessible via WebFetch** (Nature behind IDP auth, NGDC project page returns bare metadata). From NGDC: USTC, Qu Kun lab; "Spatial and single-cell characterization of human glioblastoma tumor microenvironment reveals malignant cellular communities." Open-access release 2026-02-17. 11 patient IDs in local audit (P62, P64, P77, P78, P79, P80, P83, P84, P92, P98, P101). | **Pull manually.** The paper and NGDC project page both need credentialed browser access. Try: https://www.nature.com/articles/s41593-026-02265-5 (through institutional access) and https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004942 (direct in browser, not WebFetch). | Pull manually. | Pull manually — the local audit shows ~500-800 cells per patient post-QC, which is **two orders of magnitude lower** than every 10x-based cohort here. That is consistent with **DNBelab C4** or another low-recovery platform, but confirmation requires the methods PDF. | Pull manually. | Pull manually. |
| **gse138794_guo** | GSE138794 — Wang L … Diaz AA, *Cancer Discovery* 2019, PMID 31554641 (same Diaz lab as mathewson_lupien) | "scRNA-seq / snRNA-seq of 28 gliomas and **scATAC-seq for 8 cases**" (series summary). IRB: UCSF Neurosurgical Tissue Bank. Raw data not public (dbGaP) — only processed matrices on GEO. Local audit has 3 cases: SF11215, SF11956, SF11979. | Frozen tissues for snATAC. | **"Frakenstein" protocol** (Luciano Martelotto, Melbourne VCCC — 10x customer-developed protocol). Mechanical Dounce in detergent buffer, debris removal, nuclei wash. Same lab/protocol as mathewson_lupien. | **10x Chromium Single Cell** capture chip; **~15,000 nuclei loaded per capture** for frozen; fresh tissues: 10.2 µL live cells @ 1700/µL. Methods field does not say ATAC v1 vs v1.1 — but **CellRanger ATAC v1.1.0** downstream suggests v1 chemistry. | **Illumina NovaSeq 6000**, paired-end 100 bp. | **CellRanger ATAC v1.1.0**, hg38. Processed data on GEO are sparse matrix (mtx) + barcodes.tsv + peaks.bed per sample; peak calling via CellRanger default. |
| **tcga_scatac** (ex-sundaram_gbm) | TCGA-ATAC-Seq-2024 (GDC) — Terekhanova NM … Chang HY/Corces MR, *Science* 2024, PMID 39236169 | **Not accessible via WebFetch** (Science and GDC publication page both paywalled/JS-walled). Study scope: pan-cancer scATAC across 8 TCGA tumor types, 74 samples total; GBM is one of them. 10 GBM samples in local audit, UUID-style IDs (TCGA barcodes). | Pull manually. The TCGA-ATAC-Seq-2024 GDC portal (https://gdc.cancer.gov/about-data/publications/TCGA-ATAC-Seq-2024) hosts a supplementary methods PDF behind a JS gate. | Pull manually. | Pull manually. Companion ATAC-AWG 2018 bulk ATAC used ascent-chip 10x v1 Omni-ATAC; the 2024 scATAC supplement should state 10x scATAC v1.1 vs Multiome. | Pull manually. | Pull manually. Downstream in the Science paper is ArchR per supplementary methods. |
| **gbm_space** | E-MTAB-17183 — de Jong G, … Saraswat M, … Bayraktar OA (corresp., Wellcome Sanger), *bioRxiv* 2025, DOI 10.1101/2025.05.13.653495 ("A spatiotemporal cancer cell trajectory underlies glioblastoma heterogeneity") | **12 GBM tumours** (AT3-AT15 in local audit matches this), multiple regions per tumour, matched whole blood for germline. IRB: Cambridge Local REC 18/EE/0172, Declaration of Helsinki 2000. Consent pre-op; 5-ALA guided resection. | **Fresh-frozen** in OCT, isopentane bath -75 °C via dry-ice. Pre-screen: 2-3× 10 µm sections → **Qiagen RNA extraction**, Tapestation RIN; **only RIN > 7 blocks used**. Nuclei sectioning: 50 µm thick, 500-900 µm total volume per block, collected in pre-chilled homogenization tubes on dry ice. | **Mechanical Dounce** glass homogenizer in isolation buffer (3 mM MgCl₂, 10 mM NaCl, 10 mM Tris pH 7.4, 1 mM DTT, 0.1% Tween-20, 0.1% NP-40, 1% BSA, **0.01% Digitonin**) + Protector RNase Inhibitor 0.2 U/µL. 10 strokes pestle A → 10 strokes pestle B → 40 µm filter → 500 RCF → PBS/BSA storage buffer with RNase inhibitor 1 U/µL. Trypan blue count; **Percoll gradient debris removal**. | **10x Chromium Next GEM Single Cell Multiome ATAC + Gene Expression**. Target **5,000-10,000 nuclei per reaction**; 2-3 reactions per tumour site. 3′ Reagent Kits v3 User Guide for cDNA/library construction. | **Illumina NovaSeq 6000, S4 Flowcell**; min **100,000 read-pairs/nucleus/modality**. GEX: R1 28 / i7 10 / i5 10 / R2 90. ATAC: R1N 50 / i7 8 / i5 24 / R2N 49. | **Cell Ranger ARC v2.0.1** with custom GRCh38 3.0.0 pre-mRNA + ARC 2.0.1 ATAC genome. CellBender v0.2.0 for empty-droplet correction. RNA QC in Scanpy v1.9.3 (UMI>1000, genes>500, %MT<10). Scrublet v0.2.3 + MAD-based doublet filtering. **ArchR workflow for ATAC** (min fragments>1000, >4 TSS). epiAneufinder (1 MB window) for CN on scATAC. |
| **gse276177_khan_astro** | GSE276177 — Sojka C, Wang HV, Bhatia T … Sloan SA, *Nature Cell Biology* 2025, PMID 39779941, DOI 10.1038/s41556-024-01583-9 ("Mapping the developmental trajectory of human astrocytes reveals divergence in glioblastoma"). (The GEO working title "Divergence from the human astrocyte developmental trajectory in glioblastoma" is the paper's conceptual title — not "Khan et al."; the slug is a misnomer carried from an earlier pipeline revision.) | **3 patients × matched tumour + margin** = 10 Multiome libraries total (GBM20 M IDHwt, GBM25 F IDHwt, GBM38 M IDHwt). IRB: Emory School of Medicine. Local audit has 3 patients, no CNV computed (cohort flagged `CNV_not_run_for_cohort`). | **Frozen tissue**. Tumor+margin samples placed in 4 °C Hibernate-A immediately after resection; **tissue dissociation within 1 h post-resection** (per Sojka 2025 methods extract). ~20 mg tissue per prep. | **Mechanical Dounce** 2 mL homogenizer in homogenization buffer (0.26 M sucrose, 0.03 M KCl, 0.01 M MgCl₂, 0.02 M Tricine-KOH pH 7.8, 0.001 M DTT, 0.5 mM Spermidine, 0.15 mM Spermine, 0.3% NP40, cOmplete protease inhibitor). **40 µm Flowmi then 20 µm bucket-style strainer**. **Iodixanol step gradient (25% / 30% / 40%)**; nuclei band collected from 30/40% interface. 1-2 washes in ATAC-RSB-Tween (Tris-HCl pH 7.5, NaCl, MgCl₂, Tween-20). | **10x Chromium Single Cell Multiome ATAC + Gene Expression** kit (10x Cat 1000285). **16,100 nuclei loaded, 10,000 target capture per sample**. | **Illumina NovaSeq 6000**, 2 × 150 bp paired-end; target ≥ **50,000 read-pairs/nucleus**. | **cellranger-arc v2.0.0**, hg38. (Note: slightly older than gbm_space's 2.0.1.) |
| **gse165037** | GSE165037 — Raviram R, Chen C, Ren B (Ren lab UCSD), 2021 GEO (no linked PMID in the series; the matching paper is the Raviram/Ren pan-NFI GBM study — the GEO series title "Single Cell Analysis of Chromatin Accessibility Reveals Genetic and Regulatory Heterogeneity in Glioblastomas" is a working title; the full-text paper has not been consistently indexed to this series). | **5 GBM + 1 non-tumor brain** = GBM1 (IDH1-mut), GBM4 (IDHwt), GBM9 (IDHwt), GBM11 (IDH1-mut), GBM12 (IDHwt). Raw data on dbGaP (patient privacy). Local audit has GBM4, GBM9, GBM12 (all IDHwt). | **Ground frozen tumor tissue.** | **Detergent-only, no enzyme.** NPB = 5% BSA, 0.2% IGEPAL-CA630, cOmplete protease inhibitor, 1 mM DTT in PBS. (For paired snRNA: 2% BSA, 0.2% Triton-X, cOmplete PI, 1 mM DTT, 0.2 U/µL RNasin, then SH800 Sony sort.) | **Combinatorial-indexing snATAC-seq, plate-based, NOT droplet-based.** Nuclei tagmented with Tn5 in 96 wells (**first barcode**), pooled, then 20-25 nuclei sorted into 768 wells where **PCR introduces a second barcode**. (snRNA uses a separate pipeline: 10x Chromium Single Cell 3′ v2 kit.) | **Illumina HiSeq 2500, HiSeq 4000, or NextSeq 500** — paper used a mix. (NovaSeq not used.) Read-length not stated in SOFT. | ATAC reads → **bowtie2** → hg38; **Snaptools → SnapATAC**, 5 kb bins; filter cells with >1000 reads AND **promoter ratio > 20%**. snRNA → STAR → hg38, Seurat v3 (>200 reads, <5% MT). |

## 2. Protocol clustering

### A. "Modern 10x Multiome, NovaSeq 6000" (highest sensitivity, newest kits)
- **gbm_space** (10x Multiome ARC Next GEM; cellranger-arc 2.0.1; 100 K rp/nucleus; ArchR downstream)
- **gse276177_khan_astro / Sojka 2025** (10x Multiome; cellranger-arc 2.0.0; 50 K rp/nucleus)

These two should be the "cleanest" library-wise. Both use a Dounce + detergent + gradient (iodixanol or Percoll) + 10x Multiome. Expected high per-cell fragment counts, strong TSS enrichment, good FRiP. **gbm_space has roughly 2× the sequencing depth per nucleus of gse276177 and used Percoll vs iodixanol** — a small but noteworthy difference.

### B. "First-gen 10x scATAC + NovaSeq 6000 or NextSeq 500, Diaz-lab Frankenstein pipeline"
- **mathewson_lupien** (10x scATAC, CellRanger ATAC v1.1.0, NovaSeq)
- **gse138794_guo** (10x Chromium single cell capture, 15 K nuclei loaded, NovaSeq 6000 PE100, CellRanger ATAC v1.1.0)
- **guilhamon** (10x Chromium Single Cell ATAC v1.0, NextSeq 500, 50 bp PE, Signac v1.4.1)

Same **vintage of 10x chemistry** (circa 2019–2020). Guilhamon is the oldest reference genome (hg19); the Diaz-lab pair uses hg38. These cohorts should have lower per-cell fragments than the Multiome cohorts but still clean peak structure. Guilhamon's 50 bp PE on NextSeq 500 is the shallowest per-read run of all; expect lower read-pair diversity than the NovaSeq-based cohorts.

### C. "Combinatorial-indexing snATAC, HiSeq/NextSeq, plate-based" — **idiosyncratic**
- **gse165037** (Raviram/Ren lab sci-ATAC-style: Tn5 in 96w → sort → 768w PCR; bowtie2 + SnapATAC on 5 kb bins; dbGaP-controlled raw)

This is **not a droplet method** and does not generate 10x-style barcode-structured fragments. Its per-cell fragment distribution and FRiP will look substantively different from any 10x cohort. The processed data on GEO are RDS objects of SnapATAC-called cells only, no fragment files. This cohort is the biggest methods outlier in the atlas; any pan-cohort peak set or `fragments.tsv.gz` reprocessing has to accommodate the fact that this one does not have fragment-level data from CellRanger.

### D. "Unknown / to be filled" — flag before pan-malignant interpretation
- **gbm_tme_atlas_hra004942** (Qu Kun 2026 Nat Neurosci) — **suspected DNBelab C4** given per-sample cell yields of 500-800 in the local audit (vs 10x's typical 5–40 K). The paper's methods must be read to confirm.
- **tcga_scatac** (Terekhanova 2024 *Science*) — likely 10x scATAC v1.1 and ArchR per the Chang/Corces group's standard, but confirm from Science supplementary methods.

## 3. Shared pitfalls to watch

| Pitfall | Which cohorts | Notes |
|---|---|---|
| **Frozen-tissue only** (no fresh in snATAC path) | gse138794_guo, mathewson_lupien, gse165037, gbm_space, gse276177_khan_astro, guilhamon (fresh → cryopreserved first) | All of the 10x-based cohorts here used frozen input. This is standard, so no concern — but post-mortem interval / freeze delay is only documented explicitly for gse276177 (within 1 h) and gbm_space (immediate OCT in isopentane). The others do not state time-to-freeze. |
| **No FFPE in the snATAC lanes** | all cohorts | mathewson_lupien's paper discusses FFPE for spatial, but snATAC was fresh-frozen only. Nothing in the atlas is from FFPE, which is good. |
| **Different nuclei isolation chemistry** can bias TSS enrichment and % MT | all cohorts differ: Guilhamon uses NP-40 + Digitonin + Tween-20 (ATAC-RSB-Omni); gbm_space uses NP-40 + Tween-20 + Digitonin + RNase inhibitor + Percoll; gse276177 uses sucrose/Tricine + NP-40 + iodixanol; mathewson_lupien uses Sigma lysis + sucrose gradient; gse165037 uses IGEPAL-CA630 (no Tween, no Digitonin) + no gradient | Expect systematic shifts in fragment-length distribution. **Guilhamon and gbm_space (both use Digitonin, which lyses the plasma but not nuclear membrane cleanly) should be most comparable.** gse165037's IGEPAL-only prep tends to give higher nucleosomal background. |
| **Mixed sequencers** | gse165037 was sequenced on HiSeq 2500 / HiSeq 4000 / NextSeq 500 across the series (per-sample instrument not disclosed in SOFT) | Possible batch effects within gse165037 alone, before any cross-cohort comparison. |
| **Short read length** | guilhamon (NextSeq 500 50 bp PE) is the shortest; gse138794_guo and gbm_space are 100 bp / 150 bp | Shorter reads = reduced mappability in low-complexity / repeat regions, more read loss. Guilhamon's peaks will have narrower effective coverage. |
| **Non-droplet method mixed with droplet methods** | gse165037 (sci-ATAC-style combinatorial indexing) | Different barcode structure, no 10x-style fragments.tsv from CellRanger. If this cohort was fed into a 10x-centric pipeline (CellRanger ATAC / snapatac2 fragment reader), it needs a bespoke import. Check that the processing step correctly handled the RDS → fragment conversion. |
| **CellRanger ATAC version split: 1.0.0 → 1.1.0 → arc 2.0.0 → arc 2.0.1** | guilhamon used 1.0.0 (mkfastq) + 1.1.0 (count); mathewson_lupien and gse138794_guo used 1.1.0; gse276177 used arc 2.0.0; gbm_space used arc 2.0.1 | Peak-calling internals (and the default fragments.tsv column format) changed between these. If the pan-malignant matrix was built from heterogeneous fragments files, confirm the schema matches. |
| **Reference genome split** | guilhamon = **hg19**; all others = hg38 | Guilhamon's native peaks need liftover to hg38 before any pan-cohort peak set is built. If this is already done, confirm liftover success rate (typical loss ~2–5% of peaks). |
| **Suspected non-10x platform** | gbm_tme_atlas_hra004942 (very low cells/sample suggests DNBelab C4 or similar) | DNBelab C4 chromatin recovery is known to be lower-complexity than 10x; MGISEQ reads have different quality profiles than Illumina. If confirmed, this cohort may contribute lower-quality peaks to any pan-cohort union. **Flag before including in the pan-malignant top-ranked enhancer set.** |

## 4. Preliminary "good" vs "bad" preps — what to check once the v2 audit CSV is on this Mac

The audit CSV at `/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results/qc/outputs/01_per_patient_audit.csv` has `total_cells_post_qc`, `n_cells_pan_malignant`, and `frac_pan_malignant`, but **not** TSS enrichment, FRiP, nFragments, or per-sample doublet rate. The v2 CSV (still on EC2) is expected to add those.

What the current audit already tells us (cells per patient post-QC):

| Cohort | Range of total_cells_post_qc | Comment |
|---|---|---|
| gbm_space | 34,702 – 192,661 (avg ~87 K) | **Highest cell recovery** — consistent with Multiome ARC on NovaSeq S4 at 100 K rp/nucleus, Percoll + Dounce. |
| gse276177_khan_astro | 39,744 – 129,846 (avg ~96 K) | Also very high — Multiome, iodixanol gradient. |
| mathewson_lupien | 2 – 8,058 (avg ~2.4 K excl. the n=2 outlier) | Expected lower (first-gen scATAC) and **highly uneven** — pair_16 had only 2 cells and is already excluded. |
| gse138794_guo | 766 – 3,628 | Low-moderate; consistent with 10x v1 chemistry. |
| guilhamon | 629 – 1,442 | Low cells, as expected for NextSeq 500 + first-gen scATAC v1.0. Still gives usable malignant fractions (0.48–0.75). |
| gse165037 | 8 – 3,249 | Highly variable; the combinatorial-indexing method's cell recovery depends on the sort gate and plate count. GBM9 (8 cells) is already excluded. |
| gbm_tme_atlas_hra004942 | 158 – 761 | **Striking** — two orders of magnitude lower than the Multiome cohorts. Not sequencing depth alone — likely platform (DNBelab C4?) or very stringent per-cell QC. |
| tcga_scatac | 5,431 – 64,647 | Mid-range; multiple TCGA samples, variable depth. CNV not computed (CNV_not_run_for_cohort flag). |

**Hypothesis to validate once v2 lands** (needs TSSe and FRiP):
- **Likely top-quartile TSSe** (clean 10x Multiome on high-RIN frozen with gradient purification): gbm_space, gse276177_khan_astro
- **Likely middle** (10x v1/v1.1 scATAC on frozen, standard prep): mathewson_lupien, gse138794_guo, guilhamon, tcga_scatac
- **Likely lower / outlier** (combinatorial indexing, non-droplet): gse165037 — expect **lower FRiP but reasonable TSSe**; fragment-length distribution may show more sub-nucleosomal reads
- **Unknown but suspect** (very low cell yield, suspect platform): gbm_tme_atlas_hra004942 — treat as "needs manual QC inspection"

**Checks to run once v2 CSV arrives:**
1. Plot `tss_enrichment` vs cohort, as box or violin. If gbm_tme_atlas_hra004942 is >1 SD below the 10x cohorts, flag.
2. Plot `frip` vs cohort. Expect gse165037 ~0.2–0.4 (SnapATAC-style, no CellRanger peak set); 10x cohorts ~0.4–0.6.
3. Pull `nFrags_per_cell` median per cohort. Guilhamon and gbm_tme_atlas_hra004942 are the two likeliest to be low.
4. Compare `doublet_rate` across cohorts — Multiome cohorts loaded at 10 K target should have ~5–8%; first-gen scATAC at 15 K-loaded will be higher (~10%).
5. For gse165037 (combinatorial indexing), confirm that the imported peak sets are not just the per-sample SnapATAC bins — otherwise FRiP computed on a pan-cohort union-peak set will be artifactually zero.

**Preliminary "safe to include in pan-malignant union" ranking** (based on protocol fidelity alone, pending v2 QC):

| Rank | Cohort | Rationale |
|---|---|---|
| 1 | gbm_space | Multiome ARC 2.0.1, NovaSeq S4 100 K rp, Percoll + Digitonin Dounce, RIN pre-screened, strong IRB/consent |
| 2 | gse276177_khan_astro | Multiome ARC 2.0.0, NovaSeq 50 K rp, iodixanol + Digitonin-free Dounce, <1 h to processing |
| 3 | mathewson_lupien | 10x scATAC on NovaSeq, standard Diaz-lab pipeline, large sample pool |
| 4 | gse138794_guo | Same lab/pipeline as #3, older chemistry, smaller counts |
| 5 | guilhamon | 10x scATAC v1.0 on NextSeq, hg19 → need liftover, otherwise clean |
| 6 | tcga_scatac | **Pending methods confirmation** — if 10x v1.1 + ArchR, slots in with the Diaz-lab tier |
| 7 | gse165037 | Methods are solid but **non-droplet**, dbGaP-controlled, no fragment-level data; needs special handling |
| 8 | gbm_tme_atlas_hra004942 | **Pending methods confirmation and QC** — treat as provisional until platform and FRiP/TSSe are verified |

---

## 5. Appendix — verbatim methods per paper

### A1. Guilhamon P et al, *eLife* 2021 (GSE139136, PMID 33427645)

> **Single-cell ATAC-seq** — "The four tumors used were G4218 (primary GBM, IDH wt, male, 64 years), G4250 (primary GBM, IDH wt, male, 73 years), G4275 (primary GBM, IDH wt, female, 52 years), and G4349 (primary GBM, IDH wt, male, 62 years). Fragments of tumor were received fresh from the operating room, and blunt dissected into individual fragments of approximately 0.3–0.7 cm³. Each fragment was placed in 1 mL of freezing media (400 μL of NeuroCult NS-A Basal medium with proliferation supplement (StemCell Technologies; #05751) containing 20 μg/mL rhEGF (Peprotech, AF-100–15), 10 μg/mL bFGF (StemCell Technologies, #78003), and 2 μg/mL heparin (StemCell Technologies, #07980); 500 μL of 25% bovine serum albumin (BSA) (Millipore-Sigma; A9647) in Dulbecco's modified Eagle's medium, and 100 μL DMSO (Millipore-Sigma; D2650) in a 2 mL cryotube, and placed at −80 °C in a CoolCell for at least 24 hr. Samples were then stored at −80 °C until use."
>
> "Cryopreserved primary GBM samples were washed at 1000 RPM for 5 min in phosphate-buffered saline (PBS) to remove DMSO, and then transferred to 1.5 mL tubes. Samples were resuspended in cold ATAC resuspension buffer (10 mM Tris–HCl pH 7.4, 10 mM NaCl, 3 mM MgCl₂, 0.1% NP-40, 0.1% Tween-20, 0.01% Digitonin, 1% BSA in PBS) on ice and dissociated using a wide-bore P1000 pipette tip and vortexing, followed by 10 min of incubation on ice. Cells were spun down at 500 × g for 5 min at 4 °C, washed in the ATAC resuspension buffer, spun down again, and resuspended in ATAC-Tween wash buffer (10 mM Tris–HCl pH 7.4, 10 mM NaCl, 3 mM MgCl₂, 0.1% Tween-20, 1% BSA in PBS), then passed through a cell strainer top FACS tube (Falcon; #38030) to remove debris. Nuclei quality and quantity was evaluated using trypan blue on an Invitrogen Countess II device in duplicate, and a subset of nuclei was spun down in a fresh tube and resuspended in 10× sample dilution buffer. Nuclei were then used for single-cell ATAC-seq library construction using the Chromium Single Cell ATAC Solution v1.0 kit (10× Genomics) on a Chromium controller. Completed libraries were further quality checked for fragment size and distribution using an Agilent TapeStation prior to sequencing."
>
> "Single-cell ATAC-seq samples were sequenced on a NextSeq 500 (Illumina) instrument with 50 bp paired-end reads at the Centre for Health Genomics and Informatics (CHGI) at the University of Calgary. The raw sequencing data was demultiplexed using cellranger-atac mkfastq (Cell Ranger ATAC, version 1.0.0, 10× Genomics). Single-cell ATAC-seq reads were aligned to the hg19 reference genome (hg19, version 1.1.0, 10× Genomics) and quantified using cellranger-atac count function with default parameters (Cell Ranger ATAC, version 1.1.0, 10× Genomics). The resulting data were analyzed using the chromVAR (Schep et al., 2017) and Signac (Stuart et al., 2019) R packages (v1.4.1)."
>
> IRB/consent: "All tissue samples were obtained following informed consent from patients, and all experimental procedures were performed in accordance with the Research Ethics Board at The Hospital for Sick Children (Toronto, Canada), the University of Calgary Ethics Review Board, and the Health Research Ethics Board of Alberta – Cancer Committee (HREBA). Approval to pathological data was obtained from the respective institutional review boards."

Full text source: Europe PMC XML for PMC7847307.

### A2. Wang L et al, *Nature Cancer* 2022 (GSE174554, PMID 36539501) — mathewson_lupien

> **Ethical approval** — "Study protocols and sample use were approved by the University of California, San Francisco (UCSF) Institutional Review Board. All clinical samples were analyzed in a de-identified fashion. All experiments were carried out in conformity to the principles set out in the Declaration of Helsinki as well as the Department of Health and Human Services Belmont Report. Informed written consent was provided by all patients."
>
> **Tumor tissue acquisition** — "We obtained fresh-frozen and FFPE tissue specimens from patients undergoing surgical resection for glioma at UCSF. De-identified samples were provided by the UCSF Neurosurgery Tissue Bank."
>
> **Nuclei isolation** — "For snRNA-seq, nuclei were extracted from frozen tissues following the 'Frakenstein' protocol developed by L. Martelotto, Melbourne, Centre for Cancer Research, Victorian Comprehensive Cancer Centre and available from 10x Genomics … For snATAC-seq, frozen tissues were digested mechanically in a Dounce grinder with 500 µl of lysis buffer (Sigma). The lysate was strained through a 40-μm strainer, pelleted, washed and resuspended in 500 µl nuclei wash buffer (10x Genomics). Nuclei were subsequently purified via centrifugation in a sucrose-based density gradient, pelleted, washed and resuspended in tagmentation buffer (10x Genomics)."
>
> **10x Genomics-based snRNA-seq/snATAC-Seq** — "Single-nucleus capture, reverse transcription, cell lysis and library preparation for snRNA-seq were performed on the 10x Genomics platform as per manufacturer's protocol. Approximately 15,000 nuclei were loaded per capture. For snATAC-seq assay, tagmentation, nuclei capture and library prep were likewise performed via the 10x Genomics platform as per manufacturer's protocol. Sequencing was performed on an Illumina NovaSeq with 10x Genomics recommended parameters."
>
> **snATAC-seq data processing and analysis** — "The CellRanger ATAC software (v.1.1.0) was used for read alignment, deduplication and identifying transposase cut sites … The output matrix of CellRanger was further processed via the snapATAC package … We selected the highest quality barcodes for each case based on two criteria: (1) number of filtered fragments >1,000; and (2) fragments in promoter ratio >0.2 for the case. Clustering was performed using Seurat v.3 SNN-graph clustering via the 'FindClusters' routine, with gene body-accessibility scores generated by the snapATAC package as input. Transcription factor motif frequency deviations from a data-driven background model were calculated via the computeDeviations function in chromVAR (v.1.6.0) with default parameters, using only neoplastic cells as input. Differential motif deviances were computed via a t-test and controlled for multiple hypothesis testing via fdrtool. Differentially accessible regions, peaks and motif enrichments on differential peaks (relative to a genome-wide background) were computed using snapATAC's 'findDAR', 'runMACSForAll' and 'runHomer' respectively, run with default parameters. Heat maps of differential peaks were created in deepTools v.3.4.0."
>
> **CNV analysis** — "CONICSmat (v.1.0) was used to assess the presence/absence of somatic CNVs in 10x snRNA-seq data. We retained CNVs with a CONICSmat likelihood-ratio test <0.05 and a difference in Bayesian Criterion >50 … The presence/absence of somatic CNVs in 10x snATAC-seq data was likewise estimated with CONICSmat. Here, the gene activity of cells generated by snapATAC (v.1.0.0) was used as input to perform CNV analysis."

Full text source: Europe PMC XML for PMC9767870.

### A3. Lin, Chen, Li, Chen, Fang … Qu Kun, *Nature Neuroscience* 2026 (HRA004942) — gbm_tme_atlas_hra004942

**Methods not accessible via WebFetch.** The *Nature Neuroscience* article at https://www.nature.com/articles/s41593-026-02265-5 redirects to IDP auth, and the NGDC project page at https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004942 returns bare metadata only. NGDC project record confirms: title "Spatial and single-cell characterization of human glioblastoma tumor microenvironment reveals malignant cellular communities"; PI Qu Kun (qukun@ustc.edu.cn); organization University of Science and Technology of China; data accessibility open; release 2026-02-17; submitted 2023-06-26. The patient IDs in the local audit (P62, P64, P77, P78, P79, P80, P83, P84, P92, P98, P101) number 11; this matches the "11-patient" project scope.

**To fill in**:
- Download the Nature Neuroscience article's Methods PDF through an institutional subscription.
- The patient-level cell counts in the local audit (158–761 cells/patient) are striking — if the authors used **DNBelab C4** or a comparable low-throughput platform, this is the biggest cross-cohort methods outlier. Verify.

### A4. Wang L … Diaz AA, *Cancer Discovery* 2019 (GSE138794, PMID 31554641) — gse138794_guo

**Methods are not deposited in PMC full text** (abstract and SOFT metadata only). GEO SOFT metadata gives the per-sample protocol fields (reproduced verbatim here):

> **Overall design (series)** — "We performed single-cell RNA sequencing (scRNA-seq), single-nuclei RNA sequencing (snRNA-seq), single-cell assay for transposase-accessible chromatin using sequencing (scATAC-seq), and whole-exome DNA sequencing (exome-seq) of specimens from untreated human gliomas … From the scATAC-seq we elucidated cell-type specific cis-regulatory grammars and associated transcription factors."
>
> **Sample extract protocol (snATAC samples)** — "Nuclei were extracted from frozen tissues following the 'Frakenstein' protocol developed by Luciano Martelotto, Ph.D., Melbourne, Centre for Cancer Research, Victorian Comprehensive Cancer Centre, and available from 10X Genomics."
>
> "For fresh tissues, 10.2 μL of live cells, at a concentration of 1700 live-cells/μL, were loaded into the 10X Chromium Single Cell capture chip. For frozen tissues, ~15,000 nuclei were loaded per capture. Single-cell/nucleus capture, reverse transcription, cell lysis, and library preparation were performed per manufacturer's protocol. Sequencing was performed on an Illumina NovaSeq using a paired-end 100 bp protocol."
>
> **Data processing** — "The CellRanger ATAC software (version 1.1.0) was used for read alignment, deduplication, and identifying transposase cut sites. Genome build: hg38. Processed data files format and content: Matrix file (mtx) includes sparse matrix of read counts, unique cell barcode files (tsv), and feature file (tsv)."
>
> **Series note** — "Raw data not available due to privacy concerns" (dbGaP-controlled).
>
> **Instrument** — Illumina NovaSeq 6000 (per Sample_instrument_model fields).

Case count per the series summary: **28 gliomas for RNA, 8 cases for scATAC**. Local audit retained 3 cases for the pan-malignant matrix: SF11215, SF11956, SF11979. IRB/consent presumably matches mathewson_lupien (same lab, UCSF Neurosurgery Tissue Bank), but the specific Cancer Discovery 2019 methods PDF must be read for the exact wording.

**To fill in from manual pull** (Cancer Discovery full text): per-case IDH status / primary-vs-recurrent / grade; exact loading target for scATAC (frozen vs fresh ratio in the 8 cases); the full IRB language.

### A5. Terekhanova NM … Chang HY/Corces MR, *Science* 2024 (TCGA-ATAC-Seq-2024, PMID 39236169) — tcga_scatac

**Methods not accessible via WebFetch.** Both https://www.science.org/doi/10.1126/science.adk9217 and the GDC publication page at https://gdc.cancer.gov/about-data/publications/TCGA-ATAC-Seq-2024 are JS/paywall-walled. From abstract: single-cell chromatin accessibility across 8 TCGA tumor types, 74 individual samples, with cancer cells / tumor-infiltrating immune / stromal separated by scATAC. GBM is one of the 8 types; the local audit has 10 GBM-labeled samples with the pattern `scATAC_GBMx_<UUID>_X###_S##_B#_T#` (classic TCGA barcoding).

**To fill in from manual pull**: Science supplementary methods PDF has GBM-specific sample acquisition (likely 10x scATAC v1.1 kit based on the Chang/Corces group's prior work), the exact nuclei isolation protocol (likely Omni-ATAC-derived), ArchR version, peak-calling parameters (likely iterative overlap, 500 bp fixed-width peaks, MACS2 under the hood). Patients are TCGA-GBM cases — IDH status is in GDC clinical files, not the paper.

### A6. de Jong G, … Saraswat M, … Bayraktar OA, *bioRxiv* 2025 (E-MTAB-17183, DOI 10.1101/2025.05.13.653495) — gbm_space

> **Human Subjects** — "Patients with suspected GB were identified pre-operatively and consented for entry into the study. Surgery was performed at Cambridge University Hospitals NHS Foundation Trust. Written and informed consent was obtained in accordance with the guidelines in The Declaration of Helsinki 2000. Ethical approval for the use of these tissues was obtained from the Cambridge Local Research Ethics Committee (REC 18/EE/0172). All patients underwent 5-ALA guided tumour resection as per local protocols. During tumour debulking, regions of high fluorescence were identified, their spatial location recorded and the tissue samples were collected for this study. Matched whole blood was taken during surgery for germline characterisation."
>
> **Serial tissue sectioning** — "Tissue was sampled from multiple sites of each GB tumour, targeting superior, anterior, posterior and middle regions where possible. Each tissue sample was immediately washed in saline buffer and embedded in OCT medium (Scigen OCT Compound, #4586) using a dry ice-cooled bath of isopentane at −75 °C. OCT-embedded samples were sectioned using a cryostat (Leica CX3050S). The fresh frozen tissue blocks were trimmed until the tissue surface was fully exposed. Two to three 10 µm thick sections were collected to check RNA integrity … Only samples with RNA integrity number (RIN) values >7 were used for omic profiling … Sectioning for single nuclei isolation: A series of 50 µm thick sections, totalling 500 to 900 µm thick volume depending on the size of each tumour block, were collected in pre-chilled homogenization glass tubes and kept on dry ice until processing."
>
> **Single-nuclei extraction** — "Nuclei were extracted from fresh frozen tissue sections that were homogenised using a glass Dounce homogenizer (Sigma) in nuclei isolation buffer (3 mM MgCl₂, 10 mM NaCl, 10 mM Tris (buffer pH 7.4), 1 mM DTT, 0.1% Tween-20, 0.1% Nonidet P40, 1% BSA and 0.01% Digitonin) in the presence of Protector RNase Inhibitor (Roche) at 0.2 U/μl. Tissue was homogenised using 10 strokes with pestle A and then 10 strokes with pestle B. Nuclei were then filtered through a 40 μM filter, collected at 500 RCF and resuspended in 0.25 ml of storage buffer (PBS containing 1% BSA and Protector RNase Inhibitor (Roche) 1 U/μl). An aliquot of the nuclei suspension was incubated with Trypan Blue (Gibco 15250061) for counting and purified from debris using a Percoll gradient. The cleaned nuclei suspension was stained with Trypan blue and counted."
>
> **10x Genomics Chromium GEX and ATAC library preparation and sequencing** — "For the snRNA-seq experiments, two to three 10x reactions were prepared per tumour site and loaded onto the 10X chromium controller according to the manufacturer's protocol for the Chromium Next GEM Single Cell Multiome ATAC + Gene Expression assay. Post-GEM-RT cleanup, cDNA amplification and 3′ gene expression library construction were carried out as per the Chromium Single Cell 3' Reagent Kits v3 User Guide, to obtain between 5000-10,000 nuclei per reaction. Libraries were paired end-sequenced on a NovaSeq 6000 System (Illumina) using the Novaseq S4 Flowcell, targeting a minimum coverage of 100,000 read pairs per nuclei per modality. The following sequencing formats were employed for GEX and ATAC respectively: — GEX: Read 1: 28 cycles; i7 Index: 10 cycles; i5 Index: 10 cycles; Read 2: 90 cycles — ATAC: Read 1N: 50 cycles; i7 Index: 8 cycles; i5 Index: 24 cycles; Read 2N: 49 cycles."
>
> **Single nuclei multiome data processing and quality control** — "We aligned reads from each snRNA-seq and ATAC-seq library to a custom-made genome consisting of 10X Genomics' GRCh38 3.0.0 pre-mRNA reference genome and 10X Genomics Cell Ranger ARC 2.0.1 ATAC genome. To perform quantification and initial quality control, we used the default parameters in the Cell Ranger ARC software (v2.0.1; 10X Genomics). This was followed by CellBender (v0.2.0), which was applied to the Cell Ranger output to correct for background noise and identify empty droplets … Quality control of the RNA data was performed using Cell Ranger ARC filtered count matrices with Scanpy (v1.9.3). This involved removing nuclei with total gene counts <500, total counts UMI <1000, and nuclei with >10% of reads mapping to mitochondrial content. We then applied Scrublet (v0.2.3) on each library individually and filtered our data based on a two-step method adapted from previously described median absolute deviation (MAD) thresholding … ATAC data were processed according to the ArchR workflow (minimum fragments >1000 and >4 transcriptional start sites). These barcodes were then filtered by the list of barcodes passing RNA quality control filters … Finally, the RNA data was subset by barcodes passing ATAC filters to ensure symmetry between modalities for downstream applications."

Full text source: bioRxiv JATS XML at https://www.biorxiv.org/content/early/2025/05/14/2025.05.13.653495.source.xml.

### A7. Sojka C … Sloan SA, *Nature Cell Biology* 2025 (GSE276177, PMID 39779941) — gse276177_khan_astro

**Note on slug:** the pipeline slug `gse276177_khan_astro` is a working alias from an earlier revision; the published paper is **Sojka et al. Nat Cell Biol 2025, DOI 10.1038/s41556-024-01583-9**, titled "Mapping the developmental trajectory of human astrocytes reveals divergence in glioblastoma." The GEO working title "Divergence from the human astrocyte developmental trajectory in glioblastoma" matches the paper's subject but not its final title. Not a Khan-authored paper.

GEO SOFT (verbatim, same text across all 10 Multiome samples):

> **Sample_extract_protocol_ch1** — "GBM nuclei extraction: ~20 mg of tissue was dissociated with a 2-ml Dounce homogenizer in homogenization buffer (0.26 M sucrose, 0.03 M KCl, 0.01 M MgCl₂, 0.02 M Tricine-KOH pH 7.8, 0.001 M DTT, 0.5 mM Spermidine, 0.15 mM Spermine, 0.3% NP40, and cOmplete Protease inhibitor). This was followed by filtering through a 40 µm Flowmi cell strainer, then a 20 µm bucket-style cell strainer, and centrifugation for 10 min at 600 r.c.f. After the majority of the supernatant was carefully removed, the pellet was resuspended in homogenization buffer and mixed with an equal volume of 50% iodixanol solution to make a final concentration of 25% iodixanol. Next, a 30% iodixanol solution, followed by a 40% iodixanol solution, was layered under the 25% mixture and centrifuged for 20 min at 3000 r.c.f without the centrifuge brake. Post-centrifugation, a thin white nuclei band was carefully collected from the interface of the 30% and 40% iodixanol solutions. Nuclei underwent 1-2 wash steps to remove any additional debris by gently mixing nuclei in ATAC-RSB-Tween buffer (0.01 M Tris-HCl pH 7.5, 0.01 M NaCl, 0.003 M MgCl₂, 0.1% Tween-20) and centrifuging for 10 min at 600 r.c.f."
>
> "Libraries were generated using the 10x Genomics Chromium Single Cell Multiome ATAC + Gene Expression kit following the manufacturer's instructions. Per sample, **16,100 nuclei were resuspended in 1x diluted nuclei buffer (10x Genomics) with 2% BSA (Sigma) with a capture target of 10,000 nuclei**."
>
> **Sample_data_processing** — "The demultiplexing, barcoded processing, gene counting and aggregation were made using the 10x Genomics Cell Ranger ARC (cellranger-arc-2.0.0). Assembly: hg38."
>
> **Sample_instrument_model** — Illumina NovaSeq 6000.

From the Sojka 2025 Nature Cell Biology full-text extract (via WebFetch summary of PMC12210326): IRB — "obtained in compliance with policies outlined by the Emory School of Medicine IRB office"; tumor and margin were "immediately deposited in 4 °C Hibernate-A medium … and prepared for tissue dissociation within 1 hr post-resection"; sequencing "at a target depth of at least 50,000 read-pairs per nucleus using 2 × 150-bp reads on an Illumina Novaseq 6000 instrument." Cases: GBM20 (M, IDHwt), GBM25 (F, IDHwt), GBM38 (M, IDHwt).

Full text source: GEO SOFT (verbatim) + Sojka et al. 2025 Nat Cell Biol methods summary via PMC12210326.

### A8. Raviram R, Chen C, Ren B, GSE165037 (2021, no linked PMID in series) — gse165037

Methods are from GEO SOFT only (verbatim, same text across all GBM samples):

> **Sample_extract_protocol_ch1 (snATAC)** — "Nuclei were isolated from ground frozen tumor tissue and permeabilized using NPB (5% BSA, 0.2% IGEPAL-CA630, cOmplete Protease Inhibitors, 1 mM DTT in PBS)."
>
> "Nuclei were tagmented with Tn5 in 96 wells to introduce a first barcode. After pooling, 20-25 nuclei were sorted into 768 wells and a second barcode was introduced by PCR. Libraries were purified and sequenced on a HiSeq 2500, HiSeq 4000 or NextSeq 500 (Illumina)."
>
> **Sample_extract_protocol_ch1 (snRNA)** — "Ground frozen tumor tissue was resuspended in nuclei buffer (2% BSA, 0.2% Triton-X, cOmplete Protease Inhibitors (Roche), 1 mM DTT, 0.2 U/µl RNAsin (Promega) in PBS). After spin down, nuclei were resuspended in sort buffer, filtered and sorted using a SH800 sorter (Sony). Nuclei were loaded on a 10x Chromium controller and libraries were generated using the Chromium Single Cell 3' v2 Library kit (10x Genomics) according to manufacturer descriptions and sequenced on a HiSeq 4000 or NextSeq 500 (Illumina)."
>
> **Sample_data_processing** — "Reads from snRNA-seq were mapped to the hg38 genome using STAR. Clustering analysis was performed using Seurat V3 on cells that had greater than 200 reads and less than 5% mitochondrial reads."
>
> "Reads from snATAC-seq were mapped to the hg38 genome using bowtie2. Snaptools was used to generate binary count matrices in 5kb genomic bins. Clustering analysis was performed using SnapATAC for cells with greater than 1000 reads and percentage of reads in promoters greater than 20%."
>
> **Overall design note** — "Submitter states that raw data will be submitted to dbGaP due to patient privacy concerns."

Cases: GBM1 (IDH1-mut), GBM4 (IDHwt), GBM9 (IDHwt), GBM11 (IDH1-mut), GBM12 (IDHwt), Non-tumor brain. Contributors: Ramya Raviram, Clark Chen (surgeon, UCSD), Bing Ren (PI, UCSD). The combinatorial-indexing protocol is in the Preissl/Ren lab style (Cusanovich / sci-ATAC-seq lineage).

Full text source: GEO SOFT metadata (`/private/tmp/.../geo/GSE165037.soft`).

---

## 6. Follow-ups

**Three methods pulls still needed** (behind paywall / JS; can be done manually by the user in a few minutes each):

1. **Nature Neuroscience 2026, Qu Kun / HRA004942** — download the Methods PDF through institutional access and confirm the sequencing platform (strongly suspect DNBelab C4 given the ~500-cell/patient yield in the audit). Direct URL: https://www.nature.com/articles/s41593-026-02265-5.
2. **Science 2024, Terekhanova et al.** — pull the Supplementary Materials PDF from https://www.science.org/doi/10.1126/science.adk9217. Confirm whether the GBM arm used 10x scATAC v1.1 vs Multiome vs a mixed set, and whether the ArchR peak set is 500 bp fixed-width.
3. **Cancer Discovery 2019, Wang/Diaz** — the main paper's PDF (https://aacrjournals.org/cancerdiscovery/article/9/12/1708/42245) likely has the explicit per-case IDH status / primary vs recurrent, and may add detail beyond the SOFT metadata (which only states "frozen tissues ~15,000 nuclei loaded"). The lab is the same as Mathewson_lupien, and the SOFT already names CellRanger ATAC v1.1.0 / NovaSeq 6000 / PE100 / hg38.

**Audit CSV v2 (TSSe + FRiP + fragment counts)** — once pulled from EC2 to `/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results/qc/outputs/01_per_patient_audit_v2.csv`, re-run the "preliminary good vs bad preps" section 4 with real numbers in place of the hypotheses.
