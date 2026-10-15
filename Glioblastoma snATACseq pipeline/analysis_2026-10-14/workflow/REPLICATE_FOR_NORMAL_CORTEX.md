# Replicating the pipeline for normal human cortex (v5)

A protocol for someone familiar with scATAC/snATAC-seq who wants to run the same
cross-cohort meta-analysis for **normal (non-tumor) human cortex** instead of GBM,
following the v5 architecture: per-cohort union peak reference + overlap-projected
pseudobulk + Daigle Z scoring.

Shipped at <https://github.com/jpmcginnis/HumanAAVProject> → `Glioblastoma snATACseq pipeline/`.

## Why replicate for normal cortex at all?

1. **Reference control** for the GBM atlas: identify which peaks are tumor-specific
   vs ones that would also fire in healthy brain (bad AAV cloning targets if you
   want tumor-restricted expression).
2. **Normal-cortex-specific AAV enhancers**: cell-type-specific enhancers that fire
   in neurons/astrocytes/microglia/OPC/oligo in healthy brain — the Allen BRAIN
   Armamentarium's own use case (Mich et al. 2021 *Nat Methods*, Hunker striatum).
3. **Cross-condition panels**: enhancers that fire in both normal and GBM
   (pan-brain-cell-type), useful for broad biosensors or disease-agnostic reporters.

## What's the SAME as v5 GBM

- The **union peak construction logic**: concatenate each cohort's own called peaks
  (or author-published BEDs for matrix-mode cohorts) + CATLAS, merge with
  `bedtools merge -d -100` (prevents transitive chaining of densely-packed 500 bp
  peaks). `scripts/cluster_summits_v5.py` (clusters per-sample summits with a 400 bp
  window before merge) carries over verbatim.
- The **per-cohort re-quantification** chain: snapatac2 `pp.make_peak_matrix` for
  fragments-mode cohorts; sparse matrix multiplication through an overlap map for
  matrix-mode cohorts (`cells × cohort_peaks @ cohort_peaks × v5_peaks` = `cells × v5_peaks`).
- The **Daigle Z scoring** per cell type: target column mean vs background pooled
  across all other `(cohort × cell_type)` columns, with CPM normalization.
- The **replication filter**: `n_target_cohorts ≥ 2` (count of cohorts with nonzero
  signal at this peak × cell_type).
- The **peak-level provenance table**: every union peak traced back to the
  cohort(s) that contributed and their original coordinates.
- The **fragment-level pileup helper** (`scripts/pileup_region.py`) for
  interrogating any merged interval against raw per-cohort fragments.
- The **distance annotation**: nearest refGene TSS, signed distance, four-way
  category (`promoter` ≤2 kb, `proximal` ≤10 kb, `distal` ≤100 kb, `far-distal` >100 kb).
- The **AAV-focused filtering**: top-N distal only (`distance_category != "promoter"`)
  for the shortlist BEDs/CSVs. The 2 kb TSS floor is non-negotiable.

## What CHANGES

### Drop

- **CNV / malignant pool concept.** No `malignant_cnv` column, no
  `malignant_unresolved` cell type, no pan-malignant report. Normal brain has no
  neoplastic axis.
- **chr7/10 CNV gates.** Not applicable.
- **inferCNVpy** step from the ingest chain.
- **IDH-WT subset selection** (irrelevant without tumor samples).

### Change

- **Cohort list.** Replace GBM cohorts with normal-cortex sources: AMBA/Allen
  M1/MTG snATAC, CATLAS normal cortex (already in reference), MouseLight,
  ENCODE adult human cortex, PsychENCODE2 ATAC subset. Keep the per-cohort
  canonical `<firstauthor_year>` slug convention.
- **"Pan-X" reports.** Instead of pan-malignant and pan-myeloid, build:
  - pan-neuronal (all neuronal subtypes)
  - pan-glial (astrocyte ∪ oligo ∪ OPC)
  - pan-myeloid (microglia; no TAM in normal cortex)
- **Cell-type labels.** Expect finer-grained labels in normal cortex datasets
  (L2/3, L4, L5-IT, L5-ET, L6, SST, PV, VIP, LAMP5, etc.). Decide upfront
  whether to roll these up to broad classes (GABA_neuron, glutamate_neuron) or
  score per-subtype. Daigle Z works at any granularity; per-subtype gives more
  targeted AAV candidates but requires enough cohorts to carry each subtype.
- **CATLAS contribution is now larger-share.** In v5 GBM, CATLAS contributed
  25.8% of the union because tumor cell types pulled in many de-novo peaks. In
  normal cortex, CATLAS will cover more of the signal — expect ≥50% of the union
  to be CATLAS, meaning the "de-novo gain" story is smaller than in GBM.

### Watch out for

- **Platform / batch effects** still a concern: sci-ATAC vs 10x cohorts need
  separate QC (we saw this in GBM Raviram 2023). Replication filter protects
  against single-cohort noise, but be explicit about which cohorts contribute
  per cell type in your methods.
- **Fresh vs postmortem** is a bigger issue in normal cortex than in GBM
  (shorter PMI for GBM resections). Expect microglia activation signatures to
  diverge between fresh and postmortem cohorts.
- **Reactive glia.** Even "normal" cortex contains some reactive microglia and
  astrocytes near the sample edges. Not a reason to exclude, but acknowledge.

## Pipeline order

1. **Ingest** per cohort → per-sample fragments h5ad (snapatac2) OR matrix-mode
   h5ad with cells × author-peaks (for cohorts with only published count
   matrices).
2. **Cell-type label.** Use the dataset's native labels; harmonize to a common
   schema (`cell_type`, `patient_id`, `sample_id`, `cohort`).
3. **Call per-sample peaks** (macs3 CLI on raw fragments, 500 bp summit-centered,
   ArchR flags: `-q 0.01 --shift -75 --extsize 150 --nolambda --call-summits`).
   For >150M-fragment samples, use per-chromosome splitting to avoid the known
   MACS3 3.0 memory-view bug.
4. **Cluster summits per cohort** (`cluster_summits_v5.py`): per-sample 500 bp
   peaks → cluster midpoints within 400 bp → 500 bp window around median summit.
   Produces a non-overlapping, non-chaining per-cohort BED.
5. **Build union reference**: concatenate all cohort BEDs + CATLAS → sort →
   `bedtools merge -d -100`. Write `union_peak_reference.bed`.
6. **Peak provenance table**: for each union peak, record which cohorts
   contributed via `bedtools intersect`. Write `peak_provenance.parquet`.
7. **Re-quantify** per cohort → pseudobulk per `(cohort, cell_type)` →
   `enhancer_candidate_matrix.parquet` with rows = v5 peaks, columns =
   `<slug>__<cell_type>`.
8. **Score**: CPM-normalize per column; Daigle Z per cell type; replication count.
9. **Annotate** every peak with nearest gene + distance category.
10. **Build reports** (`build_v5_reports.py`): per-cell-type + pan-X + combined
    top-50 + HTML gallery.

## Cost + wall-clock expectations

Baselined off the v5 GBM run on an r6i.4xlarge (16 vCPU, 128 GB RAM):

- Peak-calling (14 samples, 2 fragments-mode cohorts): ~6-8 hr with the 2 GB
  compressed → per-chr threshold.
- Union + provenance build: ~15 min.
- Re-quantification (6 projection cohorts + 2 fragments-mode): ~1-2 hr.
- Scoring + reports: ~5 min.
- Fragment-level pileup for interrogating one peak: a few seconds per cohort.

Normal cortex has roughly similar scale (~1-2M cells across 5-8 cohorts), so
expect ~10-15 hr end-to-end, dominated by peak-calling on the largest fragment
files.

## Known traps I'd want next-person to avoid

1. **Don't let macs3 write to /tmp on a small root partition.** Set
   `TMPDIR=/data/<project>/tmp` before any macs3 invocation. In v5 GBM we lost
   several hours when / filled silently and macs3 produced empty peak files.
2. **Don't trust a success count that comes from grepping log lines.** Verify by
   `wc -l` on each output BED. Silent 0-byte failures look identical to successes
   in a line count.
3. **Catch BaseException**, not just Exception. snapatac2 raises pyo3
   `PanicException` which does not inherit from Exception; a vanilla
   `except Exception` lets the whole run die on one bad sample.
4. **Barcode suffix auto-detection is required.** Cohort-level `catlas_quantified`
   h5ads have barcodes with `-<sample_idx>` appended; per-sample h5ads don't.
   You need a routine like `detect_sample_suffix()` in `requant_v5b.py` or every
   cell lands in "unassigned."
5. **Don't `bedtools merge` 500 bp peaks with default `-d 0`.** Dense peaks chain
   into chromosome-arm intervals. Use `-d -100` (require 100 bp overlap) for the
   cross-cohort union; use `-d -400` for within-cohort if you're not using the
   midpoint-cluster approach.
6. **Don't try to read Greenwald-style 68 GB filtered-peaks h5ads via anndata
   backed mode.** Row-slice reads are slow (~5 MB/s) because the HDF5 chunks are
   gzip-compressed. Use direct h5py reads of `X.data` + `X.indices` with a 4 GB
   chunk cache and accumulate col-sums via `np.add.at`. Takes ~15 min instead
   of ~8 hours.
