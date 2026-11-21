# PeakRankR analysis

Method (Cell Genomics 5, 100879, June 2025):
> PeakRankR employs a linear model approach to calculate the combined effects
> of positive features associated with a given peak in a cell type, which are
> then used to rank the peaks.
>
> PeakRankR_Score = ω_Specificity × Specificity + ω_Sensitivity × Sensitivity + ω_Magnitude × Magnitude
>
> By default, all weights are set to 1, indicating equal importance for each feature.

## Feature implementation on this atlas

Per-cell-type pseudobulks were built by summing raw counts across all patient
columns of a cell type, then CPM-normalizing per cell type (same aggregation as
the Gini analysis). Features per peak × target cell type c:

1. **Specificity (log2 fold change over background)**: `log2((CPM_c + 0.1) / (mean_other_cells_CPM + 0.1))`
   — higher value means target cell type dominates signal.
2. **Sensitivity (count of cell types with signal)**: n cell types with CPM ≥ 1.0
   at this peak. The methodology defines this as "number of cell types that
   have a peak at the coordinates of the enhancer" — a quality / robustness
   measure (how many lineages reproduce the peak call), not an off-target
   measure. We used a CPM detection floor as a proxy for peak-callability.
3. **Magnitude (CPM of target cell type)**: pseudobulk CPM of cell type c at
   this peak.

PeakRankR_Score_c = w_S × rank(Specificity) + w_Se × rank(Sensitivity) + w_M × rank(Magnitude)
with default weights = 1. Ranks use `method='average'`. Peaks are then ranked
within each cell type by this combined score, descending.

## NOT implemented in this run (shape-based enhanced features)

The enhanced PeakRankR adds Skew, Kurtosis, and Modality features computed on
read-pileup shape at each enhancer coordinate. These require **per-cell-type
bigWig files**, which would need to be built from BAM alignments across all 8
cohorts × 10 cell types × many patients. We did not generate those bigWigs in
this run. These shape features are a potential future-work extension and would
be most impactful for distinguishing clean enhancer peaks from noisy / artifact
pileups among the final top candidates.

## Cell types scored

astrocyte, microglia, oligodendrocyte, OPC, endothelial, neuron, GABA_neuron, T_cell, TAM, malignant

## Weights

Specificity=1.0, Sensitivity=1.0,
Magnitude=1.0 (default = equal).

## Outputs

- `peakrankr_scores_per_peak.parquet` — full table: all peaks × (coords, sensitivity,
  per-cell-type CPM / specificity log2fc / magnitude / score / rank)
- `peakrankr_rank_per_celltype_slim.tsv.gz` — slim version with just ranks + scores
- `peakrankr_top500_per_celltype/peakrankr_top500_<ct>.tsv` — top-500 per cell type,
  ranked by PeakRankR_Score; requires `magnitude_target_cpm ≥ 1` (noise floor).
- `peakrankr_summary.tsv` — per-cell-type summary (median feature values at the top-500)

## Caveats

- Amplicon-contamination caveat (per Part 1 methods disclaimer): non-malignant
  candidates near chr7p11.2 EGFR, chr12q15 MDM2, chr8q24 MYC, chr9p21.3 CDKN2A
  are subject to ambient-fragment contamination and should be flagged. Specificity
  computed here will LOOK inflated for non-malignant cell types near these loci
  because contamination raises target_cpm.
- TAM is scored against Microglia (and others) as background. A peak high in
  both TAM and microglia will get a low specificity score for TAM — this reflects
  the activation-state framing (TAM = microglia variant) correctly but may hide
  peaks that are robust across both myeloid states.
