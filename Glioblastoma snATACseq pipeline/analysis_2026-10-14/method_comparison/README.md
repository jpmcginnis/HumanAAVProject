# Gini vs PeakRankR method comparison

Compares the top-500 per cell type produced by two different enhancer
ranking methodologies on this GBM atlas.

## Why they differ

- **Gini index** is a pure cell-type-specificity measure: a peak with signal
  concentrated in ONE cell type gets a high Gini (~0.9); a peak uniformly
  accessible across all cell types gets a Gini near 0. It favors cell-type-
  exclusive enhancers, including "dark matter" peaks seen in only 1-2
  lineages.
- **PeakRankR** = Specificity + Sensitivity + Magnitude (equally weighted
  ranks). Because Sensitivity rewards universal detection (robust peak call
  across many cell types), the top-500 are dominated by peaks accessible in
  ALL 10 cell types with moderate specificity enrichment in the target.
  It favors strong-signal, broadly-detected enhancers that are additionally
  enriched in the target cell type.

The two methods answer different questions: Gini finds the most private /
cell-type-selective enhancers; PeakRankR finds the most robust strong-signal
enhancers with cell-type enrichment on top.

## Outputs

- `overlap_summary.tsv` — per-cell-type n_overlap, n_gini_only, n_peakrankr_only,
  Jaccard index, median feature values per partition.
- `per_cell_type_overlap/overlap_<ct>.tsv` — full peak_id × method-flag table
  per cell type with feature values.
- `gini_vs_peakrankr_overlap.png` — stacked-bar visualization + Jaccard column.
