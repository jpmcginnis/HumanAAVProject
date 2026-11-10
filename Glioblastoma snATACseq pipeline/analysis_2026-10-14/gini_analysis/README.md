# Gini-index specificity analysis across cell types

Method (standard, as cited in e.g. CATLAS methodology):
> We ranked regions based on the Gini index of their accessibility profile
> over all cell types to obtain the highest and most specific peaks per cell
> type. Using the Gini index only provides one value per region, so we
> assigned that value to the cell type with the highest peak value and gave
> a zero score to all the other cell types in that region.

## Implementation

- Loaded the union pseudobulk matrix (`enhancer_candidate_matrix.parquet`) with
  one column per `(cohort, patient, cell_type)` tuple.
- Pooled per-cell-type: for each of the 10 cell types below, summed
  raw counts across all its patient columns, then CPM-normalized per cell type.
- Cell types used: astrocyte, microglia, oligodendrocyte, OPC, endothelial, neuron, GABA_neuron, T_cell, TAM, malignant.
- `unassigned` columns excluded. Malignant pooled (not split by Neftel subtype).
- For each peak, computed Gini across the 10-element CPM vector.
  All-zero peaks assigned Gini = 0.
- `assigned_ct` = the cell type with the highest CPM for that peak.
- Full matrix: `gini_scores_per_peak.parquet`
  - Columns: peak_id, chrom, start, end, gini_score, assigned_ct, max_cpm,
    cpm_<ct> (10 cols), gini_<ct> (10 cols — assigned score for winner, 0 for others).
- Per-cell-type top-500 tables: `gini_top500_per_celltype/gini_top500_<ct>.tsv`
  - Requires `max_cpm >= 1.0` for inclusion (noise filter).
  - Ranked by gini_score, descending.
- Summary: `gini_summary.tsv`

## Notes / caveats

- Gini is a specificity metric only. Highly ranked peaks are those that stand
  out in ONE cell type's accessibility profile relative to the other nine.
- Pooled-malignant assignment means Neftel subtype specificity is NOT resolved
  in this run. For malignant subtypes, re-run with the Neftel subtypes
  (malignant_AC, malignant_MES, malignant_NPC, malignant_OPC) as separate columns.
- CNV-amplified-locus caveat (per Part 1 methods disclaimer): non-malignant
  candidates near chr7p11.2 EGFR, chr12q15 MDM2, chr8q24 MYC, chr9p21.3 CDKN2A
  are subject to ambient-fragment contamination and should be flagged. The
  Gini metric is sensitive to this because contamination inflates the
  non-malignant cell type's CPM at the amplicon.
- The TAM↔microglia distinction (one-direction Gini only promotes the winner)
  means a peak high in both will land in whichever is higher; downstream
  analysis of myeloid should look at both.
