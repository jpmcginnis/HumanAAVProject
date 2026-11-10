# ArchR-baseline differential accessibility

Classic ArchR `getMarkerFeatures`-style filter: log2FC > 2.0 (fold > 4.0x)
AND BH-FDR < 0.05, applied per cell type.

## Two variants (both produced here)

### (A) Pooled-background (`per_celltype_pooled_bg/`) — ArchR default
- Background = all patient-pseudobulks from OTHER cell types, pooled.
- log2FC = log2((mean target CPM + 0.1) / (mean pooled background CPM + 0.1)).
- Welch's t-test on per-patient CPMs: target cell type patients vs pooled other-cell-type patients.
- BH-FDR correction across all 1,958,753 peaks per cell type.

### (B) Max-of-any-other-cell-type (`per_celltype_max_bg/`) — strict "any other" interpretation
- Background = the SINGLE highest other-cell-type mean CPM per peak.
- log2FC = log2((mean target CPM + 0.1) / (max OTHER cell type mean CPM + 0.1)).
- This is the "at least 4x any other cell type" filter.
- Statistical test (Welch) still against pooled background for computational tractability;
  the strictness here comes from the log2FC definition, not a separate test.
- Each survivor row carries the `top_competitor_ct` — which other cell type was the
  highest competitor at that peak.

## Interpretation

- log2FC > 2.0 means fold change > 4.0x. 'At least twice' literally
  would be log2FC > 1 (fold > 2x) — if that's intended, re-run with LOG2FC_THR=1.
- The pooled-background definition is LAX because background means dilute per-peak
  highs in one or two competitor cell types; many peaks that are high in one other
  cell type still have low pooled-bg mean.
- The max-background definition is STRICT because any single competitor at or above
  target disqualifies the peak.
- A peak surviving BOTH filters is cell-type-exclusive by a very strong definition.

## Outputs

- `archr_combined_summary.tsv` — per-cell-type survivor counts under both definitions
- `per_celltype_pooled_bg/archr_pooledbg_<ct>.tsv` — survivor peaks under pooled-bg, ranked by log2FC
- `per_celltype_max_bg/archr_maxbg_<ct>.tsv` — survivor peaks under max-bg, ranked by log2FC

## Caveats

- Amplicon-contamination caveat (per Part 1 methods disclaimer): non-malignant
  survivors near chr7p11.2 EGFR, chr12q15 MDM2, chr8q24 MYC, chr9p21.3 CDKN2A
  are subject to ambient-fragment contamination; apply the `amplicon_adjacent` flag.
- TAM↔microglia share myeloid baseline — many TAM survivors will have microglia as
  the top_competitor_ct and vice versa; this is a feature of the activation-state
  framing, not a bug.
- The Welch t-test assumes approximately normal residuals, which is marginal for
  per-patient CPMs; log1p-transformed CPMs give more defensible normality — a
  more conservative FDR could be obtained by log1p-transforming first.
