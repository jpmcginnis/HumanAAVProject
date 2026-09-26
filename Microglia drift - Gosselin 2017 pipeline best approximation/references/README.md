# References

Gene lists derived directly from Gosselin 2017 supplementary tables.

## `gosselin_conserved_signature_477.csv`

The **477-gene conserved (mouse-human orthologous) microglia identity
signature.** This is the gene set Gosselin's headline "33 % dropped
> 2-fold in vitro" is computed over.

- **Source:** Gosselin 2017 supplementary Table S5, sheet "MG gene
  sig - similar expressed".
- **Column:** `human_gene` (single column, HGNC symbols).
- **Provenance:** Gosselin defined a mouse microglia signature of 900
  genes that were > 10-fold higher in microglia than in other tissue
  macrophages. Of those 900, 477 had human orthologs expressed at
  comparable levels in his human microglia dataset. This 477-gene set
  is the "conserved" signature.

## `gosselin_881_signature.csv`

The **881-gene human microglia-vs-cortex signature.** Genes expressed
> 10-fold higher in Gosselin's human microglia than in his human
cortex tissue (FDR<0.05).

- **Source:** Gosselin 2017 supplementary Table S2, sheet "Human
  microglia gene signature".
- **Column:** `gene` (single column, HGNC symbols).
- **Provenance:** Gosselin's primary microglia identity signature at
  the human level. The 477 conserved signature is *not* a subset of
  this — the two sets are defined against different reference
  populations (cortex tissue vs other tissue macrophages) and overlap
  in only 70 genes.

## `gosselin_881_per_gene_logfc_invitro_vs_exvivo.csv`

Per-gene log2 fold change (7-day in-vitro vs ex-vivo) computed from
Gosselin's own published TPM values, restricted to the 881-gene human
microglia signature. Used for the optional per-gene concordance
scatter in `pipeline/03_drift_summary.py`.

- **Source:** derived from Gosselin 2017 supplementary Table S2 TPM
  sheet: mean(log2(TPM+1) across in-vitro-7d microglia samples) −
  mean(log2(TPM+1) across ex-vivo microglia samples).
- **Columns:** `gene`, `log2FC_invitro_vs_exvivo`,
  `exvivo_mean_log2TPM`, `invitro7d_std_mean_log2TPM`.

## Excel-corrupted gene names

Gosselin's supplementary Excel tables were affected by the well-known
Excel gene-name issue (Ziemann et al. 2016). Some HGNC symbols were
auto-converted to datetimes when the tables were saved:

- MARCH1..MARCH11 → 2017-03-01..2017-03-11
- SEPT1..SEPT15   → 2017-09-01..2017-09-15
- DEC1            → 2017-12-01

The reference files here have been corrected to restore the original
HGNC symbols. If you re-extract from the raw Gosselin xlsx yourself,
use `validation/dump_gosselin_counts.py` which applies the same
correction.

## Citation

Gosselin D, Skola D, Coufal NG, et al. **An environment-dependent
transcriptional network specifies human microglia identity.**
*Science* 356, eaal3222 (2017). DOI: 10.1126/science.aal3222
