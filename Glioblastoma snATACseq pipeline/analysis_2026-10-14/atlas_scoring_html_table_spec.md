# GBM snATAC Atlas — Enhancer Scoring and HTML Table Specification

**Purpose:** This document specifies the enhancer candidate scoring methodology and the HTML results table structure for the harmonized GBM snATAC atlas preprint and companion resource. It is written to be self-contained — another session, another lab member, or a future reviewer should be able to read this document and understand exactly what every column means, how every number was computed, and why the design choices were made.

**Scope:** Atlas of 48 GBM patients across 6 published cohorts, harmonized to a union peak set. Scoring framework is cell-type-agnostic and applies to every cell type in the atlas taxonomy (malignant subclasses and non-malignant brain cell types alike).

---

## 1. Design philosophy

### 1.1 The field has no consensus metric

No single enhancer-nomination metric dominates in the field. Published approaches include ArchR log2FC vs pooled background, Aerts/CREsted Gini + pairwise ratio, PeakRankR weighted linear sum, Ben-Simon per-enhancer Z-score across subclass dimension, Gillis per-bin genome-wide Z-score, and COSG cosine-similarity marker scoring. When Johansen et al. 2025 Cell Genomics benchmarked five of these teams on common Allen Armamentarium screening data, F1 converged at ~0.57–0.58 across the top three methods. **The specificity metric choice is not the field bottleneck. Biology is the ceiling.**

This gives us permission to design a scoring function around the specific failure modes that matter for our use case (AAV-mediated CNS gene delivery in GBM), defend it on first principles, and expose component columns so downstream users can re-weight for their own applications.

### 1.2 No marker-gene anchoring

We deliberately do not use proximity to pre-specified marker genes as a nomination filter. Three reasons:

1. **Empirical.** Ben-Simon et al. 2025 Cell directly tested marker-gene proximity against in vivo enhancer performance across 802 enhancers and found it does not predict labeling specificity.
2. **Biological.** The median mammalian enhancer-promoter distance is ~50–100 kb with a long tail beyond 1 Mb; narrow proximity windows (±50 kb Hunker convention) miss most of the regulatory architecture.
3. **Context-specific.** GBM malignant cells recapitulate neurodevelopmental expression states (Neftel 2019) rather than expressing novel marker genes. Anchoring on "GBM marker genes" like SOX2 or OLIG2 would preferentially nominate enhancers also active in normal NPCs, OPCs, and astrocytes — precisely the off-target cells that limit therapeutic specificity.

We instead evaluate every accessible peak in the atlas on cell-type selectivity, strength, and reproducibility metrics without any pre-filter based on marker-gene proximity. `nearest_gene` and `dist_to_tss_signed` are reported as context only, not used in ranking.

### 1.3 Multi-component scoring with exposed components

Rather than collapse quality into a single scalar, we report four orthogonal quality components and three composite scores derived from them. Reviewers and downstream users see all four components and can derive their own composites. The three default composites represent three different weighting philosophies (balanced, safety-weighted, multiplicative).

---

## 2. The four quality components

Every peak in the atlas is scored on these four orthogonal dimensions. Each is a continuous value, not a thresholded gate.

### 2.1 Specificity: `log2FC_top_vs_next`

**Formula:** `log2(a_top / a_second)`

Where:
- `a_top` = pseudobulk mean accessibility in the cell type with the highest accessibility for this peak (identified in `top_celltype`)
- `a_second` = pseudobulk mean accessibility in the next-highest cell type in the full atlas taxonomy, regardless of whether that second cell type is malignant or non-malignant

**Interpretation:** How much more accessible the enhancer is in its top cell type than in its closest competitor. Measured in log2 units for conventional readability (log2FC > 1 ≡ ratio > 2, log2FC > 2 ≡ ratio > 4).

**Rationale for the choice of denominator:** The second-highest cell type (rather than pooled background, as in ArchR's default) is the correct denominator for AAV therapeutic design because the dominant failure mode is dual-population labeling — if an AAV lands in both the target cell type and one other cell type at similar levels, you get a two-population label regardless of what the "average background" looks like. The pairwise top-vs-next metric is monotonically equivalent to a log2FC with the second-highest cell type as reference and is more conservative than ArchR's convention for the dual-cell-type failure mode. It matches the Aerts CREsted DAR training convention (ratio > 2).

**Scale:** Theoretically unbounded above; in practice ranges roughly 0 to 8 in atlas data. Pseudocounts are added to both numerator and denominator to prevent division by zero (standard offset of 0.01 CPM).

### 2.2 Strength: `strength_mean_evaluable_wn` (PRIMARY, non-circular)

**Formula:**
```
penetrance_raw[patient, celltype] = (cells of celltype in patient with ≥1 Tn5 fragment in peak)
                                  / (total cells of celltype in patient)

strength_mean_evaluable_wn = mean over ALL EVALUABLE patients of penetrance_raw
                            × (500 / peak_width_bp)
```

**Interpretation:** Across all patients with sufficient cells of the target type to be evaluable (≥100 cells), what fraction of those cells, on average, have at least one Tn5 insertion in the peak. Measures how reliably the enhancer is "on" within the target population — WITHOUT conditioning on whether the enhancer is detected in a given patient.

**Why width-normalized:** Peak width varies based on peak-caller choice and biological context; a 2kb peak will accumulate more Tn5 cuts per cell at the same per-bp accessibility as a 200bp peak. Normalizing to a 500bp reference width makes penetrance values comparable across peaks.

**Why "across all evaluable patients" rather than "across detected patients" (non-circularity rationale):** An earlier version of this spec defined strength as the median over "detected" patients only, where "detected" was itself defined by a per-patient penetrance threshold. This is circular: the strength metric was being conditioned on a filter derived from the same quantity it was supposed to summarize. The current definition takes the mean across ALL patients with enough cells to be evaluated (the evaluability gate at ≥100 cells is a pure data-quality floor; it does not reference the peak's accessibility). This yields a strength value that reflects both how reliably the peak is accessible when biologically present AND how commonly it is biologically present at all — a cleaner, non-circular measure that supports the F1-primary ranking (section 2.5) without hidden dependencies.

**Legacy column retained for comparison:** `penetrance_median_detected_LEGACY` reports the earlier "median over detected patients" value for diagnostic comparison; it is NOT used in the primary ranking.

**Scale:** 0 to 1 (fraction). After width normalization, values can exceed 1 for peaks much narrower than 500bp if accessibility is very high; this is numerically acceptable and preserves the ranking.

### 2.3 Reproducibility within atlas: `patient_detection`

**Formula:**
```
evaluable[patient, celltype] = (patient has ≥100 cells of celltype)

detected[patient, celltype] = evaluable AND (≥1% of patient's cells of that type
                              have ≥1 Tn5 fragment overlapping the peak)

patient_detection_rate = n_patients_detected / n_patients_evaluable
```

**Display format in HTML table:** `73% (30)` — percentage first, with the evaluable denominator in parentheses.

**Sort key:** The numeric `patient_detection_rate` (0 to 1).

**Tooltip on hover:** "Detected in {n_detected} of {n_evaluable} evaluable patients. Evaluable = ≥100 cells of target type. Detected = ≥1% of those cells have the peak accessible."

**Interpretation:** How reproducibly this enhancer is detectable across patients in the atlas.

**Why the two gates:**
- **Evaluability gate (≥100 cells of target type):** A pure data-quality floor. A patient with only 30 cells of a given type cannot support a reliable per-patient detection call either way; they should drop out of the denominator entirely rather than count as a "no." This matters for rare cell types where some patients will have poor representation.
- **Detection gate (≥1% of evaluable cells have the peak):** A permissive threshold that says "there is some evidence this peak is biologically present in this patient." Permissive-for-inclusion is appropriate for enhancer discovery, where you want to catch patients where the peak is real but at low penetrance within the cell population. The 1% floor combined with the 100-cell evaluability floor means detection requires at minimum 1 cell with a fragment (1% of 100 = 1).

**Why both numerator and denominator are shown:** A peak at `100% (3)` reads very differently from `100% (46)`. The denominator calibrates reader trust in the rate without any additional explanation.

**Why not a hard threshold filter:** These gates define what counts as a "yes" or "no" per patient, but the resulting rate is a continuous rank variable, not a filter. Peaks at low detection rates aren't excluded; they rank lower.

### 2.4 Cross-cohort replication: `n_cohorts`

**Formula:** Count of the 6 pooled published cohorts in which the peak clears a minimum detection threshold (defined per the atlas construction methodology, typically ≥5% of cells in that cohort with ≥1 fragment at the peak).

**Interpretation:** Independent replication of the peak's existence across datasets collected by different groups with different protocols and sample handling. Controls for batch effects, technical artifacts, and single-cohort anomalies in a way within-cohort patient replication cannot.

**Scale:** Integer, 1 to 6.

**Why separate from `patient_detection`:** Within-cohort patient replication and cross-cohort study replication are distinct forms of evidence. A peak detected in 30/30 patients of a single cohort (patient_detection = 100%, n_cohorts = 1) is a less robust finding than one detected in 15/45 patients across 4 cohorts (patient_detection = 33%, n_cohorts = 4). Both axes matter and should contribute independently to the composite.

---

### 2.5 Primary ranking: `primary_score` = rank-percentile sum (default sort column)

**Formula:**
```
# Within the distal peak set (|dist_to_tss_signed| ≥ 2000 bp):
spec_pct     = percentile rank of log2FC_top_vs_next in appropriate taxonomy
strength_pct = percentile rank of strength_mean_evaluable_wn (per cell type)
primary_score = (spec_pct + strength_pct) / 2
```

**Interpretation:** Equal-weight nonparametric combination of specificity and strength, each converted to a 0-1 percentile rank across all distal atlas peaks. A primary_score of 0.95 means the peak is in the 95th percentile averaged across specificity and strength. Nonparametric: no sigmoid, no z-score reference distribution choice, no arbitrary thresholds; absolute metric scale irrelevant.

**Why rank-percentile sum is the primary ranking (replaces earlier `f1_primary` default):**
- The earlier F1 formulation used `spec_sigmoid = sigmoid(log2FC_top_vs_next)`, which maps `log2FC = 0 → 0.5` — giving zero-specificity peaks "half credit" on the specificity axis. Combined with high-strength peaks, this allowed broadly-accessible-but-unspecific peaks to rank highly, which is the exact wrong behavior for AAV scaffold selection.
- Percentile ranks don't have the sigmoid-centering issue: a peak with log2FC near 0 gets a near-zero percentile rank, so primary_score penalizes zero-specificity peaks heavily on the specificity axis.
- Equal-weight sum keeps the "both dimensions matter equally" design intent without biasing by variance or distributional assumptions.

**Taxonomy separation (important):**
- **Primary tabs** (astrocyte, microglia, oligodendrocyte, OPC, endothelial, neuron, GABA_neuron, T_cell, TAM, malignant pooled): `log2FC_vs_primary` computed across the 10 primary cell types (pooled malignant counts as one lineage; subtypes excluded).
- **Secondary tabs** (malignant_AC/MES/NPC/OPC): `log2FC_vs_subtype` computed across 13 cell types (9 non-malignant + 4 malignant Neftel subtypes, pooled malignant excluded).
- Each taxonomy avoids pooled-vs-subtype self-competition: an earlier version computed log2FC across all 14 cell types including both pooled and subtypes, so malignant top-ranked peaks had log2FC ≈ 0 (the "second cell type" was always one of its own subtypes with ~98-99% of the pooled CPM).

**Distal filter (applied BEFORE top-N selection):**
- Only peaks with |dist_to_tss_signed| ≥ 2000 bp are ranked and shown in the atlas (95.9% of atlas peaks survive).
- Rationale: promoter-proximal peaks are not useful AAV enhancer scaffolds and would otherwise dominate the top ranks via high accessibility.
- Distance band column: `2-10kb` (2000 ≤ |dist| ≤ 10000 bp) or `>10kb` (|dist| > 10000 bp). HTML has a dropdown to further restrict to one band.

**Scale:** 0 to 1. Both inputs are percentile ranks in [0, 1], so the mean is in [0, 1]. In practice the top-1000 per cell type have primary_score ≥ 0.92 (mostly 0.97-1.0).

**Legacy columns retained:** `f1_primary_LEGACY`, `composite_balanced`, `composite_safety_first`, `composite_multiplicative`, and `penetrance_median_detected_LEGACY` are preserved as toggleable columns for diagnostic comparison. None are the default ranking.

---

## 3. The three composite scores (LEGACY — retained as alternative rankings)

The composites below predate the F1_primary default (section 2.5) and are retained as optional alternative rankings. They are not the default sort column in v3. All composites use **z-scores computed against the full-atlas distribution of each component** (not against the top-ranked subset). This preserves the real variability of each metric and ensures that z-scores retain their meaning as "how exceptional is this peak relative to the atlas as a whole."

### 3.1 `composite_balanced` — Default sort column

**Formula:**
```
z_spec     = (log2FC_top_vs_next     − mean_atlas) / sd_atlas
z_strength = (penetrance_median_detected − mean_atlas) / sd_atlas
z_reprod   = (patient_detection_rate − mean_atlas) / sd_atlas
z_cohort   = (n_cohorts              − mean_atlas) / sd_atlas

composite_balanced = z_spec + z_strength + z_reprod + z_cohort
```

**Interpretation:** Equal-weight sum of z-scores across the four quality dimensions. No tuning assumptions; each dimension contributes equal variance.

**When to use:** As the default general-purpose ranking. The honest default — no assumptions about which failure mode matters most.

### 3.2 `composite_safety_first` — Therapeutic-use-case default

**Formula:**
```
composite_safety_first = 2 × z_spec + z_strength + z_reprod + z_cohort
```

**Interpretation:** Doubles the specificity contribution. Reflects the therapeutic reality that off-target expression in non-target cells is typically a worse failure mode than suboptimal expression in on-target cells, particularly for cytotoxic payloads, gene editors, or other cargo with real clinical risk from off-target activity.

**When to use:** When ranking for therapeutic application contexts.

### 3.3 `composite_multiplicative` — All-four-dimensions-matter

**Formula:**
```
composite_multiplicative = log2FC_top_vs_next
                         × penetrance_median_detected
                         × patient_detection_rate
                         × sqrt(n_cohorts)
```

**Interpretation:** Multiplicative rather than additive, so a peak that scores poorly on any single dimension cannot be rescued by excellence on others. Reflects "we need all four to be meaningful, not just three." The `sqrt(n_cohorts)` dampens the cohort term so a peak in 4 cohorts doesn't dominate purely on replication count.

**When to use:** When making binary clone/don't-clone decisions where any single-dimension failure is disqualifying. Less appropriate as a pure ranking metric because it can crush peaks under-replicated across cohorts even if genuinely strong.

**Why not the default:** Without a tier-threshold system doing gatekeeping work upfront, the multiplicative composite buries peaks exclusive to early cohorts even when they're strong on the other three dimensions. The additive `composite_balanced` handles this more gracefully for a pure ranking use case.

---

## 4. HTML table specification

### 4.1 Structure: top 1000 distal peaks per cell type

The HTML table includes **the top 1000 DISTAL peaks per cell type** (|dist_to_tss| ≥ 2 kb), **ranked by `primary_score` descending within each cell type.** For a taxonomy with N cell types, this yields N tabs, each showing 1000 peaks selected specifically for that cell type's strongest candidates after the distal filter and after taxonomy-correct log2FC computation.

**Rationale:** A single top-1000-across-the-atlas table would be dominated by whichever cell type has the strongest absolute scores (likely abundant malignant populations), hiding candidates for rarer or less-accessible cell types. Per-cell-type ranking ensures users looking for enhancers specific to any cell type in the atlas can find their best candidates.

**Implementation:** For each cell type C in the taxonomy:
1. Filter the full atlas peak set to peaks where `top_celltype == C`
2. Compute all four quality components and all three composites for these peaks (using full-atlas distribution z-scores for the composites)
3. Sort by `composite_balanced` descending
4. Retain top 1000
5. Emit as a tab in the HTML table (or a filterable table where `top_celltype` can be selected)

### 4.2 Column list

Left-to-right order, with hidden-by-default columns noted:

| # | Column | Format / Definition |
|---|---|---|
| 1 | `peak_id` | `chr{N}:{start}-{end}`, GRCh38 |
| 2 | `chrom` | Hidden by default, shown on expand |
| 3 | `start` | Hidden by default |
| 4 | `end` | Hidden by default |
| 5 | `peak_width_bp` | `end − start` |
| 6 | `nearest_gene` | HGNC symbol of nearest TSS to peak midpoint. Context only, not used in ranking |
| 7 | `dist_to_tss_signed` | Signed distance (bp); negative = upstream on gene's sense strand |
| 8 | `top_celltype` | Cell type with highest pseudobulk accessibility across the full atlas taxonomy |
| 9 | `log2FC_top_vs_next` | Specificity. `log2(a_top / a_second)` where `a_second` is the next-highest cell type regardless of malignant/non-malignant status |
| 10 | `penetrance_median_detected` | Strength. Width-normalized median within-patient penetrance across detected patients |
| 11 | `patient_detection` | Displayed as `73% (30)`. Sorts by percentage (stored as 0–1 numeric). Tooltip explains gates |
| 12 | `n_cohorts` | Cross-cohort replication count (1 to 6) |
| 13 | `distance_band` | '2-10kb' or '>10kb' (|dist| < 2 kb filtered out upstream). Dropdown filter in HTML |
| 14 | `spec_pct` | Percentile rank of log2FC_top_vs_next within distal peaks (0-1). Nonparametric 0-1 specificity |
| 15 | `strength_mean_evaluable_wn` | Non-circular strength: mean per-patient Poisson penetrance across all evaluable patients × (500/width) |
| 16 | `strength_pct` | Percentile rank of strength_mean_evaluable_wn within distal peaks (0-1) |
| 17 | `primary_score` | **Default sort column.** `(spec_pct + strength_pct) / 2`. Equal-weighted rank-percentile sum |
| 18 | `f1_primary_LEGACY` | LEGACY. Harmonic mean of sigmoid(log2FC) and strength. Replaced by primary_score |
| 19 | `composite_balanced` | LEGACY. Equal-weight z-score sum of 4 components |
| 20 | `composite_safety_first` | LEGACY. 2× specificity-weighted composite |
| 21 | `composite_multiplicative` | LEGACY. log2FC × penetrance × detection_rate × sqrt(n_cohorts) |
| 22 | `penetrance_median_detected_LEGACY` | LEGACY strength (circular — kept for diagnostic comparison only) |

### 4.3 Table-level behavior

- **Default sort:** `primary_score` descending within each `top_celltype` (f1_primary_LEGACY and composite_balanced retained as toggleable alternatives)
- **All numeric columns sortable** by clicking column header
- **Row expansion** reveals:
  - Per-cell-type pseudobulk accessibility (full accessibility vector that drove the log2FC)
  - Per-patient penetrance values across all evaluable patients
  - Per-cohort replication status (which of the 6 cohorts the peak was detected in)
- **Filter widget on `top_celltype`** (required if implemented as one table rather than tabs)
- **CSV export** splits `patient_detection` back into separate columns:
  - `n_patients_detected`
  - `n_patients_evaluable`
  - `patient_detection_rate` (0 to 1 numeric)

### 4.4 Rendering stack recommendations

DataTables.js or Tabulator both support sortable numeric columns with custom display formatting, header tooltips, row expansion, and CSV export out of the box. Pick based on team preference.

---

## 5. CSV export structure (for computational downstream use)

The downloadable CSV should expose maximum detail for programmatic re-ranking by users with different priorities. Columns:

```
peak_id, chrom, start, end, peak_width_bp,
nearest_gene, dist_to_tss_signed,
top_celltype,
log2FC_top_vs_next,
penetrance_median_detected,
n_patients_detected, n_patients_evaluable, patient_detection_rate,
n_cohorts,
composite_balanced, composite_safety_first, composite_multiplicative
```

Additionally, a supplementary wide-format CSV should expose:
- Per-cell-type pseudobulk accessibility for every cell type in the taxonomy
- Per-patient penetrance values for every patient × cell-type combination
- Per-cohort detection flags

These supplementary files let downstream users recompute any component with modified thresholds or definitions.

---

## 6. Preprint framing

### 6.1 How to describe the scoring in Methods

> "No consensus metric exists for cell-type-specific enhancer nomination from snATAC-seq data. Published approaches range from log2 fold-change against pooled background (ArchR convention, Johansen et al. 2025 benchmark baseline), to Gini index with pairwise ratio filters (Aerts CREsted), to per-enhancer Z-scoring across cell-type dimension (Ben-Simon et al. 2025 Allen Armamentarium flagship), to weighted linear sums of specificity, magnitude and sensitivity components (PeakRankR). When benchmarked head-to-head these methods converge at F1 ≈ 0.57–0.58 (Johansen et al. 2025), indicating no single formulation dominates. We designed a composite nomination score around four orthogonal quality dimensions: specificity (log2 fold-change of top cell type vs. second-highest), strength (width-normalized median within-patient penetrance across detected patients), within-atlas reproducibility (fraction of evaluable patients where the peak is detected, with evaluability requiring ≥100 cells of the target type and detection requiring ≥1% of those cells to have ≥1 Tn5 fragment in the peak), and cross-cohort replication (count of independent published datasets in which the peak clears the atlas detection threshold). Three composite rankings are reported: a balanced equal-weight z-score sum, a safety-weighted variant doubling the specificity contribution for therapeutic applications, and a multiplicative composite enforcing all-four-dimensions-matter structure. All four components are released as separate columns so downstream users can derive application-specific rankings."

### 6.2 How to describe the design choice against marker-gene anchoring

Frame as a biology-driven design choice specific to the system, not as a methodological critique of other groups:

> "Prior enhancer-AAV nomination efforts in normal brain have relied on proximity to cell-type marker genes to prioritize candidate elements (Graybuck et al. 2021; Mich et al. 2021; Hunker et al. 2025). This approach is efficient when the goal is to label cell populations defined by canonical lineage markers but is structurally mismatched to malignancy. GBM malignant cells recapitulate neurodevelopmental expression states without novel marker genes (Neftel et al. 2019); the regulatory elements that distinguish tumor from normal are precisely those active within a shared marker landscape. We adopted an unbiased pan-genome nomination strategy that evaluates every accessible peak on cell-type selectivity metrics without filtering candidates by proximity to pre-specified marker genes. This design choice is also supported by the empirical observation that marker-gene proximity does not predict in vivo enhancer performance (Ben-Simon et al. 2025)."

---

## 7. Field references

Key papers to cite when defending methodological choices:

- **Ben-Simon et al. 2025 Cell** [DOI 10.1016/j.cell.2025.05.002] — Allen Armamentarium flagship; 802 enhancers screened; empirically disconfirms marker-gene proximity and cross-species conservation as predictors; TCP Z-score framework.
- **Johansen et al. 2025 Cell Genomics** [DOI 10.1016/j.xgen.2025.100879] — Five-team benchmark of nomination methods; F1 convergence at 0.57–0.58.
- **Graybuck et al. 2021 Neuron** [DOI 10.1016/j.neuron.2021.03.011] — Early Allen mscRE enhancer selection; Wilcoxon marker-gene DE + manual pileup + phastCons.
- **Mich et al. 2021 Cell Reports** [DOI 10.1016/j.celrep.2021.108754] — Human neocortex enhancer catalog; Homer `-region` + subclass-specific flag + manual pileup.
- **Hunker et al. 2025 Neuron** [DOI 10.1016/j.neuron.2025.04.035] — Striatal enhancer toolbox; marker-gene proximity in BOCA + CATlas bulk ATAC.
- **Vormstein-Schneider et al. 2025 Neuron** [DOI 10.1016/j.neuron.2025.05.002] — Fishell interneuron enhancer toolbox; COSG cosine-similarity ranking + manual bigwig inspection.
- **Neftel et al. 2019 Cell** — GBM developmental state framework (NPC/OPC/AC/MES-like), biological basis for the marker-gene-anchoring critique in the GBM context.

---

## 8. Known open items when resuming

- Tier threshold calibration is no longer needed (top-1000-per-cell-type ranking replaces tier gates)
- Z-score reference distribution choice confirmed: full-atlas peak distribution, not top-ranked subset
- Width normalization reference confirmed: 500bp
- Composite default confirmed: `composite_balanced`
- Table scope confirmed: top 1000 peaks per cell type, not top 1000 overall
- Still to build: pandas pipeline producing the input dataframe, HTML rendering scaffold (DataTables.js or Tabulator), supplementary wide-format CSVs for per-cell-type and per-patient detail
