# GBM snATAC Atlas v4 — Enhancer Candidate Ranking

**Build date:** 2026-10-17
**Change from v3:** Removed the Poisson adjustment from the strength metric. Everything else carried over.
**Status:** Draft (same per-cell penetrance situation as v3 — but this version doesn't need the exact per-cell swap-in because the plain fragments-per-cell metric is already a direct measurement, not an approximation).
**Reference spec:** `../../atlas_scoring_html_table_spec.md`
**v3 artifacts preserved alongside:** `../atlas_v3/` for comparison.

---

## 1. Scope and taxonomy (unchanged from v3)

- 48 patients × 8 published cohorts × 1.52 M nuclei
- 1,958,753 peaks in the harmonized union set (CATLAS ∪ per-cohort MACS3)
- 95.9% distal (|dist_to_tss| ≥ 2 kb); only distal peaks ranked
- 14 ranking groups: 10 primary cell types + 4 malignant Neftel subtypes

Primary tabs: astrocyte, microglia, oligodendrocyte, OPC, endothelial, neuron, GABA_neuron, T_cell, TAM, malignant (pooled). Secondary tabs: malignant_AC/MES/NPC/OPC.

---

## 2. Primary ranking: `primary_score` (rank-percentile sum)

**Formula (unchanged from v3):**
```
spec_pct      = percentile rank of log2FC_top_vs_next   (within distal peaks, appropriate taxonomy)
strength_pct  = percentile rank of strength_frags_per_cell_wn   (within distal peaks)
primary_score = (spec_pct + strength_pct) / 2
```

Equal-weighted nonparametric combination of specificity and strength. Both inputs are percentile ranks in [0, 1]; mean is in [0, 1]. No model assumptions.

**Default sort column** in the HTML, highlighted in blue.

---

## 3. Strength: `strength_frags_per_cell_wn` (v4 — plain, no Poisson)

**Formula:**
```
for each evaluable patient p:
  frags_per_cell_p = (fragments in peak for patient p's cells of this type) / (cells of this type in patient p)
strength_frags_per_cell_wn = (mean over all evaluable patients of frags_per_cell_p) × (500 / peak_width_bp)
```

**Plain-language interpretation:** "Average number of Tn5 fragments per cell at this peak, across patients with ≥100 cells of the target type, normalized to a 500 bp reference width."

- **Evaluable** = patient has ≥100 cells of target cell type (data-quality floor only; does not reference accessibility)
- **Mean (not median)** across patients so high-ceiling peaks pull up — `patient_detection_rate` is a separate column that measures replication across patients
- **Width-normalized** (× 500 / peak_width_bp) so a 2 kb peak and a 200 bp peak are comparable (longer peaks capture more Tn5 insertions proportionally; one-sentence geometric correction, not a model assumption)
- **Range** typically 0 to ~10, median around 0.001-0.004 across cell types

### Why v4 drops the Poisson transform

v3 used `1 - exp(-fragments/cells)` to convert fragments per cell into an approximate per-cell penetrance (fraction of cells with ≥1 fragment). This assumed fragments distribute Poisson across cells, which:

- Is biologically wrong on diploid DNA (max 2 fragments per cell per peak — not Poisson's unbounded tail)
- Accidentally works at low signal (`1 - exp(-λ) ≈ λ` for small λ)
- Systematically underestimates at high signal
- Requires a methods paragraph explaining the approximation

For rank-based scoring, the Poisson transform is also mathematically unnecessary — percentile ranks are invariant to monotone transforms, so dropping it gives essentially the same primary_score ranks within each cell type (tiny shifts because `mean(f(x))` ≠ `f(mean(x))` across patients when f is nonlinear, but preserved qualitatively).

v4 reports the strength as a direct measurement — "average fragments per cell, width-normalized" — with zero model assumption to justify.

### Also: the exact per-cell EBS swap-in becomes much less critical

In v3 the plan was to swap in exact per-cell penetrance values when the archived EBS snapshot restore completes. With the v4 plain-strength definition, there's nothing to swap — fragments and cell counts are already in hand. The EBS restore is still useful for cleaner per-patient penetrance display and for future work, but no longer gates a v5.

---

## 4. Specificity: `log2FC_top_vs_next` with taxonomy separation (unchanged from v3)

Two separate specificity taxonomies to avoid pooled-vs-subtype self-competition:

- **log2FC_vs_primary** — across 10 primary lineages (9 non-malignant + pooled malignant; subtypes excluded). Used for all primary tabs.
- **log2FC_vs_subtype** — across 13 lineages (9 non-malignant + 4 malignant Neftel subtypes; pooled excluded). Used for the 4 malignant subtype secondary tabs.

Rationale: pooled malignant IS the sum of its 4 subtypes. If all 14 cell types were competing, pan-malignant peaks would have log2FC ≈ 0 because the "second" cell type would be one of the pooled's own subtypes.

---

## 5. Distal filter and distance bands (unchanged from v3)

- Filter: |dist_to_tss_signed| ≥ 2000 bp applied before top-N selection (95.9% survive)
- `distance_band`: `2-10kb` (perienhancer) or `>10kb` (far-distal)
- HTML dropdown: "All distal" / "2-10 kb only" / ">10 kb only"

Of top-1000 per cell type: 77-93% are >10kb, 7-23% are 2-10kb.

---

## 6. Replication: `patient_detection_rate` and `n_cohorts` (unchanged from v3)

- **patient_detection_rate** = fraction of evaluable patients where per-patient width-normalized frags-per-cell exceeds 0.01 (was Poisson penetrance ≥ 0.01 in v3; equivalent threshold on the plain metric)
- **n_cohorts** = count of cohorts (0-8) where cohort-level width-normalized frags-per-cell exceeds 0.05
- Both reported alongside strength to let users see "how strong when active" AND "how often active" independently

---

## 7. Amplicon-contamination caveat (permanent, unchanged from v3)

At known CNV-amplified loci (chr7p11.2 EGFR, chr12q15 MDM2, chr8q24 MYC, chr9p21.3 CDKN2A deletion), non-malignant pseudobulks show elevated apparent accessibility from ambient-fragment contamination. `AMP` red badge marks any non-malignant top-ranked peak within ±1 Mb of these loci. Diagnostic that confirmed the mechanism (v3 session): microglia chr7p signal correlated with per-patient malignant-cell fraction at Spearman ρ = +0.70, p < 0.001.

---

## 8. Nearest-gene caveat (unchanged from v3)

The `nearest_gene` column names the nearest annotated TSS to the peak midpoint. **Informational only** — not a regulatory target assignment.

---

## 9. HTML atlas features (unchanged from v3)

- 14 tabs: 10 primary (solid border) + 4 secondary malignant subtype (dashed border)
- Default tab: malignant (pooled)
- Default sort: primary_score descending
- **Sticky column header** — viewport-adaptive via JS (measures header + tabs height on load and resize)
- **Header tooltips** on every column
- **Distance band dropdown**
- **Column visibility toggles** for chrom/start/end (hidden by default)
- **Row expansion** (top-100 per tab inline): click any peak_id row to reveal:
 - **Dual-scale per-cell-type bars** — blue bars scaled to this peak's own top CPM (cell-type specificity); orange bars scaled to each cell type's atlas-wide max CPM (absolute strength)
 - **Per-patient fragments-per-cell** across all evaluable patients, sorted descending, color-coded
 - **Per-cohort detection chips** (✓ / ✗) for all 8 cohorts
- **AMP badge** marks amplicon-adjacent non-malignant top-ranked peaks

### Per-cell-type atlas-wide max CPM (used for the orange absolute-strength bars)

| cell type | atlas-wide max CPM |
|---|---:|
| GABA_neuron | 389.1 |
| endothelial | 306.9 |
| microglia | 250.9 |
| astrocyte | 242.7 |
| malignant_MES | 230.2 |
| malignant_AC | 205.5 |
| malignant_OPC | 188.5 |
| T_cell | 186.1 |
| malignant_NPC | 183.9 |
| malignant (pooled) | 182.9 |
| oligodendrocyte | 167.9 |
| TAM | 147.8 |
| neuron | 134.6 |
| OPC | 123.3 |

---

## 10. Top-rank sanity-check (v4 vs v3)

Removing the Poisson transform shifts ranks very slightly (because `mean(1 - exp(-λ_p))` ≠ `1 - exp(-mean(λ_p))` across patients). Top log2FC values per tab:

| cell type | v3 top log2FC | v4 top log2FC |
|---|---:|---:|
| T_cell | 4.79 | 5.98 |
| malignant | 4.53 | 4.53 |
| endothelial | 4.34 | 4.16 |
| GABA_neuron | 4.16 | 4.16 |
| oligodendrocyte | 3.41 | 3.41 |
| OPC | 3.35 | 3.35 |
| malignant_NPC | 3.46 | 3.46 |
| TAM | 1.95 | 3.80 |
| astrocyte | 2.65 | 2.65 |
| microglia | 2.22 | 2.22 |
| neuron | 2.33 | 2.33 |
| malignant_MES | 2.27 | 2.27 |
| malignant_OPC | 1.47 | 1.47 |
| malignant_AC | 1.40 | 1.40 |

Most cell types: identical top peak. T_cell and TAM top peaks shifted — honest signal from removing the per-patient Poisson compression.

---

## 11. Outputs

| File | Purpose |
|---|---|
| `atlas_v4.html` | Interactive HTML atlas (7.8 MB, self-contained; sticky header, tooltips, dropdown, row expansion with dual-scale bars) |
| `atlas_v4_slim.csv` | 14,000 rows (top-1000 × 14 tabs); 19 columns; use `tab_cell_type` to filter |
| `atlas_v4_supplementary_wide.parquet` + `.csv.gz` | Full wide-format with per-cell-type CPMs and intermediates |
| `per_cell_type_top1000/top1000_<ct>.tsv` | 14 per-cell-type TSVs ranked by primary_score |
| `v4_master_scores.parquet` | Full-atlas master table (1.96 M peaks × many columns including both log2FC taxonomies, per-cell-type strength, percentile ranks, primary_score) |
| `primary_score_summary.tsv` | Per-cell-type summary with median primary_score, median log2FC, 2-10kb / >10kb breakdown |
| `README.md` | This file |

v3 artifacts preserved in `../atlas_v3/` for comparison. Legacy composites (composite_balanced, F1_primary, etc.) are NOT carried into v4 — they were previous-generation defaults that have no reason to persist alongside the current primary_score. If you need them for diagnostic comparison, pull from `../atlas_v3/`.

---

## 12. Primary ranking evolution (full history)

Preserved for context — four primary rankings were tried in the v3/v4 arc:

1. **composite_balanced** (v3 initial — dropped) — equal-weight z-score sum of 4 components. Replaced because strength component was circular.
2. **F1_primary** (v3 interim — dropped) — harmonic mean of sigmoid(log2FC) and non-circular strength. Replaced because sigmoid(log2FC=0)=0.5 gave zero-specificity peaks half credit.
3. **primary_score with Poisson-based strength** (v3 final) — rank-percentile sum. Replaced in v4 because the Poisson transform was an unnecessary model assumption.
4. **primary_score with plain fragments-per-cell strength** (v4 CURRENT) — rank-percentile sum of specificity and plain strength. No model assumptions anywhere in the ranking pipeline.

---

## 13. Methodology references

- Ben-Simon et al. 2025 Cell (Allen Armamentarium, TCP Z-score framework)
- Johansen et al. 2025 Cell Genomics (five-team benchmark; F1 convergence 0.57-0.58)
- ArchR (Granja 2021 Nat Genet) — baseline log2FC vs pooled background
- PeakRankR (Cell Genomics 2025) — weighted linear sum convention
- Aerts CREsted — Gini + pairwise ratio
- Neftel et al. 2019 Cell (GBM developmental state framework — basis for the marker-gene-anchoring critique and the taxonomy-separation rationale)
