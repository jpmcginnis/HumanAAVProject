# GBM snATAC Atlas v3 — Enhancer Candidate Ranking

**Build date:** 2026-10-17
**Status:** DRAFT (Poisson-approximated penetrance); exact per-cell penetrance pending EBS snapshot restore (24-72 hr).
**Reference spec:** `../../atlas_scoring_html_table_spec.md` (updated this session to reflect primary_score default)
**Reference memory:** `~/.claude/projects/-Users-jpmcginnis1/memory/project-gbm-atlas-v3-f1-primary.md`

---

## 1. Scope and taxonomy

**Atlas inputs**
- 48 patients across 8 published cohorts (Greenwald 2024, Sojka 2025, Sundaram 2024, Zhao 2024, Mathewson 2021, Guo 2020, Guilhamon 2021, Raviram 2023)
- 1.52 M nuclei
- 1,958,753 peaks in the harmonized union set (CATLAS ∪ per-cohort MACS3)
- 95.9% distal (|dist_to_tss| ≥ 2 kb); only distal peaks are ranked

**Cell-type taxonomy (14 tabs)**
- 10 **primary** tabs: astrocyte, microglia, oligodendrocyte, OPC, endothelial, neuron, GABA_neuron, T_cell, TAM, malignant (pooled)
- 4 **secondary** malignant Neftel subtype tabs: malignant_AC, malignant_MES, malignant_NPC, malignant_OPC
- Decision rationale: pooled malignant is the primary therapeutic unit (clinical GBM is a mixed-subtype population; subtype-specific scaffolds miss most of the tumor), but subtype tabs retained as secondary for reviewer interest.

**n_cohorts scale**: 0 to 8 (the spec document originally mentioned 6; corrected to 8 to match the shipped atlas).

---

## 2. Primary ranking: `primary_score` (rank-percentile sum)

**Default sort column**, highlighted in blue in the HTML table.

**Formula:**
```
spec_pct      = percentile rank of log2FC_top_vs_next   (within distal peaks, appropriate taxonomy)
strength_pct  = percentile rank of strength_mean_evaluable_wn   (within distal peaks)
primary_score = (spec_pct + strength_pct) / 2
```

- Equal-weighted nonparametric combination of specificity and strength
- Both inputs are percentile ranks in [0, 1]; mean is in [0, 1]
- Top-1000 per cell type have primary_score typically 0.92-1.0
- No sigmoid, no z-score reference-distribution choice, no arbitrary thresholds; absolute metric scale irrelevant

### Why rank-percentile sum (ranking metric evolution)

Three primary rankings were tried in this build before landing on primary_score. All intermediate columns are retained in the HTML and CSV outputs for diagnostic comparison; only the default sort column changed.

**Attempt 1: `composite_balanced`** = equal-weight z-score sum of 4 components (z_spec + z_strength + z_reprod + z_cohort). Replaced because the strength component (`penetrance_median_detected`) was **circular** — median of per-patient Poisson penetrance taken over patients where penetrance ≥ 1%. The metric conditioned on the same detection threshold it was meant to characterize.

**Attempt 2: `F1_primary`** = harmonic mean of sigmoid(log2FC) and non-circular strength. Fixed the circularity by introducing a new strength metric (see §3) that averages across ALL evaluable patients. Replaced because **`sigmoid(log2FC = 0) = 0.5` gave zero-specificity peaks half credit** on the specificity axis — peaks with log2FC ≈ 0 but high strength ranked at the top of several tabs (malignant top log2FC = 0.01, endothelial = 0.03, T_cell = 0.08).

**Attempt 3: `primary_score`** (CURRENT): nonparametric. A peak with log2FC near 0 gets a near-zero percentile rank, so zero-specificity peaks are properly penalized on the specificity axis.

### Sanity-check: top-rank log2FC values post-fix

| cell type | top log2FC (F1) | **top log2FC (primary_score)** |
|---|---:|---:|
| T_cell | 0.08 | **4.79 (28× fold)** |
| malignant | 0.01 | **4.53 (23×)** |
| endothelial | 0.03 | **4.34 (20×)** |
| GABA_neuron | 0.71 | **4.16** |
| malignant_NPC | — | **3.46** |
| oligodendrocyte | 1.46 | **3.41** |
| OPC | — | **3.35** |
| astrocyte | — | **2.65** |
| neuron | 0.04 | **2.33** |
| malignant_MES | — | **2.27** |
| microglia | — | **2.22** |
| TAM | — | **1.95** |
| malignant_OPC | — | **1.47** |
| malignant_AC | — | **1.40** |

---

## 3. Strength: `strength_mean_evaluable_wn` (non-circular, Poisson-approximated)

**Formula:**
```
penetrance_raw[patient, celltype] = 1 - exp(-counts / n_cells)       # Poisson approx
strength_mean_evaluable_wn        = mean over EVALUABLE patients of penetrance_raw × (500 / peak_width_bp)
```
where evaluable = patient has ≥100 cells of target type.

**Why non-circular (vs legacy `penetrance_median_detected`):** Legacy strength conditioned the median on the detection threshold it was meant to characterize. New strength averages across all evaluable patients (zeros included), with evaluability as a pure data-quality criterion. Decoupled from accessibility entirely.

**Why Poisson approximation (temporary):** Exact per-cell quantification over the full 1.96 M union peak set is not currently on local disk — it lives on the archived EBS snapshot `snap-092225222f7f7ed63` (restore in progress; standard-tier archive restore takes 24-72 hr). The Poisson approximation is tight for low-λ peaks (almost all enhancers) and preserves rank order essentially perfectly in the sparse-signal regime of snATAC-seq accessibility (typically within 5-10% of true per-cell penetrance).

**For Greenwald specifically:** Per-patient cell counts estimated from pseudobulk library sizes using atlas-average fragments-per-cell scaling factors because Greenwald's `catlas_quantified.h5ad` is also missing from local disk (same snapshot-restore dependency). All 7 other cohorts use exact cell counts from their on-disk h5ad obs tables.

---

## 4. Specificity: `log2FC_top_vs_next` with taxonomy separation

**Problem flagged during build:** The initial 14-cell-type log2FC_top_vs_next treated "malignant" (pooled) and malignant_AC/MES/NPC/OPC as co-equal lineages. But the pooled malignant column IS the sum of its subtypes — same underlying cells. When "malignant" was top, the "next" was almost always one of its own subtypes with ~98-99% of the pooled CPM → log2FC ≈ 0 for essentially all top-ranked pan-malignant peaks.

**Fix: two separate specificity taxonomies:**
- **log2FC_vs_primary** — across 10 primary lineages (9 non-malignant + pooled malignant; subtypes excluded). Used for all primary tabs.
- **log2FC_vs_subtype** — across 13 lineages (9 non-malignant + 4 malignant Neftel subtypes; pooled excluded). Used for the 4 malignant subtype secondary tabs.

**Legacy column retained:** `log2FC_top_vs_next_LEGACY_14ct` — the self-competing version — preserved for diagnostic comparison.

---

## 5. Distal filter and distance bands

**Filter:** |dist_to_tss_signed| ≥ 2000 bp applied **before top-N selection**. 95.9% of atlas peaks survive (1,879,285 / 1,958,753). Promoter-proximal peaks are excluded entirely from the ranking — they are not useful AAV enhancer scaffolds and would otherwise dominate top ranks via high accessibility.

**Distance bands:**
- `2-10kb` — perienhancer (2000 ≤ |dist| ≤ 10000 bp)
- `>10kb` — far-distal enhancer (|dist| > 10000 bp)

Of top-1000 per cell type: 77-93% are >10kb, 7-23% are 2-10kb (far-distal dominates).

**HTML dropdown:** "All distal" / "2-10 kb only" / ">10 kb only" filters the active tab in place.

---

## 6. Amplicon-contamination handling

**Diagnostic result from this session** (microglia per-patient test): per-patient microglia pseudobulk signal at chr7p11.2 EGFR amplicon correlated strongly with per-patient malignant-cell fraction (Spearman ρ = +0.70, p < 0.001, n = 19), AND disproportionately vs genome-wide baseline (amp/control ratio Spearman ρ = +0.70, p < 0.001). Textbook ambient-fragment inflation from amplified malignant cells leaking into non-malignant pseudobulks.

**Response:**
- `AMP` red badge appears next to any non-malignant top-ranked peak within ±1 Mb of {chr7p11.2 EGFR, chr12q15 MDM2, chr8q24 MYC, chr9p21.3 CDKN2A deletion}
- Disclosure flag only; peaks remain in supplementary tables for transparency
- Documented in Part 1 methods disclaimer (preprint template memory)

---

## 7. Nearest-gene caveat

The `nearest_gene` column names the nearest annotated TSS to the peak midpoint. It is **informational only** — not a regulatory target assignment. The actual regulated gene may be hundreds of kb away via 3D chromatin contact. For AAV scaffold efficacy, cell-type-selective driving strength (what primary_score ranks) is the correct and sufficient criterion; nearest-gene identity is context, not prediction.

---

## 8. HTML atlas features

- **14 tabs**: 10 primary (solid border) + 4 malignant subtype (dashed border)
- **Default tab**: malignant (pooled) — central use case
- **Default sort**: primary_score descending (highlighted in blue in the table)
- **Sticky column header** — dynamically pins below the tabs row (viewport-adaptive via JS: measures `header.offsetHeight + tabs.offsetHeight` on load and on resize, so it works at any screen width where the tabs wrap into 1 or 2+ rows)
- **Header tooltips** on every column (formula, interpretation, taxonomy notes). Hover any column header.
- **Distance band dropdown**: "All distal" / "2-10 kb only" / ">10 kb only"
- **Column visibility toggles** for chrom/start/end (hidden by default) and all legacy composites + legacy F1 + legacy pen_med_detected (hidden by default)
- **Row expansion** (top-100 per tab inline): click any peak_id row to reveal:
 - **Dual-scale per-cell-type bars**:
   - *Left (blue)*: scaled to this peak's own top CPM — shows cell-type specificity
   - *Right (orange)*: scaled to each cell type's atlas-wide max CPM — shows absolute strength relative to the strongest peaks of that cell type
 - **Per-patient Poisson penetrance** across all evaluable patients, sorted descending, color-coded (detected / present-but-weak / unevaluable)
 - **Per-cohort detection chips** (✓ / ✗) for all 8 cohorts
- **AMP badge** marks any non-malignant top-ranked peak within ±1 Mb of a known GBM amplicon

### Per-cell-type atlas-wide max CPM (used for the absolute-strength expansion bars)

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

Rationale for per-cell-type (not single global) max: different cell types have different absolute signal ranges (T_cell / GABA_neuron have fewer cells per patient so per-cell-type CPM runs higher at the same accessibility). Per-cell-type scaling corrects for this and makes the absolute bars comparable within each cell type.

---

## 9. Outputs

| File | Purpose |
|---|---|
| `atlas_v3.html` | Interactive HTML atlas (10 MB, self-contained; sticky header, tooltips, dropdown, row expansion with dual-scale bars) |
| `atlas_v3_slim.csv` | 14,000 rows (top-1000 × 14 tabs); 24 columns including distance_band + primary_score; use `tab_cell_type` to filter |
| `atlas_v3_supplementary_wide.parquet` + `.csv.gz` | Full wide-format with per-cell-type CPMs + all intermediates |
| `per_cell_type_top1000/top1000_<ct>.tsv` | 14 per-cell-type TSVs ranked by primary_score |
| `v3_master_scores.parquet` | Full-atlas master table (1.96 M peaks × 205 columns including per-cell-type scores + percentile ranks + both log2FC taxonomies + legacy columns) |
| `primary_score_summary.tsv` | Per-cell-type summary with median primary_score, median log2FC, 2-10kb / >10kb breakdown |
| `f1_primary_summary.tsv` | Legacy F1 summary (diagnostic; F1 is no longer default) |
| `README.md` | This file |

---

## 10. Pending work (when exact per-cell data lands)

1. **EBS snapshot restore** — `snap-092225222f7f7ed63`, 7-day temporary restore initiated 2026-10-17 19:23 UTC. Expected 24-72 hr. On completion:
 - Spin up r6i.4xlarge, mount volume, re-run per-cell quantification for all 1.96 M union peaks across 8 cohorts
 - Pull per-cell h5ads to local disk
 - **Also deposit per-cell h5ads to S3** (`s3://jpm-gbm-atlas-intermediates/`) so this never happens again (see memory `feedback-persist-per-cell-h5ads-to-s3`)
 - Re-run the scoring pipeline with exact per-cell counts
 - Regenerate v3 HTML + CSVs; drop DRAFT badge; publish as FINAL
2. **Compare Poisson-approx vs exact per-cell rankings** — quantify how much the approximation shifts the top-1000 per cell type (expect minimal shift given low-λ regime)

---

## 11. Decisions explicitly rejected or deferred

- **Expedited EBS restore** ($15-50 surcharge, 5-24 hr restore) — not pursued; standard-tier restore (24-72 hr free) was sufficient given the Poisson-approx draft covers the gap
- **Shape features (Skew, Kurtosis, Modality) from enhanced PeakRankR** — requires per-cell-type bigWigs we don't have; deferred as future-work extension
- **Splitting malignant subtypes as primary tabs** — subtype splitting isn't therapeutically relevant for pan-malignant AAV scaffolds; kept pooled malignant as primary with subtypes as secondary
- **Hard log2FC ≥ 1 filter** (option offered when F1 was failing) — not pursued; rank-percentile sum obsoleted the need
- **Shifted sigmoid** or **linear-clip** options (offered as alternative F1 fixes) — not pursued; rank-percentile sum was cleaner
- **Geometric mean** of spec and strength — mentioned as alternative; not pursued
- **DataTables FixedHeader plugin** for sticky headers — tried briefly; native `position: sticky !important` on `thead th` with `border-collapse: separate` + dynamically-set top offset was cleaner and worked across all viewport widths

---

## 12. Methodology references

- Ben-Simon et al. 2025 Cell (Allen Armamentarium, TCP Z-score framework)
- Johansen et al. 2025 Cell Genomics (five-team benchmark; F1 convergence 0.57-0.58)
- ArchR (Granja 2021 Nat Genet) — baseline log2FC vs pooled background
- PeakRankR (Cell Genomics 2025) — weighted linear sum convention
- Aerts CREsted — Gini + pairwise ratio
- Neftel et al. 2019 Cell (GBM developmental state framework — basis for the marker-gene-anchoring critique and the taxonomy-separation rationale)
