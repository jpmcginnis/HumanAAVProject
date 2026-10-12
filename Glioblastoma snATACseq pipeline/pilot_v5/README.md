# pilot_v5 — Oct-11 scoring revisions

Scripts that superseded parts of `src/` for the atlas release deposited on
2026-10-11. These are the files that actually generated the per-report
parquets in `analysis_2026-10-11/`. The original `src/` pipeline modules
describe the architecture; these files are the current implementation of
the scoring step.

| file | role |
|---|---|
| `requant.py` | Pseudobulk matrix builder — projects each cohort's h5ad onto the union peak reference and emits patient × cell_type pseudobulks. |
| `requant_greenwald.py` | Serial Greenwald-specific pseudobulker (10.5B-nnz h5ad). |
| `requant_greenwald_parallel.py` | 4-worker parallel version of the above (chunk-level gzip decompress in parallel). |
| `greenwald_finish.py` | 8-worker finisher that merges partial pseudobulks after a parallel scan interrupted. |
| `remap_greenwald_coarse.py` | Collapses Greenwald's 99 fine labels to the 14-label coarse schema used across the other 7 cohorts, retaining the 4 Neftel malignant states as separate buckets. |
| `cnv_patch.py` | hg38 chromosome-arm (centromere-split) coordinates + helpers for arm-conditional Z-score, CNV-class annotation, and per-patient chr7-elevated detection. |
| `score.py` | Per-report scorer: patient-level pseudobulks → `daigle_z_raw` + `daigle_z_arm_adjusted` + CNV columns + patient/cohort replication + off-target metrics. |

## Relationship to `src/`

- `src/scoring.py`, `src/marker_peak_scoring.py`, `src/ranking.py` — ORIGINAL pipeline scoring step. Operates against CATLAS-space pseudobulks with cohort-level replication.
- `pilot_v5/score.py` — CURRENT scoring step. Operates against union-peak pseudobulks with patient-level replication and arm-conditional CNV adjustment.
- Both are kept so the architecture is traceable; the Oct-11 deposit uses `pilot_v5/`.

A future release will refold `pilot_v5/` back into `src/` (likely as `src/scoring_v2.py`) with `src/scoring.py` deprecated or deleted. For now the two live side by side to preserve traceability from the current parquets back to the exact code that produced them.
