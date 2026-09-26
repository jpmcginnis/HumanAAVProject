# Validation

This folder reruns the pipeline on Gosselin 2017's own published human
microglia data. It exists so any user can confirm that the pipeline
produces the numbers Gosselin published *before* applying it to their
own samples.

## Files

- **`dump_gosselin_counts.py`** — Extracts the read-count matrix from
  Gosselin 2017 supplementary Table S2 (`aal3222_gosselin_tables2.xlsx`,
  Read Counts sheet) and emits it as a plain CSV. Fixes the well-known
  Excel-corrupted gene-name issue (MARCH1→2017-03-01,
  SEPT1→2017-09-01, DEC1→2017-12-01) that persists in the published
  Table S2. Also emits sample-metadata CSVs for the N=5 paired
  ex-vivo → 7-day in-vitro human microglia design.
- **`run_validation.sh`** — End-to-end script. Takes the path to the
  Gosselin Table S2 xlsx as a single argument, runs `dump_gosselin_counts.py`
  then `pipeline/02_limma_voom.R` then `pipeline/03_drift_summary.py`,
  and prints the head-to-head numbers.
- **`expected_numbers.md`** — What the validation should produce, with
  tolerances. Reproduces Gosselin's primary 33 % drift headline to
  within 1 percentage point; whole-transcriptome DOWN count to within
  7 genes out of 3,702.

## Where to get the Gosselin supplementary data

- Paper: [Gosselin et al. 2017, *Science* eaal3222](https://www.science.org/doi/10.1126/science.aal3222)
- Supplementary Table S2 (Human RNA-seq): downloadable from the
  Science supplementary materials page for that paper. The file we
  need is `aal3222_gosselin_tables2.xlsx`.

## Usage

```bash
# From the repo root, after installing R + limma + edgeR + Python deps:
bash validation/run_validation.sh /path/to/aal3222_gosselin_tables2.xlsx
```

Output lands in `validation/output/`. Compare its `drift_headline.txt`
and `limma_voom_summary.csv` to `expected_numbers.md`.
