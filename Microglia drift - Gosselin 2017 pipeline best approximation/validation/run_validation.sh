#!/usr/bin/env bash
# End-to-end pipeline validation on Gosselin's own human microglia data.
# Reproduces Gosselin's published headline numbers using the same
# pipeline (pipeline/02_limma_voom.R + pipeline/03_drift_summary.py)
# that end users run on their own data.
#
# Requires the Gosselin 2017 supplementary Table S2 xlsx:
#   aal3222_gosselin_tables2.xlsx (downloadable from Science's supp).
#
# Usage:
#   validation/run_validation.sh /path/to/aal3222_gosselin_tables2.xlsx
#
# Produces validation/output/ with:
#   * gosselin_read_counts.csv         (from Table S2)
#   * gosselin_metadata_N5_avg_all_IV.csv
#   * limma_voom_de.csv                (per-gene DE)
#   * limma_voom_summary.csv           (WT DE counts)
#   * drift_headline.txt               (headline vs Gosselin)
#   * drift_headline_bar.png           (comparison plot)

set -euo pipefail

XLSX="${1:-}"
if [[ -z "$XLSX" ]]; then
  echo "Usage: $0 /path/to/aal3222_gosselin_tables2.xlsx" >&2
  exit 1
fi

REPO="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$REPO/validation/output"
mkdir -p "$OUT"

echo "=== [validation 1/3] extract Gosselin counts + metadata ==="
python "$REPO/validation/dump_gosselin_counts.py" \
    --xlsx "$XLSX" --outdir "$OUT"

echo
echo "=== [validation 2/3] run limma-voom on Gosselin data ==="
Rscript "$REPO/pipeline/02_limma_voom.R" \
    --counts   "$OUT/gosselin_read_counts.csv" \
    --metadata "$OUT/gosselin_metadata_N5_avg_all_IV.csv" \
    --outdir   "$OUT"

echo
echo "=== [validation 3/3] compute drift headline vs Gosselin ==="
python "$REPO/pipeline/03_drift_summary.py" \
    --de        "$OUT/limma_voom_de.csv" \
    --signature "$REPO/references/gosselin_conserved_signature_477.csv" \
    --outdir    "$OUT" \
    --label     "Gosselin data\n(this pipeline)"

echo
echo "=== Validation complete. Compare to expected numbers: ==="
cat "$REPO/validation/expected_numbers.md"
