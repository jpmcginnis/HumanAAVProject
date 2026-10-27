#!/bin/bash
# Atlas-build smoke test.
# Builds a tiny synthetic 2-donor × 3-class pseudobulk, runs the scoring chain, and
# verifies that an engineered Microglia-specific peak lands as the top Microglia hit.
#
# Takes ~15 sec. If it fails, something in nb_specificity_recompute.py,
# nb_counts_recompute.py, fix_strength_to_gbm.py, or the overall peak-indexing is broken.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
TMP="$HERE/_smoke_tmp"
mkdir -p "$TMP"
cd "$TMP"
rm -f *.parquet

echo "[1/5] synth pseudobulk"
python3 "$HERE/make_synthetic_pseudobulk.py"
mv "$HERE/synthetic_direct_nnz1k.parquet" "$TMP/synth_direct_nnz1k.parquet"

echo "[2/5] run scoring chain (lightweight — direct per-peak aggregation, no external parquets)"
python3 "$HERE/run_scoring_on_synth.py" "$TMP/synth_direct_nnz1k.parquet" "$TMP/synth_scored.parquet"

echo "[3/5] verify top Microglia peak is the engineered one"
python3 - "$TMP/synth_scored.parquet" <<'PY'
import polars as pl, sys
d = pl.read_parquet(sys.argv[1])
mg = d.filter(pl.col("major_class")=="Microglia").sort("strength_wn_gbm", descending=True)
top = mg.head(1)
print(top.select(["peak_num","mean_frac","log2FC_fresh","strength_wn_gbm","n_donors","n_evaluable_class","patient_detection_rate_gbm"]))
peak = int(top["peak_num"][0])
mean_frac = float(top["mean_frac"][0])
log2fc = float(top["log2FC_fresh"][0])
nd = int(top["n_donors"][0])
nev = int(top["n_evaluable_class"][0])
assert peak == 50, f"expected engineered peak 50, got {peak}"
assert mean_frac > 0.6, f"expected mean_frac>0.6, got {mean_frac}"
assert log2fc > 3.0, f"expected log2FC>3, got {log2fc}"
assert nd == 2, f"expected 2 donors detected, got {nd}"
assert nev == 2, f"expected 2 evaluable, got {nev}"
print("OK — engineered Microglia peak_50 is top-ranked with correct stats")
PY

echo "[4/5] check patient_detection_rate math"
python3 - "$TMP/synth_scored.parquet" <<'PY'
import polars as pl, sys
d = pl.read_parquet(sys.argv[1])
rate = d.filter((pl.col("major_class")=="Microglia") & (pl.col("peak_num")==50))["patient_detection_rate_gbm"][0]
assert abs(rate - 1.0) < 1e-6, f"expected 100% detection, got {rate}"
print(f"OK — patient_detection_rate = {rate:.3f}")
PY

echo "[5/5] cleanup"
rm -rf "$TMP"
echo "SMOKE TEST PASSED"
