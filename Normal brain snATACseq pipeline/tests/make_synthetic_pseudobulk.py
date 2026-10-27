#!/usr/bin/env python3
"""Build a tiny synthetic pseudobulk parquet matching the real atlas schema.
2 donors × 3 cell-type classes × 200 peaks. One peak is engineered to be
MICROGLIA-specific+prevalent so the atlas scoring chain can be verified
end-to-end against a known-correct answer.

Output: tests/synthetic_direct_nnz1k.parquet
"""
import polars as pl
import numpy as np
from pathlib import Path

OUT = Path(__file__).parent / "synthetic_direct_nnz1k.parquet"
STUDY = "synthetic-smoketest"

rng = np.random.default_rng(42)
donors = ["smoke_donor_1", "smoke_donor_2"]
classes = ["Microglia", "Astrocyte", "Oligodendrocyte"]
n_peaks = 200

# Build the (peak × donor × class) accessibility rows.
# Peak 50 is engineered: ~85% open in Microglia across both donors, ~5% in others.
# All other peaks: random fraction_open with no class bias.
rows = []
for donor in donors:
    for cls in classes:
        n_cells_total = int(rng.integers(500, 2000))  # ≥100 passes evaluability gate
        for peak_num in range(n_peaks):
            if peak_num == 50 and cls == "Microglia":
                frac = 0.85 + rng.normal(0, 0.03)
            elif peak_num == 50:
                frac = 0.05 + rng.normal(0, 0.01)
            else:
                frac = rng.beta(1, 10)  # most peaks mostly closed
            frac = float(np.clip(frac, 0.0, 1.0))
            n_access = int(round(frac * n_cells_total))
            if n_access < 3:
                continue  # min cells per peak filter mirrors real pipeline
            rows.append({
                "peak_num": peak_num,
                "major_class": cls,
                "donor_id": donor,
                "n_cells_access_direct": n_access,
                "n_cells_total_direct": n_cells_total,
                "fraction_open_direct": n_access / n_cells_total,
                "study_id": STUDY,
            })

df = pl.DataFrame(rows)
df.write_parquet(OUT)
print(f"wrote {OUT}")
print(f"rows: {df.shape[0]:,}")
print(f"donors: {df['donor_id'].n_unique()}, classes: {df['major_class'].n_unique()}")
# Sanity: expected top Microglia peak
mg = df.filter((pl.col("major_class")=="Microglia") & (pl.col("peak_num")==50))
print(f"engineered microglia peak_50: {mg['fraction_open_direct'].to_list()}")
