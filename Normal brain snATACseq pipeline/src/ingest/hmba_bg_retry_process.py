#!/usr/bin/env python3
"""HMBA-BG retry — process h5ads that FAILed download during the main continuation run.
After re-downloading them via a serial `aws s3 cp` loop into a separate directory, this
script bypasses `process_donor`'s own download step and feeds local paths directly to
`process_h5ad`, writing to a SEPARATE parquet (`hmba-bg-retry_direct_nnz1k.parquet`) so it
can't race with the still-running main continuation on the same file. The two parquets are
merged later by the scoring step (which globs `*_direct_nnz1k.parquet`).

Deletes each h5ad as soon as it's processed to recover disk.

Reconstructed from conversation transcript 2026-10-27.
"""
from __future__ import annotations
import gc, sys, time
from collections import defaultdict
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
from hmba_bg_process import (
    load_union_peaks, process_h5ad,
    NNZ_FLOOR, MIN_CELLS_PER_PEAK, STUDY_ID, OUT_DIR, AUDIT,
)

RETRY_DIR = Path("/tmp/hmba_retry")  # local dir holding the retry-downloaded h5ads
OUT_PQ = OUT_DIR / "hmba-bg-retry_direct_nnz1k.parquet"
DONOR_SUMMARY = AUDIT / "hmba-bg-retry_donor_summary.tsv"


def main():
    print("[1/3] loading union peaks")
    peak_index, all_peak_nums = load_union_peaks()
    n_peaks_global = len(all_peak_nums)

    h5ads = sorted(RETRY_DIR.glob("*.h5ad"))
    print(f"[2/3] processing {len(h5ads)} retry h5ads")

    all_rows = []
    donor_rows = []
    for i, local in enumerate(h5ads, 1):
        hashname = local.stem
        size_gb = local.stat().st_size / 1e9
        t0 = time.time()
        print(f"[{i}/{len(h5ads)}] {hashname} ({size_gb:.2f} GB)", flush=True)
        try:
            major, M, nnz, n_cells = process_h5ad(local, peak_index, n_peaks_global, hashname)
        except Exception as e:
            print(f"  ERROR {hashname}: {e}", flush=True)
            continue
        pass_mask = nnz >= NNZ_FLOOR
        kept = int(pass_mask.sum())
        print(f"  cells={n_cells} pass_nnz1k={kept} time={time.time()-t0:.0f}s", flush=True)
        local.unlink()  # free disk immediately

        donor_id = f"hmba-bg-{hashname[:12]}"
        group_counts = defaultdict(lambda: np.zeros(n_peaks_global, dtype=np.int32))
        group_n_cells = defaultdict(int)
        group_pre_cells = defaultdict(int)
        for maj in np.unique(major):
            group_pre_cells[maj] += int((major == maj).sum())
        if kept > 0:
            kept_idx = np.where(pass_mask)[0]
            for maj in np.unique(major[kept_idx]):
                mask = (major[kept_idx] == maj)
                idx = kept_idx[mask]
                n_access = np.asarray(M[idx].sum(axis=0)).ravel().astype(np.int32)
                group_counts[maj] += n_access
                group_n_cells[maj] += int(mask.sum())
        for maj, arr in group_counts.items():
            n_total = group_n_cells[maj]
            if n_total == 0: continue
            active = np.where(arr >= MIN_CELLS_PER_PEAK)[0]
            if len(active) == 0: continue
            pn_active = all_peak_nums[active]
            for pi, pn in zip(active, pn_active):
                cnt = int(arr[pi])
                all_rows.append({
                    "peak_num": int(pn), "major_class": maj, "donor_id": donor_id,
                    "n_cells_access_direct": cnt, "n_cells_total_direct": int(n_total),
                    "fraction_open_direct": cnt / n_total, "study_id": STUDY_ID,
                })
            donor_rows.append({
                "donor_id": donor_id, "major_class": maj,
                "n_cells_pre": group_pre_cells.get(maj, 0),
                "n_cells_pass_nnz1k": n_total,
                "pct_pass": 100.0 * n_total / max(group_pre_cells.get(maj, 0), 1),
                "study_id": STUDY_ID, "h5ad_hash": hashname,
            })
        gc.collect()

    print(f"[3/3] writing {len(all_rows):,} rows to {OUT_PQ}")
    OUT_PQ.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(all_rows).write_parquet(OUT_PQ)
    DONOR_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    pl.DataFrame(donor_rows).write_csv(DONOR_SUMMARY, separator="\t")
    print(f"donors recovered: {len({r['donor_id'] for r in donor_rows})}")


if __name__ == "__main__":
    main()
