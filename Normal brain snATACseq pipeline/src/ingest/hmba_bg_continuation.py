#!/usr/bin/env python3
"""HMBA-BG batch runner — iterates the queue of 188 donors (ascending size), appends new
rows to `hmba-bg_direct_nnz1k.parquet` after every donor + checkpoints `*_donor_summary.tsv`.
A kill mid-run loses at most the in-flight donor.

Downloads via `aws s3 cp --no-sign-request` and deletes raw h5ad after processing.

Reconstructed from conversation transcript 2026-10-27.
"""
from __future__ import annotations
import gc, sys, time, traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from hmba_bg_process import (
    load_union_peaks, process_donor,
    OUT_DIR, AUDIT, DL_DIR, STUDY_ID, MIN_CELLS_PER_PEAK,
)
import numpy as np
import polars as pl

QUEUE = Path(__file__).parent / "hmba_queue.txt"  # ascending-size queue, "<size>\t<hash>" lines
OUT = OUT_DIR / "hmba-bg_direct_nnz1k.parquet"
SUM = AUDIT / "hmba-bg_donor_summary.tsv"
LOG = Path(__file__).parent / "hmba-bg_continuation.log"

CORE_SUM_COLS = [
    "donor_id", "major_class", "n_cells_pre", "n_cells_pass_nnz1k", "pct_pass", "study_id", "h5ad_hash",
]


def log(msg: str) -> None:
    ts = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def read_queue():
    with open(QUEUE) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) != 2: continue
            yield int(parts[0]), parts[1]


def current_donor_set() -> set[str]:
    if not OUT.exists(): return set()
    return set(pl.scan_parquet(OUT).select("donor_id").unique().collect()["donor_id"].to_list())


def main():
    AUDIT.mkdir(parents=True, exist_ok=True); OUT_DIR.mkdir(parents=True, exist_ok=True)
    DL_DIR.mkdir(parents=True, exist_ok=True)

    log("loading union peaks")
    peak_index, all_peak_nums = load_union_peaks()
    n_peaks_global = len(all_peak_nums)

    already = current_donor_set()
    log(f"existing donors in parquet: {len(already)}")

    queue = list(read_queue())
    todo = [(sz, h) for (sz, h) in queue if f"hmba-bg-{h[:12]}" not in already]
    log(f"queue total={len(queue)} to_process={len(todo)}")

    t_start = time.time(); done = 0; skipped = 0
    for i, (size, hashname) in enumerate(todo, 1):
        log(f"[{i}/{len(todo)}] {hashname} ({size/1e9:.2f} GB)")
        try:
            donor_id, gc_, gn_, gp_ = process_donor(hashname, peak_index, n_peaks_global)
        except Exception as e:
            log(f"  FAIL: {e}")
            skipped += 1
            traceback.print_exc()
            continue
        new_rows = []
        new_sum = []
        for maj, arr in gc_.items():
            n_total = gn_[maj]
            if n_total == 0: continue
            active = np.where(arr >= MIN_CELLS_PER_PEAK)[0]
            if len(active) == 0: continue
            pn_active = all_peak_nums[active]
            for pi, pn in zip(active, pn_active):
                cnt = int(arr[pi])
                new_rows.append({
                    "peak_num": int(pn), "major_class": maj, "donor_id": donor_id,
                    "n_cells_access_direct": cnt, "n_cells_total_direct": int(n_total),
                    "fraction_open_direct": cnt / n_total, "study_id": STUDY_ID,
                })
            new_sum.append({
                "donor_id": donor_id, "major_class": maj,
                "n_cells_pre": gp_.get(maj, 0),
                "n_cells_pass_nnz1k": n_total,
                "pct_pass": 100.0 * n_total / max(gp_.get(maj, 0), 1),
                "study_id": STUDY_ID, "h5ad_hash": hashname,
            })
        # Append to parquet by concatenating existing + new, then writing. Simple + atomic.
        if new_rows:
            new_df = pl.DataFrame(new_rows)
            if OUT.exists():
                existing = pl.read_parquet(OUT)
                combined = pl.concat([existing, new_df.select(existing.columns)])
            else:
                combined = new_df
            combined.write_parquet(OUT)
        # Append to donor summary TSV
        if new_sum:
            sum_df = pl.DataFrame(new_sum).select(CORE_SUM_COLS)
            header = not SUM.exists()
            with open(SUM, "a") as f:
                if header:
                    f.write("\t".join(CORE_SUM_COLS) + "\n")
                for row in sum_df.iter_rows():
                    f.write("\t".join(str(x) for x in row) + "\n")
        done += 1
        elapsed = (time.time() - t_start) / 60
        rate = done / elapsed * 60 if elapsed > 0 else 0
        eta = (len(todo) - i) / rate if rate > 0 else 0
        total_rows = pl.scan_parquet(OUT).select(pl.len()).collect().item()
        log(f"  OK: {donor_id} +{len(new_rows):,} rows (parquet={total_rows:,} rows); "
            f"done={done} skipped={skipped} elapsed={elapsed:.1f}min rate={rate:.1f}/hr "
            f"eta_rest={eta:.1f}hr")
        gc.collect()

    log(f"finished. done={done} skipped={skipped}")


if __name__ == "__main__":
    main()
