#!/usr/bin/env python3
"""Add amplicon_adjacent boolean column to every GBM atlas candidate table.

Rule (per project-gbm-preprint-results-template memory):
  amplicon_adjacent = TRUE if the peak lies within ±1 Mb of any known GBM
  CNV-amplified locus. Four loci are flagged:
    chr7:54,500,000-56,000,000   (EGFR amplicon)
    chr12:68,000,000-70,000,000  (MDM2 amplicon)
    chr8:127,000,000-129,000,000 (MYC amplicon)
    chr9:21,000,000-22,500,000   (CDKN2A deletion)
  ±1 Mb buffer applied symmetrically.

Why: at amplified loci, non-malignant pseudobulks show elevated apparent
accessibility consistent with ambient-fragment contamination from amplified
malignant cells in the same tissue. Flag is informational — rows remain in
all tables — but non-malignant AAV candidates adjacent to these loci are
to be deprioritized for AAV lead selection, as surfaced per-reviewer.

Usage:
    python3 add_amplicon_adjacent_flag.py <report_dir>

Default: processes every .csv and .parquet under the dir that has
chrom + start + end columns. Adds the column if not already present.
Preserves original files via <name>.bak before overwriting.
"""
import argparse
import os
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

# ±1 Mb buffer regions
AMPLICONS = [
    {"name": "EGFR",   "chrom": "chr7",  "start":  53_500_000, "end":  57_000_000},  # 54.5-56 ± 1 Mb
    {"name": "MDM2",   "chrom": "chr12", "start":  67_000_000, "end":  71_000_000},  # 68-70 ± 1 Mb
    {"name": "MYC",    "chrom": "chr8",  "start": 126_000_000, "end": 130_000_000},  # 127-129 ± 1 Mb
    {"name": "CDKN2A", "chrom": "chr9",  "start":  20_000_000, "end":  23_500_000},  # 21-22.5 ± 1 Mb
]


def classify_peak(chrom: str, start: int, end: int):
    """Return (amplicon_adjacent: bool, amplicon_adjacent_locus: str)."""
    if not isinstance(chrom, str):
        return False, ""
    for amp in AMPLICONS:
        if chrom != amp["chrom"]:
            continue
        # Interval overlap: peak [start, end) overlaps buffer [amp.start, amp.end)
        if int(start) < amp["end"] and int(end) > amp["start"]:
            return True, amp["name"]
    return False, ""


def find_coord_columns(df: pd.DataFrame):
    """Return (chrom_col, start_col, end_col) if present, else (None, None, None)."""
    cols = {c.lower(): c for c in df.columns}
    chrom = cols.get("chrom") or cols.get("chr") or cols.get("chromosome")
    start = cols.get("start") or cols.get("peak_start") or cols.get("startpos")
    end = cols.get("end") or cols.get("peak_end") or cols.get("endpos") or cols.get("stop")
    return chrom, start, end


def try_parse_peak_id(df: pd.DataFrame):
    """Fall-back: parse chrom/start/end from peak_id if coord columns missing.
    Supports formats like 'chr5:71864459-71866513' or 'chr5-71864459-71866513'.
    Returns (series_chrom, series_start, series_end) or (None, None, None).
    """
    if "peak_id" not in df.columns:
        return None, None, None
    pat = re.compile(r"^(chr[\w]+)[:_-](\d+)[-_](\d+)$")
    parts = df["peak_id"].astype(str).str.extract(pat)
    if parts.isna().all().all():
        return None, None, None
    chrom = parts[0]
    start = pd.to_numeric(parts[1], errors="coerce").astype("Int64")
    end = pd.to_numeric(parts[2], errors="coerce").astype("Int64")
    return chrom, start, end


def process_file(path: Path, dry_run: bool = False) -> dict:
    """Process a single CSV or parquet. Returns a stats dict."""
    stats = {"path": str(path), "status": "skipped", "n_flagged": 0, "total_rows": 0, "reason": ""}

    try:
        if path.suffix == ".csv":
            df = pd.read_csv(path, low_memory=False)
        elif path.suffix == ".parquet":
            df = pd.read_parquet(path)
        else:
            stats["reason"] = "unsupported extension"
            return stats
    except Exception as e:
        stats["reason"] = f"read failed: {e}"
        return stats

    stats["total_rows"] = len(df)

    if "amplicon_adjacent" in df.columns:
        stats["status"] = "already_has_column"
        stats["reason"] = "amplicon_adjacent column already exists; skipping"
        return stats

    chrom_col, start_col, end_col = find_coord_columns(df)
    if not (chrom_col and start_col and end_col):
        # Try parse from peak_id
        chrom_s, start_s, end_s = try_parse_peak_id(df)
        if chrom_s is None:
            stats["reason"] = "no chrom/start/end columns and no parseable peak_id"
            return stats
        chrom_vals = chrom_s
        start_vals = start_s
        end_vals = end_s
    else:
        chrom_vals = df[chrom_col]
        start_vals = df[start_col]
        end_vals = df[end_col]

    # Vectorized flag computation (fast)
    flags = []
    loci = []
    for ch, st, en in zip(chrom_vals, start_vals, end_vals):
        try:
            is_adj, locus = classify_peak(str(ch), int(st), int(en))
        except (ValueError, TypeError):
            is_adj, locus = False, ""
        flags.append(is_adj)
        loci.append(locus)

    df["amplicon_adjacent"] = flags
    df["amplicon_adjacent_locus"] = loci
    stats["n_flagged"] = sum(flags)

    if dry_run:
        stats["status"] = "dry_run"
        return stats

    # Backup and write
    bak = path.with_suffix(path.suffix + ".bak")
    if not bak.exists():
        shutil.copy2(path, bak)

    if path.suffix == ".csv":
        df.to_csv(path, index=False)
    elif path.suffix == ".parquet":
        df.to_parquet(path, index=False)

    stats["status"] = "updated"
    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("report_dir", help="root directory to scan for candidate tables")
    ap.add_argument("--dry-run", action="store_true", help="report what would change without writing")
    ap.add_argument("--exts", default=".csv,.parquet", help="comma-separated extensions to process")
    args = ap.parse_args()

    root = Path(args.report_dir)
    if not root.exists():
        print(f"ERROR: {root} does not exist", file=sys.stderr)
        sys.exit(1)

    exts = set(args.exts.split(","))
    all_files = []
    for ext in exts:
        all_files.extend(root.rglob(f"*{ext}"))
    # Dedupe + sort
    all_files = sorted(set(all_files))
    # Skip .bak files
    all_files = [p for p in all_files if ".bak" not in p.name]

    print(f"Scanning {len(all_files)} files under {root}")
    print(f"Dry run: {args.dry_run}")
    print()

    totals = {"updated": 0, "already_has_column": 0, "skipped": 0, "dry_run": 0, "total_flagged": 0}
    for p in all_files:
        stats = process_file(p, dry_run=args.dry_run)
        status = stats["status"]
        rel = p.relative_to(root)
        if status == "updated" or status == "dry_run":
            print(f"  [{status:12s}] {stats['n_flagged']:>5d} / {stats['total_rows']:>7,} flagged  {rel}")
            totals["total_flagged"] += stats["n_flagged"]
            totals[status] = totals.get(status, 0) + 1
        elif status == "already_has_column":
            print(f"  [has_column ] already has amplicon_adjacent — skipped  {rel}")
            totals["already_has_column"] += 1
        else:
            print(f"  [skipped    ] {stats['reason']}  {rel}")
            totals["skipped"] += 1

    print()
    print(f"Summary: updated={totals['updated']}  already_had={totals['already_has_column']}  skipped={totals['skipped']}  dry_run={totals['dry_run']}")
    print(f"Total amplicon-adjacent rows flagged across all files: {totals['total_flagged']:,}")


if __name__ == "__main__":
    main()
