"""Per-cohort fragment pileup for any region of the atlas.

Usage:
    python pileup_region.py chr5:44411500-44412700 --out pileup.png

Behavior:
    - Resolves each cohort's raw tabix-indexed fragment files on a remote host
      (via find over SSH)
    - If an index .tbi lives in a sibling dir (GDC packaging), symlinks it
      into place for tabix to find, then cleans up
    - Queries fragments overlapping the region per cohort
    - Accumulates per-bp coverage
    - Optionally plots stacked tracks to PNG

Configuration — the raw-fragment store lives on the pipeline host (not in this
repo; the data/matrix/ parquets at Zenodo are the derivatives). Point this
script at your host by setting these environment variables:

    ATLAS_PILEUP_SSH   full ssh command prefix, e.g.
                       'ssh -i ~/.ssh/my-key.pem ubuntu@1.2.3.4'
    ATLAS_PILEUP_CONDA (optional) conda activation command run before each
                       tabix query, e.g. 'source ~/miniconda3/etc/profile.d/conda.sh && conda activate atacseq'
    ATLAS_PILEUP_RAW_ROOT  (optional) root dir on remote host containing the
                           per-cohort subdirs listed in COHORT_FRAG_ROOTS below;
                           defaults to /data/projects/atacseq/raw

Only `sundaram_2024` and `sojka_2025` shipped raw fragments publicly in their
deposits; the other 6 cohorts released processed h5ad matrices only, so this
script's output is limited to those 2 cohorts regardless of how it is configured.
"""
from __future__ import annotations
import argparse
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

SSH = os.environ.get("ATLAS_PILEUP_SSH", "").strip()
CONDA_ACTIVATE = os.environ.get("ATLAS_PILEUP_CONDA", "")
RAW_ROOT = os.environ.get("ATLAS_PILEUP_RAW_ROOT", "/data/projects/atacseq/raw").rstrip("/")

COHORT_FRAG_ROOTS = {
    "sundaram_2024":  (f"{RAW_ROOT}/sundaram_gbm_gdc",        "scATAC_GBMx_*.fragments.tsv.gz"),
    "sojka_2025":     (f"{RAW_ROOT}/gse276177_khan_astro",    "GSM*_fragments.tsv.gz"),
    "greenwald_2024": (f"{RAW_ROOT}/gbm_space",               "*fragments*.tsv.gz"),
    "guo_2020":       (f"{RAW_ROOT}/gse138794_guo",           "*fragments*.tsv.gz"),
    "zhao_2024":      (f"{RAW_ROOT}/gbm_tme_atlas_hra004942", "*fragments*.tsv.gz"),
    "mathewson_2021": (f"{RAW_ROOT}/mathewson_lupien",        "*fragments*.tsv.gz"),
    "guilhamon_2021": (f"{RAW_ROOT}/guilhamon",               "*fragments*.tsv.gz"),
    "raviram_2023":   (f"{RAW_ROOT}/gse165037",               "*fragments*.tsv.gz"),
}


def parse_region(region: str) -> tuple[str, int, int]:
    m = re.match(r"^(chr[\w]+):([\d,]+)-([\d,]+)$", region)
    if not m:
        sys.exit(f"region must be 'chrN:start-end', got {region!r}")
    return m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))


def _ssh(cmd: str) -> str:
    if not SSH:
        sys.exit(
            "ATLAS_PILEUP_SSH is not set.\n"
            "This script queries raw fragment files on a remote host (not in this repo —\n"
            "the atlas deposit at Zenodo contains the derived pseudobulk matrices only).\n"
            "Export your ssh prefix before running, e.g.\n"
            "  export ATLAS_PILEUP_SSH='ssh -i ~/.ssh/my-key.pem ubuntu@your.host'\n"
            "See the module docstring for ATLAS_PILEUP_CONDA and ATLAS_PILEUP_RAW_ROOT."
        )
    full = f"{CONDA_ACTIVATE}; {cmd}" if CONDA_ACTIVATE else cmd
    r = subprocess.run(SSH.split() + [full], capture_output=True, text=True)
    return r.stdout


def _list_fragment_files(root: str, pattern: str) -> list[str]:
    out = _ssh(f'find "{root}" -type f -name "{pattern}" ! -name "*.tbi" ! -name "*.tbi.gz" 2>/dev/null | sort -u')
    return [p.strip() for p in out.splitlines() if p.strip()]


def tabix_region(frag_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """tabix query via SSH. Handles the GDC packaging where .tbi lives in a sibling dir."""
    # tabix expects frag.tsv.gz.tbi next to frag.tsv.gz. In GDC packaging, both are
    # inside sibling dirs with .tbi suffix. We symlink first, query, then remove.
    stem = frag_path
    tbi_candidates = [
        f"{stem}.tbi",                                  # expected
        f"{stem}.tbi.gz",                               # some pkg
        f"{Path(stem).parent}.tbi/{Path(stem).name}.tbi",  # GDC: sibling dir
    ]
    # Build a one-shot bash that creates a symlink in a safe tmp dir, runs tabix, cleans up
    cmd = f'''
      base="$(basename {shlex.quote(stem)})"
      tmp=$(mktemp -d)
      ln -sf {shlex.quote(stem)} "$tmp/$base"
      for cand in {" ".join(shlex.quote(c) for c in tbi_candidates)}; do
        if [ -e "$cand" ]; then
          # if gzipped tbi, decompress a copy
          case "$cand" in
            *.tbi.gz) zcat "$cand" > "$tmp/$base.tbi" ;;
            *) ln -sf "$cand" "$tmp/$base.tbi" ;;
          esac
          break
        fi
      done
      if [ -e "$tmp/$base.tbi" ]; then
        tabix "$tmp/$base" "{chrom}:{start}-{end}" 2>/dev/null
      fi
      rm -rf "$tmp"
    '''
    out = _ssh(cmd)
    rows = []
    for ln in out.splitlines():
        parts = ln.split("\t")
        if len(parts) < 3:
            continue
        try:
            cnt = int(parts[4]) if len(parts) > 4 and parts[4].isdigit() else 1
            rows.append((parts[0], int(parts[1]), int(parts[2]), parts[3] if len(parts) > 3 else "", cnt))
        except ValueError:
            continue
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "barcode", "count"])


def pileup(df: pd.DataFrame, chrom: str, start: int, end: int) -> np.ndarray:
    cov = np.zeros(end - start, dtype=np.int32)
    for _, r in df.iterrows():
        s = max(r["start"], start) - start
        e = min(r["end"], end) - start
        if e > s:
            cov[s:e] += r["count"]
    return cov


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("region")
    ap.add_argument("--out", type=Path, default=None, help="save stacked pileup plot to PNG")
    ap.add_argument("--cohorts", default="", help="comma-sep slugs; default = all 8")
    ap.add_argument("--pad", type=int, default=200, help="padding each side of region for context")
    args = ap.parse_args()

    chrom, start, end = parse_region(args.region)
    start = max(0, start - args.pad)
    end = end + args.pad
    cohorts = args.cohorts.split(",") if args.cohorts else list(COHORT_FRAG_ROOTS)

    print(f"Region: {chrom}:{start}-{end}  ({end-start} bp)")
    per_cov = {}
    per_n = {}
    for c in cohorts:
        root, pat = COHORT_FRAG_ROOTS[c]
        files = _list_fragment_files(root, pat)
        if not files:
            print(f"  {c}: NO fragment files under {root}")
            continue
        print(f"  {c}: {len(files)} files; querying…")
        cov = np.zeros(end - start, dtype=np.int32)
        n = 0
        for f in files:
            df = tabix_region(f, chrom, start, end)
            if len(df):
                cov += pileup(df, chrom, start, end)
                n += len(df)
        per_cov[c] = cov
        per_n[c] = n
        print(f"    {n:,} fragments; max depth {cov.max()}")

    if args.out and per_cov:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("matplotlib not installed; skipping plot")
            return 0
        fig, axes = plt.subplots(len(per_cov), 1, figsize=(10, 1.3 * len(per_cov)),
                                 sharex=True, constrained_layout=True)
        if len(per_cov) == 1:
            axes = [axes]
        x = np.arange(start, end)
        for ax, (c, cov) in zip(axes, per_cov.items()):
            ax.fill_between(x, cov, step="pre", linewidth=0)
            ax.set_ylabel(f"{c}\n({per_n[c]:,} frags)", fontsize=7, rotation=0, ha="right", va="center")
            ax.tick_params(labelsize=7)
        axes[-1].set_xlabel(f"{chrom} position")
        fig.suptitle(f"{chrom}:{start}-{end} per-cohort fragment pileup", fontsize=10)
        fig.savefig(args.out, dpi=150)
        print(f"wrote {args.out}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
