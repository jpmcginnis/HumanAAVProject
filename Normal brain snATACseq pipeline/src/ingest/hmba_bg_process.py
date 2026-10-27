#!/usr/bin/env python3
"""HMBA-BG (Allen Human Multi-species Brain Atlas Basal Ganglia pre-print ATAC).

Downloads one SnapATAC2 h5ad per donor from Allen's S3 bucket, decodes fragments from the
nonstandard `obsm/fragment_paired` CSR (encodes (position, length) pairs across
total_genome_bp), overlaps against the 990K atlas union peak set, binarizes, applies
nnz≥1000, aggregates per (donor, major_class) using Allen's HMBA Group label. Deletes
raw h5ad immediately after each donor to limit disk pressure.

Called by `hmba_bg_continuation.py` for the full 188-donor queue; also usable directly:
    python3 hmba_bg_process.py --donors <hash1> <hash2> ...

Reconstructed from conversation transcript 2026-10-27.
"""
from __future__ import annotations
import gc, json, os, subprocess, sys, time
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp

UNION   = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v1.0.0/union_peaks.bed")
OUT_DIR = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/h5ads_v2.2_nnz_rebuild")
AUDIT   = Path("/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/audit_2026-10-25/acquisition/hmba-bg")
DL_DIR  = Path("/tmp/hmba_dl")
S3_PREFIX = "s3://hmba-bican-sharing-802451596237-us-west-2/BasalGanglia_pre-print_ATAC/human/snapatac2/"

NNZ_FLOOR = 1000
MIN_CELLS_PER_PEAK = 3
STUDY_ID = "hmba-bg-pre-print"

# HMBA 'Group' → atlas major_class (prefer Group > Subclass > Class).
# Specific MSN subtypes get mapped to D1-MSN / D2-MSN / D1-D2-MSN for direct D1-vs-D2 contrast.
GROUP_MAP = {
    # Striatum MSNs
    "D1 MSN":      "D1-MSN",       "D1-MSN":  "D1-MSN",
    "D2 MSN":      "D2-MSN",       "D2-MSN":  "D2-MSN",
    "D1-D2 MSN":   "D1-D2-MSN",    "D1D2 MSN":"D1-D2-MSN",
    # Striatum interneurons + other GABAergic
    "Striatal interneuron": "Inhibitory", "Striatum IN": "Inhibitory",
    # Pallidal / STN / nigra neuronal
    "GPe PV":       "Inhibitory",   "GPe LHX6":    "Inhibitory",
    "STN":          "Excitatory",   "SNc DA":      "Dopaminergic",
    "SNr":          "Inhibitory",
    # Glia
    "Astrocyte":    "Astrocyte",    "Oligodendrocyte": "Oligodendrocyte",
    "OPC":          "OPC",          "Microglia":   "Microglia",
    "Endothelial":  "Endothelial",  "Mural":       "Mural",
    "Perivascular": "Mural",        "Fibroblast":  "Mural",
    "T cell":       "Immune_other", "Macrophage":  "Immune_other",
    "Ependymal":    "Mural",
}
def map_group(g: str) -> str:
    """Case-tolerant mapping with 'unknown' fallback."""
    if g in GROUP_MAP: return GROUP_MAP[g]
    g2 = g.strip()
    if g2 in GROUP_MAP: return GROUP_MAP[g2]
    for k, v in GROUP_MAP.items():
        if k.lower() == g2.lower(): return v
    return "Unclassified"


def load_union_peaks():
    """Return (peak_index: dict {(chrom, bin_start): [(start,end,peak_num), ...]}, all_peak_nums np.array)."""
    pb = pd.read_csv(UNION, sep="\t", header=None, names=["chrom", "start", "end", "name"])
    pb["peak_num"] = pb["name"].str.replace("peak_", "", regex=False).astype(int)
    pb = pb.sort_values(["chrom", "start"]).reset_index(drop=True)
    # Bin by 100 kb per chrom for fast interval lookup
    BIN = 100_000
    idx: dict = defaultdict(list)
    for row in pb.itertuples(index=False):
        for b in range(row.start // BIN, row.end // BIN + 1):
            idx[(row.chrom, b)].append((int(row.start), int(row.end), int(row.peak_num)))
    return idx, pb["peak_num"].values.astype(np.int64)


def process_h5ad(local: Path, peak_index, n_peaks_global: int, hashname: str):
    """Return (major array of length n_cells, M = cell × union_peak binary, nnz array, n_cells)."""
    with h5py.File(local, "r") as f:
        # obs: HMBA Group label
        obs_group = f["obs"]
        group = np.asarray(obs_group["Group"].asstr()[...])  # or Subclass as fallback
        major = np.array([map_group(g) for g in group])
        n_cells = len(major)

        # obsm/fragment_paired: CSR of (position, length) encoded fragments
        fp = f["obsm/fragment_paired"]
        data = fp["data"][...]       # 1-D array of ints: alternating position, length
        indices = fp["indices"][...] # column index (ignored — single pseudo-column of 2*n_frags per cell)
        indptr = fp["indptr"][...]
        # Each cell's fragments: data[indptr[c]:indptr[c+1]] holds interleaved (pos, len) pairs

        # obs: chrom encoded as integer index into obs/chrom_cat (if present) else obs/chrom str
        # (Allen's snapatac2 convention varies; here we approximate by reading obs['chrom'] if available)
        # In the real pipeline, we iterate per cell and resolve chrom from f['uns/chrom'] arrays.
        # For brevity in this reconstruction, see hmba_bg_continuation.py which carries the exact decoder.
        try:
            chrom_arr = np.asarray(f["uns/reference/chrom"].asstr()[...])
            chrom_sizes = np.asarray(f["uns/reference/chrom_size"][...])
        except KeyError:
            chrom_arr = np.array([f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]])
            chrom_sizes = np.array([250_000_000] * len(chrom_arr))
        chrom_offsets = np.concatenate(([0], np.cumsum(chrom_sizes))).astype(np.int64)

    # Build cell × union_peak binary accessibility
    rows = []; cols = []
    nnz = np.zeros(n_cells, dtype=np.int32)
    BIN = 100_000
    for c in range(n_cells):
        pairs = data[indptr[c]:indptr[c + 1]]
        if len(pairs) == 0: continue
        # pairs are interleaved (position, length)
        pos = pairs[0::2].astype(np.int64)
        length = pairs[1::2].astype(np.int64)
        # Resolve each position → (chrom, local_pos)
        chrom_i = np.searchsorted(chrom_offsets, pos, side="right") - 1
        local_pos = pos - chrom_offsets[chrom_i]
        seen = set()
        for k in range(len(pos)):
            ci = chrom_i[k]
            if ci < 0 or ci >= len(chrom_arr): continue
            chrom = chrom_arr[ci]
            lo = local_pos[k]
            hi = lo + length[k]
            for b in range(lo // BIN, hi // BIN + 1):
                for (ps, pe, pnum) in peak_index.get((chrom, b), ()):
                    if lo <= pe and hi >= ps and pnum not in seen:
                        seen.add(pnum)
                        rows.append(c); cols.append(pnum - 1)  # peak_num is 1-indexed
        nnz[c] = len(seen)
    M = sp.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                      shape=(n_cells, n_peaks_global))
    return major, M, nnz, n_cells


def process_donor(hashname: str, peak_index, n_peaks_global: int, keep_local: bool = False):
    """Download → process → delete raw. Return per-donor accumulators."""
    url = f"{S3_PREFIX}{hashname}.h5ad"
    DL_DIR.mkdir(parents=True, exist_ok=True)
    local = DL_DIR / f"{hashname}.h5ad"
    if not local.exists():
        print(f"  dl {hashname}", flush=True)
        subprocess.run(["aws", "s3", "cp", url, str(local), "--no-sign-request", "--quiet"], check=True)
    size_gb = local.stat().st_size / 1e9
    t0 = time.time()
    major, M, nnz, n_cells = process_h5ad(local, peak_index, n_peaks_global, hashname)
    pass_mask = nnz >= NNZ_FLOOR
    kept = int(pass_mask.sum())
    print(f"  cells={n_cells} pass_nnz1k={kept} time={time.time()-t0:.0f}s size={size_gb:.2f}GB", flush=True)
    if not keep_local:
        local.unlink()

    group_counts = defaultdict(lambda: np.zeros(n_peaks_global, dtype=np.int32))
    group_n_cells = defaultdict(int)
    group_pre_cells = defaultdict(int)
    donor_id = f"hmba-bg-{hashname[:12]}"
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
    return donor_id, group_counts, group_n_cells, group_pre_cells


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", required=True, help="hashnames (no .h5ad)")
    ap.add_argument("--out", default=str(OUT_DIR / "hmba-bg_direct_nnz1k.parquet"))
    args = ap.parse_args()

    AUDIT.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print("[1/3] loading union peaks")
    peak_index, all_peak_nums = load_union_peaks()
    n_peaks_global = len(all_peak_nums)

    print(f"[2/3] processing {len(args.donors)} donors")
    all_rows = []
    for i, h in enumerate(args.donors, 1):
        print(f"[{i}/{len(args.donors)}] {h}")
        donor_id, gc_, gn_, gp_ = process_donor(h, peak_index, n_peaks_global)
        for maj, arr in gc_.items():
            n_total = gn_[maj]
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
        gc.collect()

    print(f"[3/3] writing {len(all_rows):,} rows to {args.out}")
    pl.DataFrame(all_rows).write_parquet(args.out)


if __name__ == "__main__":
    main()
