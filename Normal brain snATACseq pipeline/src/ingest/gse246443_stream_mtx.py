#!/usr/bin/env python3
"""GSE246443 MDD dlPFC — stream the 31.5 GB uncompressed cell_peak_counts.mtx cell-by-cell
from the tar.gz without extracting (saves disk). ~50 min.

Input is a MatrixMarket .mtx inside raw_data_GEO_cluster_peak/ in the GEO deposit tarball:
  825,492 peaks × 200,071 cells, 2,072,292,226 nz entries
The file is sorted by column 2 (cell_idx), so processing one cell at a time with modest
memory works. Cell metadata (Subject, condition, Celltypes, etc.) comes from the cluster
metadata.csv — we only ingest cells with condition='control'.

Reconstructed from conversation transcript 2026-10-27.
"""
import csv, subprocess, sys, time, io
import numpy as np
import pandas as pd
from pathlib import Path
import pyranges as pr
import polars as pl

TAR     = "/tmp/geo246443_work/cluster_peak.tar.gz"
MEMBER  = "raw_data_GEO_cluster_peak/cluster_peak_counts.mtx"
META    = "/tmp/geo246443_work/raw_data_GEO_cluster_peak/cluster_metadata.csv"
FEAT    = "/tmp/geo246443_work/raw_data_GEO_cluster_peak/cluster_features.txt"
UNION   = "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/references/union_peaks.bed"
OUT_PQ  = "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0/h5ads_v2.2_nnz_rebuild/geo-GSE246443_direct_nnz1k.parquet"
STUDY_ID = "geo-GSE246443"

NNZ_FLOOR = 1000
MIN_CELLS_PER_PEAK = 3
BATCH = 50  # chunk size for progress logging

# Author's broad Celltypes column → atlas major_class
CT_MAP = {
    "ExN": "Excitatory", "InN": "Inhibitory",
    "Ast": "Astrocyte", "Oli": "Oligodendrocyte",
    "OPC": "OPC", "Mic": "Microglia", "End": "Endothelial",
    # Mix1/2 unmapped
}

canon = set([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])

print(f"[{time.strftime('%H:%M:%S')}] load union peaks", flush=True)
union = pd.read_csv(UNION, sep="\t", header=None, names=["chrom", "start", "end", "name"])
union = union[union.chrom.isin(canon)].reset_index(drop=True)
union["peak_num"] = union["name"].str.replace("peak_", "", regex=False).astype(int)
union["ord"] = np.arange(len(union))

print(f"[{time.strftime('%H:%M:%S')}] load features + overlap author → union", flush=True)
feats = pd.read_csv(FEAT, sep="\t", header=None, names=["name"])
feat_parts = feats["name"].str.extract(r"^(chr[0-9XY]+)[-](\d+)[-](\d+)$")
feat_parts.columns = ["chrom", "start", "end"]
feat_parts["start"] = feat_parts["start"].astype("Int64").astype(int)
feat_parts["end"] = feat_parts["end"].astype("Int64").astype(int)
feat_parts["author_idx"] = np.arange(1, len(feats) + 1)
usable = feat_parts.dropna(subset=["start", "end"])
usable = usable[usable.chrom.isin(canon)].copy()
a = pr.PyRanges(usable.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
u = pr.PyRanges(union.rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"}))
j = a.join(u, how=None).df
author_to_union = {}
for ai, uo in zip(j["author_idx"].values, j["ord"].values):
    author_to_union.setdefault(int(ai), []).append(int(uo))
print(f"  mapped {len(author_to_union):,} author peaks → union", flush=True)

print(f"[{time.strftime('%H:%M:%S')}] load cell metadata", flush=True)
md = pd.read_csv(META)
md.columns = ["barcode"] + list(md.columns[1:])
md["major_class"] = md["Celltypes"].map(CT_MAP)
# Keep only MDD controls (condition == 'control') with mapped class
md["atlas_use"] = (md["condition"] == "control") & md["major_class"].notna()
md["cell_idx"] = np.arange(1, len(md) + 1)
md["donor_id"] = md["Subject"]  # already unique per subject × sex
cell_donor = md["donor_id"].astype(str).values
cell_class = md["major_class"].astype(str).values
cell_keep = md["atlas_use"].values
print(f"  controls to ingest: {cell_keep.sum():,} / {len(md):,}", flush=True)
print(f"  PASS donors: {md.loc[cell_keep,'donor_id'].nunique()}", flush=True)

# Group ids: (donor, class) → batch index. 0 = not ingested.
groups = sorted(set(zip(cell_donor[cell_keep], cell_class[cell_keep])))
group_to_idx = {g: i + 1 for i, g in enumerate(groups)}
cell_group = np.zeros(len(md) + 2, dtype=np.int32)
for i in np.where(cell_keep)[0]:
    cell_group[i + 1] = group_to_idx[(cell_donor[i], cell_class[i])]
n_groups = len(groups)
print(f"  n (donor, class) groups: {n_groups:,}", flush=True)

# Per-group accessibility accumulator: dense (n_groups+1) × len(union)
accum = np.zeros((n_groups + 1, len(union)), dtype=np.int32)
n_cells_per_group = np.zeros(n_groups + 1, dtype=np.int32)
print(f"  allocated accumulator: {accum.nbytes/1e9:.2f} GB", flush=True)

print(f"[{time.strftime('%H:%M:%S')}] stream mtx", flush=True)
proc = subprocess.Popen(["tar", "-xOzf", TAR, MEMBER],
                        stdout=subprocess.PIPE, bufsize=1024 * 1024)
rd = io.TextIOWrapper(proc.stdout, encoding="utf-8", errors="replace")

cur_cell = -1
cur_union_peaks = set()
cells_proc = cells_kept = 0
dims_read = False
nnz = 0
line_count = 0
last_report = time.time()
t_start = last_report

for line in rd:
    if not dims_read:
        if line.startswith("%"): continue
        parts = line.split()
        if len(parts) == 3:
            n_rows, n_cols, nnz = (int(x) for x in parts)
            print(f"  mtx dims: {n_rows:,} peaks × {n_cols:,} cells, {nnz:,} nz", flush=True)
            dims_read = True
            continue
    parts = line.split()
    if len(parts) != 3: continue
    peak_idx, cell_idx, cnt = int(parts[0]), int(parts[1]), int(parts[2])
    line_count += 1
    if cell_idx != cur_cell:
        if cur_cell > 0:
            gi = cell_group[cur_cell] if cur_cell < len(cell_group) else 0
            if gi > 0 and len(cur_union_peaks) >= NNZ_FLOOR:
                for uo in cur_union_peaks:
                    accum[gi, uo] += 1
                n_cells_per_group[gi] += 1
                cells_kept += 1
            cells_proc += 1
        cur_cell = cell_idx
        cur_union_peaks = set()
    if cnt > 0:
        for uo in author_to_union.get(peak_idx, ()):
            cur_union_peaks.add(uo)
    if time.time() - last_report > 30:
        pct = 100.0 * line_count / max(nnz, 1)
        rate = line_count / (time.time() - t_start) if time.time() > t_start else 0
        print(f"  [{time.strftime('%H:%M:%S')}] lines={line_count:,}/{nnz:,} ({pct:.1f}%) "
              f"cells_proc={cells_proc:,} cells_kept={cells_kept:,} rate={rate/1e6:.1f}M/s", flush=True)
        last_report = time.time()

# Finalize last cell
if cur_cell > 0:
    gi = cell_group[cur_cell] if cur_cell < len(cell_group) else 0
    if gi > 0 and len(cur_union_peaks) >= NNZ_FLOOR:
        for uo in cur_union_peaks:
            accum[gi, uo] += 1
        n_cells_per_group[gi] += 1
        cells_kept += 1
    cells_proc += 1
proc.stdout.close()
proc.wait()
print(f"[{time.strftime('%H:%M:%S')}] stream done — cells_proc={cells_proc:,} cells_kept={cells_kept:,}",
      flush=True)

rows = []
for gi, (donor, cls) in enumerate(groups, start=1):
    n_total = int(n_cells_per_group[gi])
    if n_total == 0: continue
    access = accum[gi, :]
    active = np.where(access >= MIN_CELLS_PER_PEAK)[0]
    if len(active) == 0: continue
    rows.append(pd.DataFrame({
        "peak_num": union["peak_num"].values[active],
        "major_class": cls,
        "donor_id": donor,
        "n_cells_access_direct": access[active].astype(np.int32),
        "n_cells_total_direct": n_total,
        "fraction_open_direct": access[active].astype(np.float64) / n_total,
        "study_id": STUDY_ID,
    }))
out = pd.concat(rows, ignore_index=True)
out.to_parquet(OUT_PQ, index=False)
print(f"wrote {OUT_PQ} ({len(out):,} rows, {out['donor_id'].nunique()} donors)")
