"""Parallel version of requant_greenwald — spawn N workers each handling a slice
of the nnz range. HDF5 chunk-level gzip is per-chunk, so parallel readers of
different byte ranges actually get CPU-parallel decompression (each on its own
Python process).

Each worker writes a partial group_sums (n_groups x n_vars int32) to a .npy file.
Main script merges them, does overlap projection, writes pseudobulk.

Writes to a DIFFERENT output file so it doesn't clobber the running serial job:
  pseudobulk_v7_greenwald_2024_fast.parquet

If this finishes first, we rename it to pseudobulk_v7_greenwald_2024.parquet and
kill the serial run.
"""
import os
import sys
import time
import subprocess
import multiprocessing as mp
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

PILOT = Path("/data/projects/atacseq/pilot_v5")
V5_BED = str(PILOT / "v5_union_reference.bed")
OUT = PILOT / "requant"

ATAC_H5AD = str(PILOT / "gbm_space.h5ad")
RNA_H5AD = "/data/projects/atacseq/raw/gbm_space/GBM_space_snRNA.h5ad"
LOG = PILOT / "requant_greenwald_fast.log"
PARTIAL_DIR = PILOT / "_greenwald_parallel_partials"
N_WORKERS = 4
WORKER_CHUNK = 100_000_000  # 100M nnz per worker read


def log(m):
    s = f"[greenwald-fast] {m}"
    print(s, flush=True)
    with open(LOG, "a") as f:
        f.write(s + "\n")


def _h5col(group, key):
    d = group[key]
    if isinstance(d, h5py.Group):
        if "codes" in d.keys() and "categories" in d.keys():
            codes = d["codes"][:]
            cats_raw = d["categories"][:]
            cats = [c.decode() if isinstance(c, (bytes, bytearray)) else str(c) for c in cats_raw]
            return np.asarray([cats[c] if 0 <= c < len(cats) else "" for c in codes])
    arr = d[:]
    if arr.dtype.kind in ("S", "O"):
        return np.asarray([x.decode() if isinstance(x, (bytes, bytearray)) else str(x) for x in arr])
    return arr


def load_rna_annotations():
    log(f"opening {RNA_H5AD}")
    with h5py.File(RNA_H5AD, "r") as f:
        obs = f["obs"]
        idx_key = obs.attrs.get("_index", "index")
        if isinstance(idx_key, bytes): idx_key = idx_key.decode()
        cell_ids = _h5col(obs, idx_key) if idx_key in obs.keys() else _h5col(obs, "cell_id")
        donor = _h5col(obs, "donor_id") if "donor_id" in obs.keys() else np.array(["unknown"] * len(cell_ids))
        lab = _h5col(obs, "annotation_granular") if "annotation_granular" in obs.keys() else _h5col(obs, "annotation_coarse")
        status = _h5col(obs, "Cell_status") if "Cell_status" in obs.keys() else None
    labels = []
    for i in range(len(cell_ids)):
        is_mal = False
        if status is not None and str(status[i]).lower() in ("malignant", "tumor", "tumour"):
            is_mal = True
        if "malignant" in str(lab[i]).lower() or "tumor" in str(lab[i]).lower():
            is_mal = True
        if is_mal:
            l = str(lab[i])
            if any(s in l for s in ["MES", "NPC", "AC", "OPC", "Hypox", "Mes", "Prolif", "G1S", "G2M"]):
                labels.append(f"malignant_{l}")
            else:
                labels.append("malignant")
        else:
            labels.append(str(lab[i]))
    return pd.DataFrame({"cell_id": cell_ids, "donor_id": donor, "label": labels}).set_index("cell_id")


def worker(args):
    """Process one slice of nnz [i0, i1). Produces partial group_sums as .npy."""
    worker_id, nnz_start, nnz_end, cell_group_ids_bytes, n_groups, n_vars, indptr_bytes, partial_path = args
    cell_group_ids = np.frombuffer(cell_group_ids_bytes, dtype=np.int32)
    indptr = np.frombuffer(indptr_bytes, dtype=np.int64)
    partial = np.zeros((n_groups, n_vars), dtype=np.int32)
    t0 = time.time()
    chunks = 0
    i = nnz_start
    with h5py.File(ATAC_H5AD, "r", rdcc_nbytes=2 << 30, rdcc_nslots=100003) as f:
        data_ds = f["X"]["data"]
        indices_ds = f["X"]["indices"]
        while i < nnz_end:
            j = min(i + WORKER_CHUNK, nnz_end)
            d_chunk = data_ds[i:j]
            idx_chunk = indices_ds[i:j]
            entry_pos = np.arange(i, j, dtype=np.int64)
            cell_ids = np.searchsorted(indptr, entry_pos, side="right") - 1
            cell_ids = np.clip(cell_ids, 0, len(cell_group_ids) - 1)
            g_ids = cell_group_ids[cell_ids]
            # fast scatter via scipy coo
            coo = sp.coo_matrix((d_chunk, (g_ids, idx_chunk)), shape=(n_groups, n_vars), dtype=np.int32)
            partial += coo.toarray()
            chunks += 1
            if chunks % 5 == 0:
                print(f"[worker-{worker_id}] nnz {j:,}/{nnz_end:,} ({100*(j-nnz_start)/(nnz_end-nnz_start):.1f}% of slice)  elapsed {time.time()-t0:.0f}s", flush=True)
            i = j
    np.save(partial_path, partial)
    print(f"[worker-{worker_id}] DONE partial saved to {partial_path} ({time.time()-t0:.0f}s)", flush=True)


def main():
    t0 = time.time()
    log(f"=== greenwald FAST start {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    PARTIAL_DIR.mkdir(exist_ok=True, parents=True)

    rna_anno = load_rna_annotations()
    log(f"  rna_anno {len(rna_anno)} rows")

    with h5py.File(ATAC_H5AD, "r") as f:
        obs = f["obs"]
        cell_id = _h5col(obs, "cell_id")
        donor = _h5col(obs, "donor_id")
        var = f["var"]
        chrom = _h5col(var, "seqnames")
        start = _h5col(var, "start").astype(int)
        end = _h5col(var, "end").astype(int)
        n_cells = len(cell_id)
        n_vars = len(chrom)
        indptr = f["X"]["indptr"][:]
        nnz_total = f["X"]["data"].shape[0]
    log(f"  ATAC: {n_cells:,} cells × {n_vars:,} peaks  nnz={nnz_total:,}")

    atac_df = pd.DataFrame({"cell_id": cell_id, "donor_id": donor}).reset_index(drop=True)
    joined = atac_df.join(rna_anno[["label"]], on="cell_id")
    labels = joined["label"].fillna("unassigned").values
    donors = atac_df["donor_id"].astype(str).values

    group_key = np.array([f"{d}||{l}" for d, l in zip(donors, labels)])
    unique_groups = np.unique(group_key)
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    cell_group_ids = np.array([group_to_idx[g] for g in group_key], dtype=np.int32)
    n_groups = len(unique_groups)

    group_counts = pd.Series(group_key).value_counts()
    keep_groups = group_counts[group_counts >= 50].index.tolist()
    log(f"  {n_groups} total patient×label groups; {len(keep_groups)} with ≥50 cells")

    # Split nnz range into N_WORKERS slices
    slice_size = nnz_total // N_WORKERS
    slices = []
    for w in range(N_WORKERS):
        s = w * slice_size
        e = (w + 1) * slice_size if w < N_WORKERS - 1 else nnz_total
        slices.append((s, e))
    log(f"  spawning {N_WORKERS} workers over nnz slices: {slices}")

    args_list = []
    cgi_bytes = cell_group_ids.tobytes()
    indptr_bytes = indptr.tobytes()
    for w, (s, e) in enumerate(slices):
        partial_path = PARTIAL_DIR / f"partial_{w}.npy"
        args_list.append((w, s, e, cgi_bytes, n_groups, n_vars, indptr_bytes, partial_path))

    with mp.Pool(N_WORKERS) as pool:
        pool.map(worker, args_list)

    log(f"  all workers done, merging partials")
    group_sums = np.zeros((n_groups, n_vars), dtype=np.int64)
    for w in range(N_WORKERS):
        partial = np.load(PARTIAL_DIR / f"partial_{w}.npy")
        group_sums += partial.astype(np.int64)
    log(f"  merged.  nnz-equivalent total: {group_sums.sum():,}")

    # Overlap map
    cohort_bed = PILOT / "_greenwald_fast_src.bed"
    with open(cohort_bed, "w") as f:
        for i in range(n_vars):
            c = chrom[i]
            if not (isinstance(c, str) and c.startswith("chr")): continue
            f.write(f"{c}\t{start[i]}\t{end[i]}\t{i}\n")
    v5_df = pd.read_csv(V5_BED, sep="\t", header=None, names=["chrom","start","end","sources"])
    v5_df["v5_idx"] = np.arange(len(v5_df))
    v5_idx_bed = PILOT / "_greenwald_fast_v5idx.bed"
    v5_df[["chrom","start","end","v5_idx"]].to_csv(v5_idx_bed, sep="\t", header=False, index=False)
    out = PILOT / "_greenwald_fast_overlap.tsv"
    subprocess.run(
        f"sort -k1,1 -k2,2n {cohort_bed} > {cohort_bed}.s && "
        f"sort -k1,1 -k2,2n {v5_idx_bed} > {v5_idx_bed}.s && "
        f"bedtools intersect -a {cohort_bed}.s -b {v5_idx_bed}.s -wa -wb > {out}",
        shell=True, check=True,
    )
    rows, cols = [], []
    with open(out) as f:
        for line in f:
            p = line.rstrip().split("\t")
            try: rows.append(int(p[3])); cols.append(int(p[-1]))
            except: pass
    for p in [cohort_bed, v5_idx_bed, Path(f"{cohort_bed}.s"), Path(f"{v5_idx_bed}.s"), out]:
        if p.exists(): p.unlink()
    M = sp.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_vars, len(v5_df)))
    log(f"  overlap nnz={M.nnz:,}")

    group_sums_sp = sp.csr_matrix(group_sums)
    proj = (group_sums_sp @ M).toarray()
    log(f"  projection shape {proj.shape}")

    cols = {}
    for gi, g in enumerate(unique_groups):
        if g not in keep_groups: continue
        d, lab = g.split("||", 1)
        cols[f"greenwald_2024__{d}__{lab}"] = proj[gi, :]
    df = pd.DataFrame(cols)
    out_pq = OUT / "pseudobulk_v7_greenwald_2024_fast.parquet"
    df.to_parquet(out_pq, index=False)
    log(f"wrote {out_pq}  shape={df.shape}  total elapsed={time.time()-t0:.0f}s")

    # cleanup partials
    for w in range(N_WORKERS):
        p = PARTIAL_DIR / f"partial_{w}.npy"
        if p.exists(): p.unlink()


if __name__ == "__main__":
    sys.exit(main() or 0)
