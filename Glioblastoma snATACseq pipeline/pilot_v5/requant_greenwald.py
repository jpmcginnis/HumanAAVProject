"""Greenwald integration — label transfer from snRNA + sequential projection.

Steps:
  1. Read snRNA obs: (cell_id, annotation_coarse, annotation_granular,
     Cell_status, donor_id). Build dict cell_id -> (label, patient).
  2. Read ATAC filtered_peaks obs: (cell_id, Sample, donor_id). The ATAC
     cell_id and RNA cell_id share the same prefix+barcode in Multiome.
  3. Join: for each ATAC cell, lookup its (label, patient) from snRNA.
     Apply pool rule: Cell_status == "malignant" OR annotation_coarse contains
     "malignant" OR cnv_corr > threshold → "malignant"; else annotation_coarse.
  4. Build sparse group-indicator matrix G (n_cells × n_groups).
  5. Sequential pass over X.data + X.indices + indptr:
     for each chunk of nnz, compute (cell_id, var_id, count) tuples, scatter
     into group_sums[group_id_of_cell, var_id].
  6. Project group_sums @ M (where M is cohort_peaks × v5_peaks overlap).
  7. Write as patient × cell_type pseudobulk parquet.
"""
import sys
import time
import subprocess
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
LOG = PILOT / "requant_greenwald.log"

CHUNK_NNZ = 50_000_000


def log(m):
    s = f"[greenwald-v1] {m}"
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
    """Return DataFrame indexed by cell_id with (label, patient_id)."""
    log(f"opening {RNA_H5AD}")
    with h5py.File(RNA_H5AD, "r", rdcc_nbytes=1 << 30) as f:
        obs = f["obs"]
        idx_key = obs.attrs.get("_index", "index")
        if isinstance(idx_key, bytes): idx_key = idx_key.decode()
        cell_ids = _h5col(obs, idx_key) if idx_key in obs.keys() else _h5col(obs, "cell_id")
        donor = _h5col(obs, "donor_id") if "donor_id" in obs.keys() else np.array(["unknown"] * len(cell_ids))
        # Priority: annotation_granular > annotation_coarse > celltypist > scPoli
        if "annotation_granular" in obs.keys():
            lab = _h5col(obs, "annotation_granular")
            which = "annotation_granular"
        elif "annotation_coarse" in obs.keys():
            lab = _h5col(obs, "annotation_coarse")
            which = "annotation_coarse"
        elif "celltypist" in obs.keys():
            lab = _h5col(obs, "celltypist")
            which = "celltypist"
        else:
            lab = np.array(["unknown"] * len(cell_ids))
            which = "fallback"
        status = _h5col(obs, "Cell_status") if "Cell_status" in obs.keys() else None
    log(f"  used '{which}' as label; {len(cell_ids):,} cells; {len(np.unique(lab))} unique labels")
    uniq_labels = pd.Series(lab).value_counts()
    log(f"  top labels: {uniq_labels.head(15).to_dict()}")
    if status is not None:
        log(f"  Cell_status unique: {pd.Series(status).value_counts().to_dict()}")

    # Apply pool rule — a cell is malignant if Cell_status says so OR label contains 'malignant'/'tumor'
    labels = []
    for i in range(len(cell_ids)):
        is_mal = False
        if status is not None and str(status[i]).lower() in ("malignant", "tumor", "tumour"):
            is_mal = True
        if "malignant" in str(lab[i]).lower() or "tumor" in str(lab[i]).lower():
            is_mal = True
        if is_mal:
            # Prefer the granular label if it carries a Neftel state, else 'malignant'
            l = str(lab[i])
            # Keep granular (e.g. "MES-like", "NPC-like") but prefix with malignant_
            if any(s in l for s in ["MES", "NPC", "AC", "OPC", "Hypox", "Mes", "Prolif", "G1S", "G2M"]):
                labels.append(f"malignant_{l}")
            else:
                labels.append("malignant")
        else:
            labels.append(str(lab[i]))
    log(f"  after pool rule — top labels: {pd.Series(labels).value_counts().head(15).to_dict()}")

    df = pd.DataFrame({"cell_id": cell_ids, "donor_id": donor, "label": labels}).set_index("cell_id")
    return df


def load_atac_annotations_and_x_meta():
    """Return ATAC per-cell DataFrame with (cell_id, donor, Sample) + handle to X.
    Also returns n_vars, chrom/start/end for peak coords."""
    log(f"opening ATAC {ATAC_H5AD}")
    with h5py.File(ATAC_H5AD, "r", rdcc_nbytes=1 << 30) as f:
        obs = f["obs"]
        cell_id = _h5col(obs, "cell_id")
        donor = _h5col(obs, "donor_id")
        # Peak coords for overlap map
        var = f["var"]
        chrom = _h5col(var, "seqnames")
        start = _h5col(var, "start").astype(int)
        end = _h5col(var, "end").astype(int)
        shape = tuple(f["X"].attrs.get("shape"))
    log(f"  ATAC {shape[0]:,} cells × {shape[1]:,} peaks")
    atac_df = pd.DataFrame({"cell_id": cell_id, "donor_id": donor}).reset_index(drop=True)
    return atac_df, (chrom, start, end), shape


def build_overlap_map(src_coords, v5_df):
    src_bed = PILOT / "_greenwald_v1_src.bed"
    with open(src_bed, "w") as f:
        for i in range(len(src_coords[0])):
            c = src_coords[0][i]
            if not (isinstance(c, str) and c.startswith("chr")): continue
            f.write(f"{c}\t{src_coords[1][i]}\t{src_coords[2][i]}\t{i}\n")
    v5_idx_bed = PILOT / "_greenwald_v1_v5idx.bed"
    v5_df[["chrom","start","end","v5_idx"]].to_csv(v5_idx_bed, sep="\t", header=False, index=False)
    out = PILOT / "_greenwald_v1_overlap.tsv"
    subprocess.run(
        f"sort -k1,1 -k2,2n {src_bed} > {src_bed}.s && "
        f"sort -k1,1 -k2,2n {v5_idx_bed} > {v5_idx_bed}.s && "
        f"bedtools intersect -a {src_bed}.s -b {v5_idx_bed}.s -wa -wb > {out}",
        shell=True, check=True,
    )
    rows, cols = [], []
    with open(out) as f:
        for line in f:
            p = line.rstrip().split("\t")
            try: rows.append(int(p[3])); cols.append(int(p[-1]))
            except: pass
    for p in [src_bed, v5_idx_bed, Path(f"{src_bed}.s"), Path(f"{v5_idx_bed}.s"), out]:
        if p.exists(): p.unlink()
    n_src = len(src_coords[0])
    M = sp.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_src, len(v5_df)))
    log(f"  overlap: {M.nnz:,} nnz out of {n_src} src × {len(v5_df)} v5")
    return M


def main():
    t0 = time.time()
    log(f"=== greenwald v1 start {time.strftime('%Y-%m-%d %H:%M:%S')} ===")

    rna_anno = load_rna_annotations()

    atac_df, (chrom, start, end), (n_cells, n_vars) = load_atac_annotations_and_x_meta()
    # Join labels
    joined = atac_df.join(rna_anno[["label"]], on="cell_id")
    n_matched = joined["label"].notna().sum()
    log(f"  joined ATAC → RNA labels: {n_matched:,}/{n_cells:,} ({100*n_matched/n_cells:.1f}%) matched")
    labels = joined["label"].fillna("unassigned").values
    donors = atac_df["donor_id"].astype(str).values

    # Build (patient × label) groups
    group_key = np.array([f"{d}||{l}" for d, l in zip(donors, labels)])
    unique_groups = np.unique(group_key)
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}
    cell_group_ids = np.array([group_to_idx[g] for g in group_key], dtype=np.int32)
    n_groups = len(unique_groups)

    # Filter groups with < 50 cells
    group_counts = pd.Series(group_key).value_counts()
    keep_groups = group_counts[group_counts >= 50].index.tolist()
    log(f"  {n_groups} total patient×label groups; {len(keep_groups)} with ≥50 cells")
    keep_mask = pd.Series(group_key).isin(keep_groups).values
    log(f"  {keep_mask.sum():,}/{n_cells:,} cells land in a kept group")

    # Preallocate accumulator
    group_sums = np.zeros((n_groups, n_vars), dtype=np.int64)

    # Load v5 union
    v5_df = pd.read_csv(V5_BED, sep="\t", header=None, names=["chrom","start","end","sources"])
    v5_df["v5_idx"] = np.arange(len(v5_df))

    # Sequential pass over X
    log(f"starting sequential CSR scan — {n_cells:,} cells × {n_vars:,} peaks")
    with h5py.File(ATAC_H5AD, "r", rdcc_nbytes=4*(1 << 30), rdcc_nslots=100003) as f:
        data_ds = f["X"]["data"]
        indices_ds = f["X"]["indices"]
        indptr = f["X"]["indptr"][:]
        nnz_total = data_ds.shape[0]
        log(f"  total nnz={nnz_total:,}")

        i = 0
        chunks = 0
        while i < nnz_total:
            j = min(i + CHUNK_NNZ, nnz_total)
            tt = time.time()
            d_chunk = data_ds[i:j]
            idx_chunk = indices_ds[i:j]
            entry_pos = np.arange(i, j, dtype=np.int64)
            cell_ids = np.searchsorted(indptr, entry_pos, side="right") - 1
            cell_ids = np.clip(cell_ids, 0, n_cells - 1)
            g_ids = cell_group_ids[cell_ids]
            np.add.at(group_sums, (g_ids, idx_chunk), d_chunk)
            chunks += 1
            if chunks % 20 == 0 or j == nnz_total:
                pct = 100 * j / nnz_total
                log(f"  nnz {j:,}/{nnz_total:,} ({pct:.1f}%)  chunk {time.time()-tt:.1f}s")
            i = j
    log(f"  CSR scan done in {time.time()-t0:.0f}s")

    # Build overlap map
    M = build_overlap_map((chrom, start, end), v5_df)

    # Project each group to v5
    # group_sums is (n_groups, n_vars); M is (n_vars, n_v5)
    log("projecting group_sums @ M")
    group_sums_sp = sp.csr_matrix(group_sums)
    proj = group_sums_sp @ M  # (n_groups, n_v5)
    proj_dense = np.asarray(proj.todense())
    log(f"  projection shape {proj_dense.shape}")

    # Write as DataFrame with columns = greenwald_2024__<patient>__<label>
    cols = {}
    for gi, g in enumerate(unique_groups):
        if g not in keep_groups:
            continue
        donor, label = g.split("||", 1)
        cols[f"greenwald_2024__{donor}__{label}"] = proj_dense[gi, :]
    df = pd.DataFrame(cols)
    out_pq = OUT / "pseudobulk_v7_greenwald_2024.parquet"
    df.to_parquet(out_pq, index=False)
    log(f"wrote {out_pq}  shape={df.shape}  elapsed={time.time()-t0:.0f}s")


if __name__ == "__main__":
    sys.exit(main() or 0)
