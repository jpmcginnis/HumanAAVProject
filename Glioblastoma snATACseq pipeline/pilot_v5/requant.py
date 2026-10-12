"""v7 re-quantification — patient × pool-rule-label pseudobulks.

Change from v6: aggregates per (patient_id, pool_label) instead of
(cohort, pool_label). This makes patient the primary unit of replication.

Each column of the resulting matrix is named `<cohort>__<patient_id>__<label>`
so the cohort can be recovered for batch sanity.

Filters out patient × label groups with fewer than MIN_CELLS cells.
"""
import os
import sys
import glob
import time
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import h5py
import snapatac2 as snap
import anndata as ad_mod

ROOT = Path("/data/projects/atacseq")
PILOT = ROOT / "pilot_v5"
V5_BED = str(PILOT / "v5_union_reference.bed")
OUT = PILOT / "requant"
OUT.mkdir(exist_ok=True, parents=True)
PROCESSED = ROOT / "processed"

LOG = PILOT / "requant_v7.log"
MIN_CELLS = 50  # minimum cells per (patient × label) group to pseudobulk


def log(m):
    s = f"[requant-v7] {m}"
    print(s, flush=True)
    with open(LOG, "a") as f:
        f.write(s + "\n")


COHORTS = {
    "sundaram_2024":  {"type": "fragments",
                       "glob": str(PROCESSED / "sundaram_gbm" / "scATAC_GBMx_*.h5ad"),
                       "labeled": str(PROCESSED / "sundaram_gbm" / "catlas_quantified.h5ad")},
    "sojka_2025":     {"type": "fragments",
                       "glob": str(PROCESSED / "gse276177_khan_astro" / "GSM*_fragments.tsv.gz.h5ad"),
                       "labeled": str(PROCESSED / "gse276177_khan_astro" / "catlas_quantified.h5ad")},
    "zhao_2024":      {"type": "matrix",
                       "h5ad": str(PROCESSED / "gbm_tme_atlas_hra004942" / "catlas_quantified.h5ad")},
    "mathewson_2021": {"type": "matrix",
                       "h5ad": str(PROCESSED / "mathewson_lupien" / "catlas_quantified.h5ad")},
    "guo_2020":       {"type": "matrix",
                       "h5ad": str(PROCESSED / "gse138794_guo" / "catlas_quantified.h5ad")},
    "guilhamon_2021": {"type": "matrix",
                       "h5ad": str(PROCESSED / "guilhamon" / "catlas_quantified.h5ad")},
    "nagaraja_2021":  {"type": "matrix",
                       "h5ad": str(PROCESSED / "gse165037" / "catlas_quantified.h5ad")},
}


def pool_label(ct, mc):
    try:
        if int(mc) == 1:
            return "malignant"
    except (TypeError, ValueError):
        pass
    if ct == "malignant_unresolved":
        return "malignant"
    return ct


def read_v5_bed():
    df = pd.read_csv(V5_BED, sep="\t", header=None, names=["chrom", "start", "end", "sources"])
    df["peak_id"] = df["chrom"] + ":" + df["start"].astype(str) + "-" + df["end"].astype(str)
    df["v5_idx"] = np.arange(len(df))
    return df


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


def load_cohort_annotations(labeled_path):
    """Return dict barcode -> (pool_label, patient_id)."""
    with h5py.File(labeled_path, "r") as f:
        obs = f["obs"]
        idx_key = obs.attrs.get("_index", "index")
        if isinstance(idx_key, bytes): idx_key = idx_key.decode()
        barcodes = _h5col(obs, idx_key)
        cts = _h5col(obs, "cell_type")
        mcs = _h5col(obs, "malignant_cnv") if "malignant_cnv" in obs.keys() else np.zeros(len(barcodes), dtype=int)
        pids = _h5col(obs, "patient_id") if "patient_id" in obs.keys() else np.full(len(barcodes), "unknown", dtype=object)
    labels = [pool_label(ct, mc) for ct, mc in zip(cts, mcs)]
    return {b: (lab, str(pid)) for b, lab, pid in zip(barcodes, labels, pids)}


def detect_sample_suffix(persample_barcodes, cohort_barcode_set):
    sample = persample_barcodes[:200]
    best_suffix, best_match = "", 0
    for i in range(25):
        suffix = f"-{i}"
        matches = sum(1 for b in sample if (b + suffix) in cohort_barcode_set)
        if matches > best_match:
            best_match = matches; best_suffix = suffix
    plain = sum(1 for b in sample if b in cohort_barcode_set)
    if plain > best_match:
        best_match = plain; best_suffix = ""
    return best_suffix, best_match / len(sample)


def read_snap_peak_matrix_X(peak_h5ad_path):
    a = snap.read(peak_h5ad_path, backed="r")
    X = a.X[:]
    if sp.issparse(X):
        return X.tocsr(), a.n_obs, a.n_vars, list(a.obs_names)
    return sp.csr_matrix(X), a.n_obs, a.n_vars, list(a.obs_names)


def aggregate_patient_label(X, patient_ids, labels, min_cells=MIN_CELLS):
    """Returns dict (patient_id, label) -> col_sum vector."""
    out = {}
    for pat in np.unique(patient_ids):
        pmask = patient_ids == pat
        for lab in np.unique(labels[pmask]):
            if lab in ("", "nan"): continue
            mask = pmask & (labels == lab)
            n = int(mask.sum())
            if n < min_cells: continue
            key = (pat, lab)
            vec = np.asarray(X[mask, :].sum(axis=0)).flatten().astype(np.int64)
            if key in out:
                out[key] += vec
            else:
                out[key] = vec
    return out


def requant_fragments(slug, cfg, v5_df):
    log(f"=== {slug} (fragments, v7 patient × pool-label) ===")
    anno = load_cohort_annotations(cfg["labeled"])
    log(f"  loaded {len(anno)} barcode→(label, patient)")
    cohort_bc_set = set(anno.keys())
    samples = sorted(glob.glob(cfg["glob"]))
    samples = [p for p in samples if "_peakmat" not in p and "catlas_quantified" not in Path(p).name]
    log(f"  {len(samples)} per-sample fragment h5ads")

    patient_label_sums = {}  # (patient, label) -> cumulative count vector
    for i, sh in enumerate(samples):
        name = Path(sh).stem
        t0 = time.time()
        tmp_pm = f"/data/projects/atacseq/pilot_v5/_requant_v7_{slug}_{name}_pm.h5ad"
        try:
            a = snap.read(sh, backed="r")
            if a.n_obs < 100:
                log(f"  [{i+1}/{len(samples)}] skip {name}")
                continue
            log(f"  [{i+1}/{len(samples)}] {name} n_obs={a.n_obs} — make_peak_matrix")
            pm = snap.pp.make_peak_matrix(a, peak_file=V5_BED, file=tmp_pm, chunk_size=1000)
            X, n_obs, n_vars, obs_names = read_snap_peak_matrix_X(tmp_pm)
            suffix, match_rate = detect_sample_suffix(obs_names, cohort_bc_set)
            log(f"    suffix='{suffix}' match_rate={match_rate:.2%}")
            labels = []
            pids = []
            for b in obs_names:
                lab, pid = anno.get(b + suffix, ("unassigned", "unknown"))
                labels.append(lab); pids.append(pid)
            labels = np.array(labels); pids = np.array(pids)
            sample_agg = aggregate_patient_label(X, pids, labels)
            for key, vec in sample_agg.items():
                if key in patient_label_sums:
                    patient_label_sums[key] += vec
                else:
                    patient_label_sums[key] = vec
            log(f"    {len(sample_agg)} patient×label groups ≥{MIN_CELLS} cells  {time.time()-t0:.1f}s")
        except BaseException as e:
            log(f"  [{i+1}/{len(samples)}] ERROR {name}: {type(e).__name__}: {e}")
        finally:
            if os.path.exists(tmp_pm):
                os.remove(tmp_pm)

    if not patient_label_sums:
        return pd.DataFrame()
    df = pd.DataFrame({f"{slug}__{pat}__{lab}": v for (pat, lab), v in patient_label_sums.items()})
    log(f"  final: {df.shape[1]} patient×label pseudobulk cols")
    return df


def _parse_peak_coords(a_var):
    cols = {c.lower(): c for c in a_var.columns}
    chrom_k = cols.get("chrom") or cols.get("chr") or cols.get("seqnames")
    start_k = cols.get("start")
    end_k = cols.get("end")
    if chrom_k and start_k and end_k:
        return (a_var[chrom_k].astype(str).values,
                a_var[start_k].astype(int).values,
                a_var[end_k].astype(int).values)
    ids = a_var.index.astype(str)
    chrom, start, end = [], [], []
    for s in ids:
        tok = s.replace(":", "\t").replace("-", "\t").split("\t")
        chrom.append(tok[0] if tok else "")
        try: start.append(int(tok[1])); end.append(int(tok[2]))
        except: start.append(0); end.append(0)
    return np.array(chrom), np.array(start), np.array(end)


def build_overlap_map(src_coords, v5_df, slug):
    src_bed = PILOT / f"_requant_v7_{slug}_src.bed"
    with open(src_bed, "w") as f:
        for i in range(len(src_coords[0])):
            c = src_coords[0][i]
            if not (isinstance(c, str) and c.startswith("chr")): continue
            f.write(f"{c}\t{src_coords[1][i]}\t{src_coords[2][i]}\t{i}\n")
    v5_bed_idx = PILOT / f"_requant_v7_{slug}_v5idx.bed"
    v5_df[["chrom","start","end","v5_idx"]].to_csv(v5_bed_idx, sep="\t", header=False, index=False)
    out = PILOT / f"_requant_v7_{slug}_overlap.tsv"
    subprocess.run(
        f"sort -k1,1 -k2,2n {src_bed} > {src_bed}.s && "
        f"sort -k1,1 -k2,2n {v5_bed_idx} > {v5_bed_idx}.s && "
        f"bedtools intersect -a {src_bed}.s -b {v5_bed_idx}.s -wa -wb > {out}",
        shell=True, check=True,
    )
    rows, cols = [], []
    with open(out) as f:
        for line in f:
            p = line.rstrip().split("\t")
            try: rows.append(int(p[3])); cols.append(int(p[-1]))
            except: pass
    n_src = len(src_coords[0])
    M = sp.csr_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(n_src, len(v5_df)))
    for p in [src_bed, v5_bed_idx, Path(f"{src_bed}.s"), Path(f"{v5_bed_idx}.s"), out]:
        if p.exists(): p.unlink()
    log(f"  overlap nnz={M.nnz:,}")
    return M


def requant_matrix(slug, cfg, v5_df):
    log(f"=== {slug} (matrix, v7 patient × pool-label) ===")
    a = ad_mod.read_h5ad(cfg["h5ad"])
    log(f"  n_obs={a.n_obs:,} n_vars={a.n_vars:,}")
    chrom, start, end = _parse_peak_coords(a.var)
    M = build_overlap_map((chrom, start, end), v5_df, slug)
    X = a.X.tocsr() if sp.issparse(a.X) else sp.csr_matrix(a.X)
    log(f"  projecting {X.shape} @ {M.shape}")
    Xv5 = X @ M
    cts_raw = a.obs.get("cell_type", pd.Series(["unassigned"] * a.n_obs)).astype(str).values
    mcs = a.obs.get("malignant_cnv", pd.Series([0] * a.n_obs)).astype(int).values
    labels = np.array([pool_label(ct, mc) for ct, mc in zip(cts_raw, mcs)])
    pids = a.obs.get("patient_id", pd.Series(["unknown"] * a.n_obs)).astype(str).values
    agg = aggregate_patient_label(Xv5, pids, labels)
    log(f"  {len(agg)} patient×label groups ≥{MIN_CELLS} cells")
    df = pd.DataFrame({f"{slug}__{pat}__{lab}": v for (pat, lab), v in agg.items()})
    return df


def main():
    log(f"=== v7 requant start {time.strftime('%Y-%m-%d %H:%M:%S')} MIN_CELLS={MIN_CELLS} ===")
    v5_df = read_v5_bed()
    log(f"v5 union: {len(v5_df):,} peaks")

    pbs = []
    for slug, cfg in COHORTS.items():
        out_pq = OUT / f"pseudobulk_v7_{slug}.parquet"
        if out_pq.exists() and out_pq.stat().st_size > 1000:
            log(f"SKIP {slug}")
            pbs.append(pd.read_parquet(out_pq)); continue
        try:
            t0 = time.time()
            if cfg["type"] == "fragments":
                pb = requant_fragments(slug, cfg, v5_df)
            elif cfg["type"] == "matrix":
                pb = requant_matrix(slug, cfg, v5_df)
            else:
                continue
            if len(pb.columns):
                pb.to_parquet(out_pq, index=False)
                log(f"  wrote {out_pq}  cols={len(pb.columns)}  {time.time()-t0:.0f}s")
                pbs.append(pb)
        except Exception as e:
            import traceback
            log(f"  FAIL {slug}: {e}\n{traceback.format_exc()}")

    if pbs:
        peak_index = v5_df[["chrom","start","end","peak_id"]].copy()
        combined = pd.concat(pbs, axis=1)
        combined.insert(0, "peak_id", peak_index["peak_id"].values)
        combined.to_parquet(OUT / "enhancer_candidate_matrix_v7.parquet", index=False)
        log(f"wrote combined matrix {combined.shape}")

    log(f"=== v7 requant end {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    sys.exit(main() or 0)
