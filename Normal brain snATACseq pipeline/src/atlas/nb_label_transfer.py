#!/usr/bin/env python3
"""
nb_label_transfer.py
--------------------
Major-class label transfer for human brain snATAC-seq cells, scored by cosine
similarity against per-class accessibility signatures built from labeled
reference cohorts (cxg-1, cxg-2 veiga, sea-ad, GSE214979, GSE174367).

Design
======
Reference parquets at `h5ads_v2.2_nnz_rebuild/*.parquet` carry rows of
(peak_num, major_class, donor_id, n_cells_access_direct, n_cells_total_direct,
 fraction_open_direct, study_id).

- peak_num is the Design A union-peak index (1-indexed, 1..1,012,572).
- `n_cells_total_direct` is constant per (donor_id, major_class).
- Signature for class c, peak p is a donor-pooled fraction-open:
    sig(c, p) = sum_d n_cells_access(d, c, p) / sum_d n_cells_total(d, c)
  Then L2-normalized per class.

Label transfer on a new cohort takes a cells x peaks sparse matrix with peak
IDs matching the union BED. We binarize (any fragment -> 1), then score each
cell by cosine similarity to each class signature and assign argmax if the
best score exceeds a confidence floor (default 0.1), else "Unclassified".

Usage
=====
    python nb_label_transfer.py                 # validate on cxg-2 veiga
    python nb_label_transfer.py --holdout sea-ad

Validation strategy
===================
Per-cell matrices are not stored in the direct parquets, so holdout validation
uses pseudobulk rows: each (donor_id, major_class) row in the holdout parquet
becomes one test "pseudo-cell" with a vector of fraction_open values across
peaks. Each test row is scored against class signatures built from the other
reference cohorts and the argmax is compared to its true label. This is a
proxy for per-cell accuracy; actual per-cell cosine will be noisier because a
single cell is a very sparse binary draw from its class's accessibility
distribution, so pseudobulk accuracy is an optimistic upper bound.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import polars as pl
from scipy import sparse

PARQUET_DIR_DEFAULT = Path(
    "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas 10-13-2026/"
    "atlas_v2.0.0/h5ads_v2.2_nnz_rebuild"
)
UNION_BED_DEFAULT = Path(
    "/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas 10-13-2026/"
    "atlas_v1.0.0/union_peaks.bed"
)

# High-quality author-labeled references.
REFERENCE_FILES = {
    "cxg-2-veiga-adolescent": "veiga_direct_nnz1k.parquet",
    "sea-ad":                 "sea-ad_direct_nnz1k.parquet",
    "geo-GSE214979":          "GSE214979_direct_nnz1k.parquet",
    "geo-GSE174367":          "GSE174367_direct_nnz1k.parquet",
    "cxg-1-wang-kriegstein":  "cxg1_direct_nnz1k.parquet",
}

# Canonical target label set.
TARGET_CLASSES = [
    "Excitatory", "Inhibitory", "Oligodendrocyte", "OPC", "Astrocyte",
    "Microglia", "Endothelial", "Mural", "Dopaminergic", "Serotonergic",
    "Immature_neuron", "Progenitor", "Neuron_other", "Immune_other",
    "T_cell", "Meninges", "Unclassified",
]

# Classes we build signatures for = the intersect of TARGET_CLASSES and
# whatever actually appears in a labeled reference. "Unclassified" is never
# used as a signature target (we assign it from the confidence floor).
_SIGNATURE_BLOCKLIST = {"Unclassified"}


# ---------------------------------------------------------------------------
# Peak-set helpers
# ---------------------------------------------------------------------------

def load_peak_coords(union_bed: Path = UNION_BED_DEFAULT) -> pl.DataFrame:
    """Return the union peak set as a DataFrame with columns:
    peak_num (1-indexed int), chrom, start, end, peak_id ('peak_<n>')."""
    df = pl.read_csv(
        union_bed,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "start", "end", "peak_id"],
    )
    df = df.with_columns(
        peak_num=(pl.col("peak_id").str.replace("peak_", "").cast(pl.Int64))
    )
    return df.select(["peak_num", "chrom", "start", "end", "peak_id"])


def _peak_count(peaks: pl.DataFrame | None) -> int:
    if peaks is None:
        return 1012572  # Design A union length (from union_peaks.bed)
    return int(peaks["peak_num"].max())


# ---------------------------------------------------------------------------
# Signature building
# ---------------------------------------------------------------------------

def _list_reference_parquets(parquet_dir: Path) -> List[Tuple[str, Path]]:
    out = []
    for study, name in REFERENCE_FILES.items():
        p = parquet_dir / name
        if p.exists():
            out.append((study, p))
        else:
            print(f"[warn] missing reference parquet: {p}", file=sys.stderr)
    return out


def build_class_signatures(
    parquet_dir: Path | str = PARQUET_DIR_DEFAULT,
    exclude_studies: Sequence[str] = (),
    peak_total: int | None = None,
    classes: Sequence[str] | None = None,
    min_denominator_cells: int = 50,
) -> Dict[str, np.ndarray]:
    """Build per-class L2-normalized accessibility signatures.

    Parameters
    ----------
    parquet_dir : Path
        Directory with reference `*_direct_nnz1k.parquet` files.
    exclude_studies : iterable of str
        study_id values to exclude (e.g. ("cxg-2-veiga-adolescent",) for
        holdout validation).
    peak_total : int, optional
        Total number of peaks in the union set. Default inferred from
        the union_peaks.bed next to atlas_v1.0.0.
    classes : iterable of str, optional
        Which classes to build signatures for. Defaults to TARGET_CLASSES
        minus Unclassified, intersected with classes actually observed.
    min_denominator_cells : int
        Any (class, donor) with n_cells_total below this is dropped before
        aggregation. 50 is a cheap guard against noisy micro-pseudobulks.

    Returns
    -------
    dict[str, np.ndarray]
        {class_name: dense float32 vector of length peak_total, L2=1}.
    """
    parquet_dir = Path(parquet_dir)
    exclude_studies = set(exclude_studies)
    peaks = load_peak_coords() if peak_total is None else None
    P = _peak_count(peaks) if peak_total is None else int(peak_total)

    files = _list_reference_parquets(parquet_dir)
    files = [(s, p) for (s, p) in files if s not in exclude_studies]
    if not files:
        raise ValueError(
            f"No reference parquets available after excluding {exclude_studies}"
        )

    # Pool numerators across cohorts/donors.
    numerators: Dict[str, np.ndarray] = {}
    denominators: Dict[str, int] = {}

    for study, path in files:
        # Compute per-(class, donor) denominator from a lazy scan.
        donor_class_denom = (
            pl.scan_parquet(path)
              .group_by(["major_class", "donor_id"])
              .agg(pl.col("n_cells_total_direct").first().alias("n_tot"))
              .filter(pl.col("n_tot") >= min_denominator_cells)
              .collect()
        )
        kept = set(
            zip(donor_class_denom["major_class"].to_list(),
                donor_class_denom["donor_id"].to_list())
        )
        # Sum denominators per class.
        for row in donor_class_denom.iter_rows(named=True):
            if row["major_class"] in _SIGNATURE_BLOCKLIST:
                continue
            denominators[row["major_class"]] = (
                denominators.get(row["major_class"], 0) + int(row["n_tot"])
            )

        # Sum numerators per (class, peak).
        agg = (
            pl.scan_parquet(path)
              .select(["peak_num", "major_class", "donor_id",
                       "n_cells_access_direct"])
              .filter(~pl.col("major_class").is_in(list(_SIGNATURE_BLOCKLIST)))
              .collect(engine="streaming")
        )
        # Filter to kept (class, donor) pairs via a semi-join.
        kept_df = pl.DataFrame({
            "major_class": [k[0] for k in kept],
            "donor_id":    [k[1] for k in kept],
        })
        agg = agg.join(kept_df, on=["major_class", "donor_id"], how="inner")
        if agg.height == 0:
            continue
        agg = agg.group_by(["major_class", "peak_num"]).agg(
            pl.col("n_cells_access_direct").sum().alias("num")
        )

        for cls, sub in agg.partition_by("major_class", as_dict=True).items():
            cls = cls[0] if isinstance(cls, tuple) else cls
            if cls not in numerators:
                numerators[cls] = np.zeros(P, dtype=np.float64)
            peak_idx = sub["peak_num"].to_numpy() - 1  # 1-indexed -> 0-indexed
            vals = sub["num"].to_numpy()
            # Guard: ignore peak_num outside Design A (shouldn't happen).
            mask = (peak_idx >= 0) & (peak_idx < P)
            np.add.at(numerators[cls], peak_idx[mask], vals[mask])
        print(f"[sig] + {study}: classes={sorted(agg['major_class'].unique().to_list())}",
              file=sys.stderr)

    want = set(classes) if classes else (set(TARGET_CLASSES) - _SIGNATURE_BLOCKLIST)
    signatures: Dict[str, np.ndarray] = {}
    for cls, num in numerators.items():
        if cls not in want:
            continue
        denom = denominators.get(cls, 0)
        if denom == 0:
            continue
        frac = (num / float(denom)).astype(np.float32)
        norm = float(np.linalg.norm(frac))
        if norm == 0.0:
            continue
        signatures[cls] = frac / norm

    print(
        f"[sig] built signatures for {len(signatures)} classes: "
        f"{sorted(signatures)}", file=sys.stderr,
    )
    return signatures


# ---------------------------------------------------------------------------
# Label transfer
# ---------------------------------------------------------------------------

def _signatures_to_matrix(
    signatures: Dict[str, np.ndarray]
) -> Tuple[List[str], np.ndarray]:
    classes = sorted(signatures.keys())
    P = next(iter(signatures.values())).shape[0]
    M = np.empty((len(classes), P), dtype=np.float32)
    for i, c in enumerate(classes):
        v = signatures[c]
        if v.shape[0] != P:
            raise ValueError(f"Signature '{c}' has length {v.shape[0]} != {P}")
        M[i] = v
    return classes, M


def label_transfer(
    X_sparse,
    peak_ids: Iterable[int] | np.ndarray,
    signatures: Dict[str, np.ndarray],
    *,
    conf_threshold: float = 0.1,
    binarize: bool = True,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transfer class labels onto a cells x peaks sparse matrix.

    Parameters
    ----------
    X_sparse : scipy.sparse matrix (csr/coo/csc)
        Cells x peaks fragment or binary count matrix. Row-order of cells is
        preserved; column order must match `peak_ids`.
    peak_ids : array-like of int
        1-indexed Design A `peak_num` for each column of `X_sparse`. The query
        matrix may be a subset of the full union — unknown peaks are ignored.
    signatures : dict[str, np.ndarray]
        Output of `build_class_signatures`. Each array is a full-length
        (peak_total,) L2-normalized vector.
    conf_threshold : float
        Cells whose best cosine similarity is <= this are called
        "Unclassified".
    binarize : bool
        If True (default) the input is binarized (any fragment -> 1) before
        scoring. ATAC accessibility is naturally binary at the single-cell
        level so this is the right default.

    Returns
    -------
    labels : np.ndarray of str, shape (n_cells,)
    scores : np.ndarray of float32, shape (n_cells,)
        Best cosine similarity.
    full_scores : np.ndarray of float32, shape (n_cells, n_classes)
        Per-class cosine similarity (columns in sorted(class name) order).
    """
    classes, M = _signatures_to_matrix(signatures)   # (C, P)
    P_full = M.shape[1]

    peak_ids = np.asarray(peak_ids, dtype=np.int64)
    if peak_ids.shape[0] != X_sparse.shape[1]:
        raise ValueError(
            f"peak_ids length ({peak_ids.shape[0]}) must match "
            f"X_sparse.shape[1] ({X_sparse.shape[1]})"
        )

    # Reduce the signature matrix to the columns present in the query.
    col_idx = peak_ids - 1  # 1-indexed -> 0-indexed
    valid = (col_idx >= 0) & (col_idx < P_full)
    if not valid.all():
        X_sparse = X_sparse[:, valid]
        col_idx = col_idx[valid]
    M_sub = M[:, col_idx]  # (C, P_query)

    X = X_sparse.tocsr().astype(np.float32, copy=False)
    if binarize:
        # scipy>=1.11 has sign(); fall back to a simple inline if not.
        X = X.copy()
        X.data = np.ones_like(X.data, dtype=np.float32)

    # Score in chunks to control memory on big cohorts.
    n_cells = X.shape[0]
    chunk = 20_000
    full_scores = np.zeros((n_cells, len(classes)), dtype=np.float32)
    row_norm = np.sqrt(X.multiply(X).sum(axis=1)).A1.astype(np.float32)  # (n,)
    for start in range(0, n_cells, chunk):
        stop = min(start + chunk, n_cells)
        # dot = X_chunk @ M_sub.T -> (chunk, C)
        dot = X[start:stop].dot(M_sub.T)
        # Signatures are L2=1, so cosine = dot / ||x||. Guard zero-nnz cells.
        rn = row_norm[start:stop]
        with np.errstate(divide="ignore", invalid="ignore"):
            cos = np.where(rn[:, None] > 0, dot / rn[:, None], 0.0)
        full_scores[start:stop] = cos.astype(np.float32)

    best_idx = np.argmax(full_scores, axis=1)
    best_score = full_scores[np.arange(n_cells), best_idx]
    labels = np.array(classes)[best_idx]
    labels = np.where(best_score > conf_threshold, labels, "Unclassified")
    return labels, best_score.astype(np.float32), full_scores


# ---------------------------------------------------------------------------
# Holdout validation (pseudobulk)
# ---------------------------------------------------------------------------

def _pseudobulk_test_matrix(
    parquet_path: Path,
    peak_total: int,
    min_denominator_cells: int = 50,
) -> Tuple[sparse.csr_matrix, np.ndarray, List[str], List[str]]:
    """Pivot a holdout direct parquet into one row per (donor, class).

    Each row is a dense fraction_open_direct vector across all union peaks,
    returned as a sparse CSR matrix (donors*classes) x peaks for memory.

    Returns
    -------
    X : csr_matrix of float32, (n_rows, peak_total)
    peak_ids : np.ndarray[int], 1..peak_total
    true_labels : list[str], length n_rows
    donor_ids : list[str], length n_rows
    """
    print(f"[val] loading {parquet_path.name}", file=sys.stderr)
    df = pl.scan_parquet(parquet_path).collect(engine="streaming")
    # Drop micro-pseudobulks.
    denom = (
        df.group_by(["major_class", "donor_id"])
          .agg(pl.col("n_cells_total_direct").first().alias("n_tot"))
          .filter(pl.col("n_tot") >= min_denominator_cells)
    )
    df = df.join(denom.select(["major_class", "donor_id"]),
                 on=["major_class", "donor_id"], how="inner")

    # Enumerate rows.
    row_meta = (
        df.select(["major_class", "donor_id"]).unique().sort(
            ["major_class", "donor_id"]
        )
    )
    row_meta = row_meta.with_row_index("row_idx")
    df = df.join(row_meta, on=["major_class", "donor_id"], how="left")

    row_idx = df["row_idx"].to_numpy()
    col_idx = df["peak_num"].to_numpy() - 1
    vals = df["fraction_open_direct"].to_numpy().astype(np.float32)
    mask = (col_idx >= 0) & (col_idx < peak_total)
    X = sparse.coo_matrix(
        (vals[mask], (row_idx[mask], col_idx[mask])),
        shape=(row_meta.height, peak_total),
    ).tocsr()
    peak_ids = np.arange(1, peak_total + 1, dtype=np.int64)
    true_labels = row_meta["major_class"].to_list()
    donor_ids = row_meta["donor_id"].to_list()
    return X, peak_ids, true_labels, donor_ids


def _confusion(true_labels: Sequence[str], pred_labels: Sequence[str],
               classes: Sequence[str]) -> np.ndarray:
    idx = {c: i for i, c in enumerate(classes)}
    M = np.zeros((len(classes), len(classes)), dtype=np.int64)
    for t, p in zip(true_labels, pred_labels):
        if t in idx and p in idx:
            M[idx[t], idx[p]] += 1
    return M


def _print_confusion(cm: np.ndarray, classes: Sequence[str]) -> None:
    w = max(14, max(len(c) for c in classes) + 1)
    header = "true \\ pred".ljust(w) + "".join(c[:w-1].ljust(w) for c in classes) + "total"
    print(header)
    for i, c in enumerate(classes):
        row = cm[i]
        tot = int(row.sum())
        line = c.ljust(w) + "".join(str(int(x)).ljust(w) for x in row) + str(tot)
        print(line)


def validate_holdout(
    parquet_dir: Path = PARQUET_DIR_DEFAULT,
    holdout_study: str = "cxg-2-veiga-adolescent",
    conf_threshold: float = 0.1,
) -> Dict[str, object]:
    """Build signatures excluding `holdout_study`, score the holdout parquet
    pseudobulks, report accuracy and confusion matrix."""
    peaks = load_peak_coords()
    P = int(peaks["peak_num"].max())

    t0 = time.time()
    signatures = build_class_signatures(
        parquet_dir=parquet_dir,
        exclude_studies=[holdout_study],
        peak_total=P,
    )
    t_sig = time.time() - t0
    print(f"[val] signatures built in {t_sig:.1f}s", file=sys.stderr)

    # Locate the holdout parquet.
    holdout_file = parquet_dir / REFERENCE_FILES[holdout_study]
    X, peak_ids, true_labels, donors = _pseudobulk_test_matrix(holdout_file, P)
    print(f"[val] test rows: {X.shape[0]} across "
          f"{len(set(true_labels))} classes", file=sys.stderr)

    t1 = time.time()
    pred, scores, full = label_transfer(
        X, peak_ids, signatures,
        conf_threshold=conf_threshold,
        binarize=False,  # test rows are already fraction_open floats
    )
    t_pred = time.time() - t1
    print(f"[val] scored {X.shape[0]} rows in {t_pred:.1f}s", file=sys.stderr)

    classes = sorted(set(signatures.keys()) | set(true_labels) | {"Unclassified"})
    cm = _confusion(true_labels, list(pred), classes)
    overall = (np.array(true_labels) == np.array(pred)).mean()

    print("\n=== Holdout: {} ===".format(holdout_study))
    print(f"Signatures built from {sum(1 for s in REFERENCE_FILES if s != holdout_study)} "
          f"reference cohorts, excluding '{holdout_study}'.")
    print(f"Overall pseudobulk accuracy: {overall:.3f}")
    print(f"Confidence threshold: {conf_threshold}")

    print("\nPer-class accuracy (true -> argmax predicted):")
    for i, c in enumerate(classes):
        row = cm[i]
        tot = int(row.sum())
        if tot == 0:
            continue
        correct = int(cm[i, i]) if c in sorted(set(signatures.keys())) else 0
        # If true class never had a signature, treat the "correct" column as 0.
        if c not in signatures:
            correct_mark = "  (no reference signature)"
        else:
            correct_mark = ""
        acc = correct / tot if tot else float("nan")
        print(f"  {c:20s}  n={tot:3d}  acc={acc:.3f}{correct_mark}")

    print("\nConfusion matrix (rows=true class, cols=predicted class):")
    _print_confusion(cm, classes)

    # Show the predicted-class columns that appear
    print("\nMean best cosine by true class:")
    labels_arr = np.array(true_labels)
    for c in sorted(set(true_labels)):
        m = scores[labels_arr == c]
        print(f"  {c:20s}  mean={m.mean():.3f}  min={m.min():.3f}  "
              f"max={m.max():.3f}  n={len(m)}")

    return {
        "signatures": signatures,
        "true_labels": true_labels,
        "pred_labels": list(pred),
        "scores": scores,
        "donors": donors,
        "confusion": cm,
        "classes": classes,
        "accuracy": float(overall),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parquet-dir", type=Path, default=PARQUET_DIR_DEFAULT)
    ap.add_argument("--holdout", default="cxg-2-veiga-adolescent",
                    choices=list(REFERENCE_FILES.keys()))
    ap.add_argument("--conf-threshold", type=float, default=0.1)
    return ap.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    validate_holdout(
        parquet_dir=args.parquet_dir,
        holdout_study=args.holdout,
        conf_threshold=args.conf_threshold,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
