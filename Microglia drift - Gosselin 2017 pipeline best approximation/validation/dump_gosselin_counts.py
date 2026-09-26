#!/usr/bin/env python3
"""
Extract Gosselin 2017 Table S2 (Human RNA-seq read counts) into a
plain CSV that Step 2 can consume.

Also emits sample-metadata CSVs for the N=5 paired ex-vivo → 7-day
in-vitro human microglia design used in the validation, and fixes
Excel-corrupted gene names (MARCH1→2017-03-01 etc — the Ziemann 2016
issue that persists in the published Table S2).

Usage:
  python dump_gosselin_counts.py \\
      --xlsx /path/to/aal3222_gosselin_tables2.xlsx \\
      --outdir validation/gosselin_input/
"""
from __future__ import annotations
import argparse
import datetime
from pathlib import Path

import openpyxl
import pandas as pd


def _fix_gene(v):
    """Recover gene symbol from Excel-auto-corrupted datetime."""
    if isinstance(v, (datetime.datetime, datetime.date)):
        m, d = v.month, v.day
        if m == 3:  return f"MARCH{d}"
        if m == 9:  return f"SEPT{d}"
        if m == 12: return f"DEC{d}"
    return v


def load_read_counts(xlsx_path: str) -> pd.DataFrame:
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb["Read Counts"]
    rows = list(ws.iter_rows(min_row=2, values_only=True))
    header = rows[0]
    cols = [str(h) if h else f"col{i}" for i, h in enumerate(header)]
    df = pd.DataFrame(rows[1:], columns=cols).rename(columns={cols[0]: "gene"})
    df["gene"] = df["gene"].apply(_fix_gene)
    df = df[df["gene"].notna()]
    df["gene"] = df["gene"].astype(str)
    df = df.groupby("gene", as_index=True).sum(numeric_only=True)
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0)
    df = df.T.groupby(level=0).mean().T.round().astype(int)
    return df


# ---------------------------------------------------------------------------
# The N=5 paired design used in the validation.
# Every donor here has paired ex-vivo AND 7-day in-vitro microglia
# per Gosselin Table S1 (assay code A = RNA-seq ex-vivo; A* = in-vitro).
# HMG021 has no plain-IL-34 7d sample — its 7d in-vitro samples were run
# in HumanSerum and StemPro conditions. HMG025 has 7d + ACM.
# The pipeline lets limma treat multiple in-vitro samples per donor as
# within-donor replicates in the design.
# ---------------------------------------------------------------------------
N5_AVG_ALL_IV = [
    ("HMG008", "RNA_HMG008_S011_Microglia_ExVivo",                       "ex"),
    ("HMG008", "RNA_HMG008_S014_Microglia_InVitro_7d",                   "iv"),
    ("HMG008", "RNA_HMG008_S014_Microglia_InVitro_7d_MCSF",              "iv"),
    ("HMG008", "RNA_HMG008_S014_Microglia_InVitro_7d_MCSF_TGFB",         "iv"),
    ("HMG008", "RNA_HMG008_S014_Microglia_InVitro_7d_TGFB",              "iv"),
    ("HMG015", "RNA_HMG015_S022_Microglia_ExVivo",                       "ex"),
    ("HMG015", "RNA_HMG015_S022_Microglia_InVitro_7d",                   "iv"),
    ("HMG017", "RNA_HMG017_S025_Microglia_ExVivo",                       "ex"),
    ("HMG017", "RNA_HMG017_S025_Microglia_InVitro_7d",                   "iv"),
    ("HMG017", "RNA_HMG017_S025_Microglia_InVitro_7d_TGFB",              "iv"),
    ("HMG021", "RNA_HMG021_S031_Microglia_ExVivo",                       "ex"),
    ("HMG021", "RNA_HMG021_S031_Microglia_InVitro_7d_HumanSerum",        "iv"),
    ("HMG021", "RNA_HMG021_S031_Microglia_InVitro_7d_StemPro",           "iv"),
    ("HMG021", "RNA_HMG021_S031_Microglia_InVitro_7d_StemPro_HumanSerum","iv"),
    ("HMG025", "RNA_HMG025_S036_Microglia_ExVivo",                       "ex"),
    ("HMG025", "RNA_HMG025_S036_Microglia_InVitro_7d_ACM",               "iv"),
    ("HMG025", "RNA_HMG025_S036_Microglia_InVitro_7d_ACM_HumanSerum",    "iv"),
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xlsx", required=True,
                     help="path to Gosselin 2017 Table S2 xlsx")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    out = Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    df = load_read_counts(args.xlsx)
    df.to_csv(out / "gosselin_read_counts.csv")
    print(f"wrote {out}/gosselin_read_counts.csv  ({df.shape[0]} genes x "
          f"{df.shape[1]} samples)")

    meta = pd.DataFrame(N5_AVG_ALL_IV, columns=["donor", "sample", "condition"])
    meta.to_csv(out / "gosselin_metadata_N5_avg_all_IV.csv", index=False)
    print(f"wrote {out}/gosselin_metadata_N5_avg_all_IV.csv "
          f"({len(meta)} samples, {meta.donor.nunique()} donors)")


if __name__ == "__main__":
    main()
