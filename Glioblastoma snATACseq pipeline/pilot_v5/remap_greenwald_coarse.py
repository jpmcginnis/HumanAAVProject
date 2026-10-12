"""Remap Greenwald's fine-grained labels to a coarse schema compatible with
the other 7 cohorts, PLUS keep Neftel states (MES / AC / OPC / NPC) under
`malignant_<STATE>` buckets so the atlas can score per-state malignant separately
from pan_malignant.

The fine-grained pseudobulk is preserved at pseudobulk_v7_greenwald_2024_fine_backup.parquet.
This script writes pseudobulk_v7_greenwald_2024.parquet with coarse columns.

Final coarse labels: astrocyte, OPC, microglia, TAM, oligodendrocyte, neuron,
GABA_neuron, endothelial, T_cell, unassigned, malignant_AC, malignant_OPC,
malignant_NPC, malignant_MES, malignant (anything else that's clearly tumor
but no clean state assignment).
"""
import re
from pathlib import Path
import pandas as pd
import numpy as np

REQ = Path("/data/projects/atacseq/pilot_v5/requant")
FINE = REQ / "pseudobulk_v7_greenwald_2024_fine_backup.parquet"
OUT = REQ / "pseudobulk_v7_greenwald_2024.parquet"


def remap(label: str) -> str:
    L = label.lower()
    # malignant Neftel states
    if "mes" in L and ("hypox" in L or "mes-like" in L or L.startswith("mes") or "mesen" in L):
        return "malignant_MES"
    if "hypox" in L:
        return "malignant_MES"  # Hypoxic clusters are MES-like in Neftel schema
    if "ac-gliosis" in L or "ac progenitor" in L or ("ac" in L and "like" in L):
        return "malignant_AC"
    if "opc-npc" in L or "opc-like" in L or L == "opc-like":
        return "malignant_OPC"
    if "npc-neuronal" in L or "npc-like" in L:
        return "malignant_NPC"
    if any(k in L for k in ["proliferat", "cycling", "g1s", "g2m"]):
        return "malignant"
    if "undefined" in L or "ambiguous" in L:
        return "unassigned"
    # non-malignant lineages
    if "oligodendrocyte" in L or L.startswith("oligo"):
        return "oligodendrocyte"
    if "tam" in L or "macrophage" in L:
        return "TAM"
    if "microglia" in L:
        return "microglia"
    if "endothel" in L or "mural" in L or "pericyte" in L or "smc" in L:
        return "endothelial"
    if "t cell" in L or "cd4" in L or "cd8" in L or "treg" in L or "nk" in L:
        return "T_cell"
    if "neuron" in L or "exc " in L or "inh " in L or "l5" in L or "l6" in L or "l2" in L or "l3" in L or "l4" in L:
        if "gaba" in L or "inh" in L:
            return "GABA_neuron"
        return "neuron"
    if "astrocyte" in L:
        return "astrocyte"
    if "opc" in L:
        return "OPC"
    return "unassigned"


def main():
    df = pd.read_parquet(FINE)
    print(f"loaded fine-grained parquet: {df.shape}")

    # Each column is: greenwald_2024__<patient>__<fine_label>
    remap_dict = {}  # new_col_name -> list of old col names to sum
    for col in df.columns:
        parts = col.split("__", 2)
        if len(parts) != 3:
            print(f"  WARN: unparseable column {col}")
            continue
        cohort, patient, fine = parts
        coarse = remap(fine)
        new_col = f"{cohort}__{patient}__{coarse}"
        remap_dict.setdefault(new_col, []).append(col)

    print(f"remapped {len(df.columns)} fine → {len(remap_dict)} coarse columns")
    # Show counts per coarse label (across patients)
    coarse_label_counts = {}
    for new_col in remap_dict:
        lab = new_col.split("__", 2)[2]
        coarse_label_counts[lab] = coarse_label_counts.get(lab, 0) + 1
    print(f"patient-column counts per coarse label:")
    for lab, n in sorted(coarse_label_counts.items(), key=lambda x: -x[1]):
        print(f"  {lab}: {n}")

    # Sum the fine columns into coarse columns
    out_df = pd.DataFrame({
        new_col: df[old_cols].sum(axis=1).astype(np.int64).values
        for new_col, old_cols in remap_dict.items()
    })

    out_df.to_parquet(OUT, index=False)
    print(f"wrote {OUT} shape={out_df.shape}")


if __name__ == "__main__":
    main()
