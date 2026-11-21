"""v1.1 patch module — adds CNV-aware annotations and arm-conditional Z.

Called by score.py (updated to import these functions).
"""
import numpy as np
import pandas as pd

# hg38 centromere positions (approximate, UCSC cytoBand midpoints)
CENTROMERES = {
    "chr1": 122026460, "chr2": 92188146, "chr3": 90772459, "chr4": 50404875,
    "chr5": 48802000, "chr6": 59774000, "chr7": 60145000, "chr8": 45190000,
    "chr9": 43389000, "chr10": 40199000, "chr11": 51644000, "chr12": 34770000,
    "chr13": 17700000, "chr14": 17800000, "chr15": 19000000, "chr16": 36835000,
    "chr17": 25263000, "chr18": 17900000, "chr19": 24400000, "chr20": 25844000,
    "chr21": 10900000, "chr22": 15400000, "chrX": 58000000, "chrY": 10400000,
}

# EGFR focal amplicon: EGFR gene is chr7:55019017-55207338 (hg38). Use ±2 Mb buffer.
EGFR_FOCAL_START = 52_000_000
EGFR_FOCAL_END   = 58_000_000

# chr9p arm end = centromere position of chr9
CHR9P_END = CENTROMERES["chr9"]


def annotate_cnv_class_and_arm(df):
    """Add chromosome_cnv_class and chromosome_arm columns to a peak DataFrame
    (must have chrom + start columns)."""
    chrom = df["chrom"].values
    start = df["start"].values

    # chromosome arm
    arm = np.full(len(df), "", dtype=object)
    for ch, cent in CENTROMERES.items():
        m = chrom == ch
        if not m.any(): continue
        is_p = start < cent
        arm[m & is_p] = f"{ch}p"
        arm[m & ~is_p] = f"{ch}q"
    df["chromosome_arm"] = arm

    # chromosome_cnv_class
    cls = np.full(len(df), "none", dtype=object)
    # chr7_egfr_focal: chr7 and in EGFR focal window
    m_chr7 = chrom == "chr7"
    m_egfr = m_chr7 & (start >= EGFR_FOCAL_START) & (start < EGFR_FOCAL_END)
    cls[m_egfr] = "chr7_egfr_focal"
    cls[m_chr7 & ~m_egfr] = "chr7_gain_broad"
    cls[chrom == "chr10"] = "chr10_loss"
    cls[(chrom == "chr9") & (start < CHR9P_END)] = "chr9p_loss"
    df["chromosome_cnv_class"] = cls
    return df


def arm_adjusted_z(df, z_col="daigle_z"):
    """Return a new Series with per-arm Z-scored values of z_col.
    df must have 'chromosome_arm' column."""
    z = df[z_col].astype(float).values
    arm = df["chromosome_arm"].values
    out = np.zeros_like(z)
    for a in np.unique(arm):
        m = arm == a
        if m.sum() < 10:
            # too few peaks on this arm for sensible stats → leave as raw
            out[m] = z[m]
            continue
        mu = z[m].mean()
        sd = z[m].std()
        if sd > 0:
            out[m] = (z[m] - mu) / sd
        else:
            out[m] = 0.0
    return pd.Series(out, index=df.index)


def patient_fraction_elevated(target_cols, col_meta, big_cpm, df_peaks,
                              elevated_class="chr7_gain_broad",
                              reference_class="none",
                              ratio_threshold=1.3):
    """For each peak, among patients where it has signal in the target label,
    what fraction have genome-wide `elevated_class` CPM ≥ ratio_threshold ×
    their own `reference_class` CPM?"""
    # Precompute per-patient genome-wide ratio
    cnv_class = df_peaks["chromosome_cnv_class"].values
    is_elev = (cnv_class == elevated_class) | (cnv_class == "chr7_egfr_focal")
    is_ref  = cnv_class == reference_class
    patient_ratios = {}
    for col in target_cols:
        colvals = big_cpm[col].values
        elev_sum = colvals[is_elev].sum()
        ref_sum = colvals[is_ref].sum()
        if ref_sum > 0:
            # normalize for size: elev CPM / n_elev vs ref CPM / n_ref
            elev_rate = elev_sum / is_elev.sum()
            ref_rate  = ref_sum / is_ref.sum()
            patient_ratios[col] = elev_rate / max(ref_rate, 1e-9)
        else:
            patient_ratios[col] = 1.0
    # For each peak, patient columns with signal:
    target_present = (big_cpm[target_cols] > 0.5).astype(np.int8)
    patient_is_elev = pd.Series({c: int(patient_ratios[c] >= ratio_threshold) for c in target_cols})
    # Fraction of elevated patients among those with signal
    n_with_signal = target_present.sum(axis=1)
    n_elev_with_signal = (target_present * patient_is_elev.values).sum(axis=1)
    frac = n_elev_with_signal / n_with_signal.replace(0, np.nan)
    return frac.fillna(0.0), patient_ratios
