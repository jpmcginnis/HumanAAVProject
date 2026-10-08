"""QC Check #9: TSS re-classification using ENCODE cCRE promoter-like set.

Problem: current distal filter is |d| > 2 kb from nearest RefSeq TSS. Alternate
TSSs 2-5 kb from annotated ones look like enhancers but aren't. Need an
orthogonal "is this actually a promoter?" check.

Uses: ENCODE SCREEN cCRE "promoter-like signature" (PLS) set for hg38.
(FANTOM5 CAGE requires significantly larger download; cCRE PLS is the
curated promoter set from the same ENCODE project, ~41K elements.)

Output: qc/outputs/09_tss_reclassification.csv — one row per top-100 candidate
        with flags for ENCODE promoter cCRE overlap.
"""
import pandas as pd
import polars as pl
from pathlib import Path
import subprocess
import bisect

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results")
PROMOTER_BED = HERE / "reference/encode_ccre_promoter_hg38.bed"
TOP100 = HERE / "reports/pan_malignant_top100_all.csv"
TOP50 = HERE / "reports/pan_malignant_top50_distal_selective.csv"
OUT = HERE / "qc/outputs/09_tss_reclassification.csv"
OUT_SUMMARY = HERE / "qc/outputs/09_tss_reclassification_summary.txt"

# Download the ENCODE cCRE promoter-like set
if not PROMOTER_BED.exists():
    print(f"[09] downloading ENCODE SCREEN cCRE promoter-like set (hg38)")
    URL = "https://downloads.wenglab.org/V3/GRCh38-PLS.bed"
    subprocess.run(["curl", "-sL", URL, "-o", str(PROMOTER_BED)], check=True)
    print(f"  wrote {PROMOTER_BED}")

prom = pd.read_csv(PROMOTER_BED, sep="\t", header=None,
                   names=["chrom","start","end","name"], usecols=[0,1,2,3])
print(f"[09] ENCODE cCRE promoter regions: {len(prom):,}")

# Per-chrom sorted lookup
prom_by_chr = {}
for ch, grp in prom.groupby("chrom"):
    g = grp.sort_values("start").reset_index(drop=True)
    prom_by_chr[ch] = (g["start"].values, g["end"].values, g["name"].values)

def overlaps_promoter(chrom, start, end):
    if chrom not in prom_by_chr: return (False, None)
    starts, ends, names = prom_by_chr[chrom]
    i = bisect.bisect_right(starts, end)
    for k in range(i-1, -1, -1):
        if ends[k] <= start: break
        if ends[k] > start:
            return (True, names[k])
    return (False, None)

def check_df(df, label):
    results = []
    hits = 0
    for _, r in df.iterrows():
        o, name = overlaps_promoter(r["chrom"], int(r["start"]), int(r["end"]))
        if o: hits += 1
        results.append({
            "source": label,
            "peak_id": r["peak_id"],
            "chrom": r["chrom"],
            "start": r["start"],
            "end": r["end"],
            "nearest_gene": r.get("nearest_gene", ""),
            "dist_to_tss_signed": r.get("dist_to_tss_signed", None),
            "encode_promoter_overlap": o,
            "encode_promoter_name": name or "",
        })
    print(f"[09] {label}: {hits} / {len(df)} overlap ENCODE promoter cCRE")
    return pd.DataFrame(results), hits

top100 = pd.read_csv(TOP100)
top50 = pd.read_csv(TOP50)

t100_df, t100_hits = check_df(top100, "top100_all")
t50_df, t50_hits = check_df(top50, "top50_distal_selective")

out_df = pd.concat([t100_df, t50_df], ignore_index=True)
out_df.to_csv(OUT, index=False)
print(f"[09] wrote {OUT}")

with open(OUT_SUMMARY, "w") as f:
    f.write("ENCODE SCREEN cCRE promoter-like overlap — top candidates (2026-10-07)\n")
    f.write("=====================================================================\n\n")
    f.write(f"Promoter cCRE regions used: {len(prom):,}\n\n")
    f.write(f"Top 100 (unfiltered): {t100_hits} / {len(top100)} overlap ENCODE promoter\n")
    f.write(f"Top 50 (post-distal filter): {t50_hits} / {len(top50)} overlap ENCODE promoter\n\n")
    if t50_hits > 0:
        f.write("Promoter-overlap top-50 candidates (DROP from distal-enhancer shortlist):\n")
        for _, r in t50_df[t50_df["encode_promoter_overlap"]].iterrows():
            dist = r.get("dist_to_tss_signed", "?")
            dist_kb = f"{dist/1000:+.1f} kb" if isinstance(dist, (int, float)) and pd.notna(dist) else dist
            f.write(f"  {r['peak_id']}  near {r.get('nearest_gene','?')} ({dist_kb})  cCRE={r['encode_promoter_name']}\n")
    f.write(f"\nNote: distal filter in v2 report used |dist_to_tss| ≥ 2 kb to nearest\n")
    f.write(f"RefSeq transcript. cCRE promoter is an orthogonal annotation using\n")
    f.write(f"DNase + H3K4me3 + H3K27ac signatures; catches alternate TSSs the\n")
    f.write(f"RefSeq transcript distance misses.\n")
print(f"[09] wrote {OUT_SUMMARY}")
