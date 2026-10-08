"""QC Check #8: ENCODE blacklist + high-copy repeat intersect.

Reports:
  - how many of the 544K CATLAS atlas peaks overlap ENCODE blacklist hg38 v2
  - how many of the current top-50 pan-malignant shortlist overlap blacklist

Blacklist used: ENCODE blacklist hg38 v2 — the canonical set of regions
that produce spurious signal in chromatin-based assays (centromeres,
telomeres, high-copy satellite, assembly gaps).

Output: qc/outputs/08_blacklist_intersect.csv (per-peak intersect result)
        qc/outputs/08_blacklist_summary.txt
"""
import pandas as pd
import polars as pl
from pathlib import Path
import subprocess
import bisect

HERE = Path("/Users/jpmcginnis1/Desktop/GBM enhancer atlas 10-4-26/2026-10-6 dataset and results")
BLACKLIST = HERE / "reference/encode_blacklist_hg38.v2.bed"
MATRIX = HERE / "data/matrix/enhancer_candidate_matrix.parquet"
TOP50 = HERE / "reports/pan_malignant_top50_distal_selective.csv"
OUT_PEAK = HERE / "qc/outputs/08_blacklist_intersect.csv"
OUT_SUMMARY = HERE / "qc/outputs/08_blacklist_summary.txt"

# Download blacklist if not present
if not BLACKLIST.exists():
    print(f"[08] downloading ENCODE blacklist hg38 v2")
    URL = "https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz"
    subprocess.run(["curl", "-sL", URL, "-o", str(BLACKLIST) + ".gz"], check=True)
    subprocess.run(["gunzip", "-f", str(BLACKLIST) + ".gz"], check=True)
    print(f"  wrote {BLACKLIST}")

bl = pd.read_csv(BLACKLIST, sep="\t", header=None,
                 names=["chrom","start","end","reason"], usecols=[0,1,2,3])
print(f"[08] blacklist regions: {len(bl):,}  distinct chroms: {bl['chrom'].nunique()}")
print(f"[08] blacklist reasons: {bl['reason'].value_counts().to_dict()}")

# Build per-chromosome sorted interval lists for fast lookup
bl_by_chr = {}
for ch, grp in bl.groupby("chrom"):
    g = grp.sort_values("start").reset_index(drop=True)
    bl_by_chr[ch] = (g["start"].values, g["end"].values, g["reason"].values)

def overlaps_blacklist(chrom, start, end):
    if chrom not in bl_by_chr: return (False, None)
    starts, ends, reasons = bl_by_chr[chrom]
    # Find all blacklist intervals whose start < our end AND end > our start
    i = bisect.bisect_right(starts, end)
    for k in range(i-1, -1, -1):
        if ends[k] <= start: break
        if ends[k] > start:  # overlap
            return (True, reasons[k])
    return (False, None)

# === All 544K CATLAS peaks ===
print("[08] loading CATLAS atlas peaks (via matrix var)")
m = pl.read_parquet(MATRIX, columns=["peak_id","chrom","start","end"]).unique()
print(f"  atlas unique peaks: {m.height:,}")

atlas_rows = []
hits = 0
for row in m.iter_rows(named=True):
    o, reason = overlaps_blacklist(row["chrom"], row["start"], row["end"])
    if o:
        hits += 1
        atlas_rows.append({"peak_id": row["peak_id"], "chrom": row["chrom"],
                           "start": row["start"], "end": row["end"],
                           "blacklist_reason": reason})
print(f"[08] atlas peaks overlapping blacklist: {hits:,} / {m.height:,} ({100*hits/m.height:.2f}%)")

# === Top 50 shortlist ===
print("[08] checking top-50 pan-malignant shortlist")
top50 = pd.read_csv(TOP50)
top50_results = []
for _, r in top50.iterrows():
    o, reason = overlaps_blacklist(r["chrom"], int(r["start"]), int(r["end"]))
    top50_results.append({
        "peak_id": r["peak_id"],
        "chrom": r["chrom"],
        "start": r["start"],
        "end": r["end"],
        "nearest_gene": r.get("nearest_gene", ""),
        "in_blacklist": o,
        "blacklist_reason": reason or "",
    })
top50_df = pd.DataFrame(top50_results)
n_hits_top50 = int(top50_df["in_blacklist"].sum())
print(f"[08] top-50 shortlist peaks in blacklist: {n_hits_top50} / {len(top50_df)}")
if n_hits_top50 > 0:
    print(top50_df[top50_df["in_blacklist"]].to_string(index=False))

# Save output
atlas_df = pd.DataFrame(atlas_rows)
# Union: all top-50 (both pass and fail) + atlas-level hits
out = pd.concat([
    top50_df.assign(source="top50_shortlist"),
    atlas_df.assign(nearest_gene="", source="atlas_hit", in_blacklist=True),
], ignore_index=True, sort=False)
out.to_csv(OUT_PEAK, index=False)
print(f"[08] wrote {OUT_PEAK}")

with open(OUT_SUMMARY, "w") as f:
    f.write(f"ENCODE blacklist hg38 v2 — intersect results (2026-10-07)\n")
    f.write(f"=====================================================\n\n")
    f.write(f"Blacklist regions: {len(bl):,} (reasons: {bl['reason'].value_counts().to_dict()})\n")
    f.write(f"\nAtlas (544K CATLAS peaks that appear in matrix):\n")
    f.write(f"  total unique peaks: {m.height:,}\n")
    f.write(f"  overlapping blacklist: {hits:,} ({100*hits/m.height:.2f}%)\n")
    f.write(f"\nCurrent top-50 pan-malignant shortlist:\n")
    f.write(f"  overlapping blacklist: {n_hits_top50} / {len(top50_df)}\n")
    if n_hits_top50:
        f.write(f"\nBlacklist-flagged top-50 peaks (DROP from AAV cloning shortlist):\n")
        for _, r in top50_df[top50_df["in_blacklist"]].iterrows():
            f.write(f"  {r['peak_id']}  near {r.get('nearest_gene','?')}  reason={r['blacklist_reason']}\n")
print(f"[08] wrote {OUT_SUMMARY}")
