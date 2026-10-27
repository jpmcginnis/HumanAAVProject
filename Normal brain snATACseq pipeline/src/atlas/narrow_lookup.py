#!/usr/bin/env python3
"""narrow_lookup.py — given a Design A region (from nb_atlas_v2 HTML), return the overlapping
narrow Design B peaks + their hg38 sequences.

Design A peaks (v1.0.0 scored parquet) are often 1-10 kb wide — too big to clone into AAV.
Design B peaks (narrow_union_B.bed, 2.06M peaks) are median 518 bp, max 7 kb — AAV-cloneable.
This script bridges the two: pick a top-ranked region from the A HTML viewer, get back the
narrow sub-peaks within it + their actual DNA.

USAGE:
  python3 narrow_lookup.py chr5:71864459-71866513
  python3 narrow_lookup.py chr5 71864459 71866513
  python3 narrow_lookup.py --batch peaks.txt         # one peak_id per line
  python3 narrow_lookup.py chr5:71864459-71866513 --fasta        # output FASTA format

By default sequences come from the UCSC REST API (requires internet, no setup, ~1 sec per peak).
To use a local hg38 FASTA instead (faster, offline), pass --local-fasta /path/to/hg38.fa
(requires `pip install pyfaidx` and ~3 GB FASTA + ~2 MB index).

The narrow_union_B.bed is auto-downloaded from S3 on first run and cached at
~/.cache/nb_atlas/narrow_union_B.bed (47 MB).
"""
import argparse
import os
import re
import sys
import time
from pathlib import Path

CACHE_DIR = Path.home() / ".cache" / "nb_atlas"
BED_LOCAL = CACHE_DIR / "narrow_union_B.bed"
S3_BED = "s3://jpm-atacseq-archive-2026/data_processed_nb_integrated/v1.1/narrow_union_B.bed"
UCSC_URL = "https://api.genome.ucsc.edu/getData/sequence"


def ensure_bed():
    """Download narrow_union_B.bed from S3 on first run; cache locally."""
    if BED_LOCAL.exists() and BED_LOCAL.stat().st_size > 10_000_000:
        return BED_LOCAL
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[narrow_lookup] fetching narrow_union_B.bed from S3 (one-time, ~47 MB)...", file=sys.stderr)
    import subprocess
    r = subprocess.run(["aws", "s3", "cp", S3_BED, str(BED_LOCAL), "--only-show-errors"],
                       capture_output=True, text=True)
    if r.returncode != 0 or not BED_LOCAL.exists():
        print(f"[narrow_lookup] S3 fetch failed. Point to a local copy with env var NARROW_BED=/path/to/narrow_union_B.bed", file=sys.stderr)
        print(f"stderr: {r.stderr}", file=sys.stderr)
        sys.exit(2)
    return BED_LOCAL


def parse_peak(s: str):
    """Parse 'chr5:71864459-71866513' into (chrom, start, end)."""
    m = re.match(r"^(chr[\w]+):(\d+)-(\d+)$", s.strip())
    if not m:
        raise ValueError(f"could not parse peak '{s}' (expected chrN:start-end)")
    return m.group(1), int(m.group(2)), int(m.group(3))


def load_bed(bed_path: Path):
    """Load the narrow_union_B.bed into a per-chrom list of (start, end) sorted."""
    from collections import defaultdict
    by_chrom = defaultdict(list)
    with open(bed_path) as f:
        for line in f:
            toks = line.rstrip().split("\t")
            if len(toks) < 3:
                continue
            by_chrom[toks[0]].append((int(toks[1]), int(toks[2])))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def overlapping_peaks(by_chrom, chrom, start, end):
    """Return B peaks that overlap [start, end) on this chrom."""
    if chrom not in by_chrom:
        return []
    peaks = by_chrom[chrom]
    hits = []
    # Linear scan: 2M peaks total but one chrom has at most ~200K; fast enough
    for b_start, b_end in peaks:
        if b_start >= end:
            break  # sorted; past our region
        if b_end > start:
            hits.append((chrom, b_start, b_end))
    return hits


def fetch_sequence_ucsc(chrom, start, end, retries=3):
    """Fetch DNA via UCSC REST API."""
    import urllib.request
    import urllib.parse
    import json
    url = f"{UCSC_URL}?genome=hg38;chrom={chrom};start={start};end={end}"
    last_err = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=10) as r:
                data = json.loads(r.read().decode())
            seq = data.get("dna", "")
            if not seq:
                raise RuntimeError(f"UCSC returned no DNA for {chrom}:{start}-{end}: {data}")
            return seq.upper()
        except Exception as e:
            last_err = e
            time.sleep(1 + attempt)
    raise RuntimeError(f"UCSC fetch failed after {retries} tries for {chrom}:{start}-{end}: {last_err}")


def fetch_sequence_local(fa, chrom, start, end):
    """Fetch DNA from a local hg38 pyfaidx handle (half-open coords)."""
    return str(fa[chrom][start:end]).upper()


def lookup_one(chrom, start, end, by_chrom, fasta_handle=None, fasta_mode="ucsc"):
    hits = overlapping_peaks(by_chrom, chrom, start, end)
    results = []
    for h_chrom, h_start, h_end in hits:
        if fasta_mode == "ucsc":
            seq = fetch_sequence_ucsc(h_chrom, h_start, h_end)
        else:
            seq = fetch_sequence_local(fasta_handle, h_chrom, h_start, h_end)
        results.append({
            "peak_id": f"{h_chrom}:{h_start}-{h_end}",
            "chrom": h_chrom, "start": h_start, "end": h_end,
            "width_bp": h_end - h_start,
            "sequence": seq,
        })
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("peak", nargs="*", help="peak_id (chrN:start-end) or chrom start end")
    ap.add_argument("--batch", metavar="FILE", help="file with one peak_id per line")
    ap.add_argument("--fasta", action="store_true", help="output in FASTA format")
    ap.add_argument("--tsv", action="store_true", help="output in TSV format (one row per narrow peak)")
    ap.add_argument("--local-fasta", metavar="HG38.FA",
                    help="use a local hg38 FASTA (requires pyfaidx; faster for bulk)")
    args = ap.parse_args()

    # Determine input peaks
    inputs = []
    if args.batch:
        with open(args.batch) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                inputs.append(parse_peak(line))
    elif len(args.peak) == 1:
        inputs.append(parse_peak(args.peak[0]))
    elif len(args.peak) == 3:
        inputs.append((args.peak[0], int(args.peak[1]), int(args.peak[2])))
    else:
        ap.print_help()
        sys.exit(1)

    # Load BED
    bed_path = Path(os.environ.get("NARROW_BED", "")) if os.environ.get("NARROW_BED") else ensure_bed()
    print(f"[narrow_lookup] loading narrow_union_B.bed...", file=sys.stderr)
    by_chrom = load_bed(bed_path)
    total_peaks = sum(len(v) for v in by_chrom.values())
    print(f"[narrow_lookup]   {total_peaks:,} narrow B peaks loaded across {len(by_chrom)} chroms", file=sys.stderr)

    # Set up sequence source
    fasta_handle = None
    fasta_mode = "ucsc"
    if args.local_fasta:
        try:
            from pyfaidx import Fasta
            fasta_handle = Fasta(args.local_fasta)
            fasta_mode = "local"
            print(f"[narrow_lookup] using local FASTA: {args.local_fasta}", file=sys.stderr)
        except ImportError:
            print("[narrow_lookup] pyfaidx not installed; `pip install pyfaidx`. Falling back to UCSC API.", file=sys.stderr)

    # Output
    if args.tsv:
        print("input_peak_id\tnarrow_peak_id\tchrom\tstart\tend\twidth_bp\tsequence")

    for in_chrom, in_start, in_end in inputs:
        in_id = f"{in_chrom}:{in_start}-{in_end}"
        in_width = in_end - in_start
        hits = lookup_one(in_chrom, in_start, in_end, by_chrom, fasta_handle, fasta_mode)
        if args.tsv:
            for h in hits:
                print(f"{in_id}\t{h['peak_id']}\t{h['chrom']}\t{h['start']}\t{h['end']}\t{h['width_bp']}\t{h['sequence']}")
        elif args.fasta:
            for h in hits:
                print(f">{h['peak_id']}  (narrow sub-peak of {in_id}, {h['width_bp']} bp)")
                # Wrap at 60 chars (standard FASTA line length)
                seq = h['sequence']
                for i in range(0, len(seq), 60):
                    print(seq[i:i+60])
        else:
            # Pretty summary
            print(f"\n=== Input A peak: {in_id}  ({in_width:,} bp) ===")
            if not hits:
                print("  (no overlapping Design B narrow peaks — unusual, flag this)")
                continue
            print(f"  {len(hits)} Design B narrow peak(s) within this region:\n")
            for h in hits:
                print(f"  {h['peak_id']}  ({h['width_bp']} bp)")
                seq = h['sequence']
                # Print first 60 chars with ... if truncated for terminal display
                if len(seq) <= 180:
                    for i in range(0, len(seq), 60):
                        print(f"    {seq[i:i+60]}")
                else:
                    print(f"    {seq[:60]}...")
                    print(f"    ...{seq[-60:]}  (full {len(seq)} bp — use --fasta for complete seq)")
                print()


if __name__ == "__main__":
    main()
