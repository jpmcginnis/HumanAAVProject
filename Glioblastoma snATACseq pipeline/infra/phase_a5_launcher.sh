#!/bin/bash
# Phase A.5/A.6/A.7 auto-launcher. Fires when all required downloads are:
#   (1) present and above the min expected size
#   (2) pass gzip integrity check (where applicable)
#   (3) no wget still active for them
#   (4) no snakemake already running
#   (5) not already fully processed through to catlas_quantified.h5ad
# Idempotent. Called by cron every 20 minutes.

set -euo pipefail
RAW=/data/projects/atacseq/raw
PROC=/data/projects/atacseq/processed
LAUNCHER_LOG=/data/projects/atacseq/logs/phase_a5_launcher.log

log(){ echo "[$(date -u +%H:%M:%S)] $*" >> "$LAUNCHER_LOG"; }

# File:min_size pairs (bytes)
H5AD_GZ_FILES=(
  "$RAW/gbm_space/GBM_space_ATAC_filtered_peaks.h5ad.gz:20000000000"
  "$RAW/gbm_space/GBM_space_snRNA.h5ad.gz:10000000000"
)
TAR_FILES=(
  "$RAW/gse138794_guo/GSE138794_RAW.tar:500000000"
  "$RAW/gse276177_khan_astro/GSE276177_RAW.tar:20000000000"
)
AUX_FILES=(
  "$RAW/gbm_space/E-MTAB-17183.sdrf.txt"
)

# Size check
for pair in "${H5AD_GZ_FILES[@]}" "${TAR_FILES[@]}"; do
  f="${pair%%:*}"; min="${pair##*:}"
  if [ ! -f "$f" ]; then log "wait: $f missing"; exit 0; fi
  sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$sz" -lt "$min" ]; then log "wait: $f too small ($sz < $min)"; exit 0; fi
done
for f in "${AUX_FILES[@]}"; do
  [ -f "$f" ] || { log "wait: aux missing: $f"; exit 0; }
done

# No wget still active for these inputs?
if pgrep -af "wget.*(gbm_space|gse138794|GSE276177|E-MTAB-17183)" > /dev/null; then
  log "wait: wget still active"; exit 0
fi

# CRITICAL: gzip integrity for h5ad.gz files. This was the gap that let
# a corrupt snRNA through earlier. gunzip -t reads the whole file (~2 min
# on 11 GB) so we cache the result in a .verified marker.
for pair in "${H5AD_GZ_FILES[@]}"; do
  f="${pair%%:*}"
  marker="${f}.verified"
  if [ ! -f "$marker" ] || [ "$marker" -ot "$f" ]; then
    log "verifying gzip integrity: $f"
    if gunzip -t "$f" 2>/dev/null; then
      touch "$marker"
      log "  OK"
    else
      log "  CORRUPT — deleting + exiting (wget will re-pull on next cron)"
      rm -f "$f"
      # Launch a fresh wget so the file re-populates
      TS=$(date -u +%Y%m%d_%H%M%S)
      name=$(basename "$f")
      nohup setsid wget -q -O "$f" \
        "https://www.ebi.ac.uk/biostudies/files/E-MTAB-17183/${name}" \
        > "/data/projects/atacseq/logs/wget_${name}_relaunch_${TS}.log" 2>&1 &
      disown
      log "  re-pull launched (pid=$!)"
      exit 0
    fi
  fi
done

# tar integrity (quick — just checks header)
for pair in "${TAR_FILES[@]}"; do
  f="${pair%%:*}"
  marker="${f}.verified"
  if [ ! -f "$marker" ] || [ "$marker" -ot "$f" ]; then
    log "verifying tar integrity: $f"
    if tar -tf "$f" >/dev/null 2>&1; then
      touch "$marker"
      log "  OK"
    else
      log "  CORRUPT — deleting + exiting"
      rm -f "$f"
      exit 0
    fi
  fi
done

# Extract tars if not yet done
for ext_pair in \
  "$RAW/gse138794_guo/GSE138794_RAW.tar:$RAW/gse138794_guo/.extracted" \
  "$RAW/gse276177_khan_astro/GSE276177_RAW.tar:$RAW/gse276177_khan_astro/.extracted"; do
  tar_path="${ext_pair%%:*}"; marker="${ext_pair##*:}"
  if [ ! -f "$marker" ]; then
    log "extracting $(basename $tar_path)"
    cd "$(dirname $tar_path)"
    tar -xf "$(basename $tar_path)"
    touch "$marker"
    log "  done"
  fi
done

if pgrep -af 'bin/snakemake' > /dev/null; then
  log "wait: snakemake already running"; exit 0
fi

all_done=1
for ct in gbm_space gse138794_guo gse276177_khan_astro; do
  [ -f "$PROC/$ct/catlas_quantified.h5ad" ] || { all_done=0; break; }
done
if [ "$all_done" = "1" ]; then
  log "done: all new cohorts already quantified"; exit 0
fi

TS=$(date -u +%Y%m%d_%H%M%S)
RUN_LOG=/data/projects/atacseq/logs/overnight_${TS}_phase_a5.log
log "LAUNCHING Phase A.5/A.6/A.7 — all files present + verified + extracted"

cd /data/projects/atacseq/claude/gbm-enhancer-atlas
source /data/miniconda3/etc/profile.d/conda.sh && conda activate atacseq

if ! grep -q 'PHASE_A_COHORTS = \["gbm_space"\]' Snakefile; then
  python3 -c "
import re
with open('Snakefile') as f: t = f.read()
t = re.sub(r'PHASE_A_COHORTS = \[\]', 'PHASE_A_COHORTS = [\"gbm_space\"]', t, count=1)
with open('Snakefile','w') as f: f.write(t)
"
  log "patched Snakefile: added gbm_space to PHASE_A_COHORTS"
fi

rm -f matrix/enhancer_candidate_matrix.parquet matrix/candidate_scores.parquet \
      matrix/top_candidates.parquet matrix/surprises.parquet \
      matrix/dataset_issues.parquet matrix/provenance.parquet
rm -rf reports/out/*

tmux kill-session -t phase_a5 2>/dev/null || true
tmux new-session -d -s phase_a5 "bash -c '
cd /data/projects/atacseq/claude/gbm-enhancer-atlas
source /data/miniconda3/etc/profile.d/conda.sh && conda activate atacseq
export GBM_PHASE=A2
snakemake --cores 8 --keep-going --restart-times 1 --rerun-incomplete --rerun-triggers mtime > $RUN_LOG 2>&1
'"
log "launched in tmux session 'phase_a5'"
