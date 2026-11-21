#!/bin/bash
# Phase A.7: Sundaram GBM add-on. Fires ONLY when:
#   (1) Phase A.5 atlas is complete (matrix/candidate_scores.parquet exists)
#   (2) Sundaram tar is downloaded, verified, and extracted
#   (3) No snakemake currently running
#   (4) Sundaram not already integrated
# Zero-risk to Phase A.5: this launcher does NOT touch Phase A.5 or any of its outputs.
# Called by cron every 20 minutes (independent of phase_a5_launcher.sh).

set -euo pipefail
RAW=/data/projects/atacseq/raw
PROC=/data/projects/atacseq/processed
PROJ=/data/projects/atacseq/claude/gbm-enhancer-atlas
LOG_F=/data/projects/atacseq/logs/phase_a7_launcher.log
MATRIX=$PROJ/matrix
log(){ echo "[$(date -u +%H:%M:%S)] $*" >> "$LOG_F"; }

# (1) Phase A.5 complete? (candidate_scores.parquet AND index.html exist)
if [ ! -f "$MATRIX/candidate_scores.parquet" ] || [ ! -f "$PROJ/reports/out/index.html" ]; then
  log "wait: Phase A.5 not yet complete"
  exit 0
fi

# (2) Sundaram tar present + full enough?
SUNDARAM_TAR=$RAW/sundaram_gbm_gdc/full_cancer_scatacseq
if [ ! -f "$SUNDARAM_TAR" ]; then log "wait: Sundaram tar missing"; exit 0; fi
sz=$(stat -c%s "$SUNDARAM_TAR" 2>/dev/null || echo 0)
if [ "$sz" -lt 3000000000 ]; then log "wait: Sundaram tar too small ($sz)"; exit 0; fi

# (3) Not actively being downloaded?
if pgrep -af "wget.*full_cancer_scatacseq" > /dev/null; then
  log "wait: Sundaram wget still active"
  exit 0
fi

# (4) Verify gzip integrity
MARKER=$RAW/sundaram_gbm_gdc/.verified
if [ ! -f "$MARKER" ] || [ "$MARKER" -ot "$SUNDARAM_TAR" ]; then
  log "verifying Sundaram tar.gz integrity"
  if gunzip -t "$SUNDARAM_TAR" 2>/dev/null; then
    touch "$MARKER"
    log "  OK"
  else
    log "  CORRUPT — removing + exiting (next run will trigger re-download)"
    rm -f "$SUNDARAM_TAR"
    exit 0
  fi
fi

# (5) Extract GBM subset if not yet done
EXTRACTED=$RAW/sundaram_gbm_gdc/.extracted
if [ ! -f "$EXTRACTED" ]; then
  log "extracting GBM subset from Sundaram pan-cancer tar"
  cd $RAW/sundaram_gbm_gdc
  tar -xzf full_cancer_scatacseq --wildcards "Cancer_scATACseq_data/scATAC_GBMx_*" || {
    log "extraction failed — aborting"
    exit 1
  }
  touch "$EXTRACTED"
  n=$(find . -name "scATAC_GBMx_*fragments.tsv.gz" | wc -l)
  log "  extracted $n GBM fragment files"
fi

# (6) No snakemake running?
if pgrep -af 'bin/snakemake' > /dev/null; then
  log "wait: snakemake already running (Phase A.5 cleanup?)"
  exit 0
fi

# (7) Already done?
if [ -f "$PROC/sundaram_gbm/catlas_quantified.h5ad" ]; then
  log "done: Sundaram already quantified"
  exit 0
fi

# (8) Fire
TS=$(date -u +%Y%m%d_%H%M%S)
RUN_LOG=/data/projects/atacseq/logs/overnight_${TS}_phase_a7.log
log "LAUNCHING Phase A.7 — Sundaram GBM add-on"
cd "$PROJ"
source /data/miniconda3/etc/profile.d/conda.sh && conda activate atacseq

# Patch Snakefile: add sundaram_gbm to PHASE_A_COHORTS (fragment-mode, not matrix)
if ! grep -q 'sundaram_gbm' Snakefile; then
  python3 -c "
import re
with open('Snakefile') as f: t = f.read()
# Append sundaram_gbm to PHASE_A_COHORTS
t = re.sub(r'PHASE_A_COHORTS = \[(.*?)\]', r'PHASE_A_COHORTS = [\1, \"sundaram_gbm\"]', t, count=1)
# Clean up leading comma if PHASE_A_COHORTS was empty
t = t.replace('PHASE_A_COHORTS = [, \"sundaram_gbm\"]', 'PHASE_A_COHORTS = [\"sundaram_gbm\"]')
with open('Snakefile','w') as f: f.write(t)
"
  log "patched Snakefile: added sundaram_gbm to PHASE_A_COHORTS"
fi

# Clear matrix so matrix_build re-runs with 8 cohorts
rm -f "$MATRIX/enhancer_candidate_matrix.parquet" "$MATRIX/candidate_scores.parquet" \
      "$MATRIX/top_candidates.parquet" "$MATRIX/surprises.parquet" \
      "$MATRIX/dataset_issues.parquet" "$MATRIX/provenance.parquet"
rm -rf "$PROJ/reports/out/"*

tmux kill-session -t phase_a7 2>/dev/null || true
tmux new-session -d -s phase_a7 "bash -c '
cd $PROJ
source /data/miniconda3/etc/profile.d/conda.sh && conda activate atacseq
export GBM_PHASE=A2
snakemake --cores 8 --keep-going --restart-times 1 --rerun-incomplete --rerun-triggers mtime > $RUN_LOG 2>&1
'"
log "Phase A.7 launched in tmux session 'phase_a7'"
