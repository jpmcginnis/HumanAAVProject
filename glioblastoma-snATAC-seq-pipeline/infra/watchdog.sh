#!/usr/bin/env bash
# Pipeline watchdog — runs on the EC2 every 20 min via cron.
# Writes a status snapshot to /data/projects/atacseq/logs/watchdog.log.
# If snakemake is NOT running AND the latest overnight log is not "complete",
# re-launches snakemake in the `work` tmux session once (with retry counter).
#
# Install: `crontab -e` on the EC2, add:
#   */20 * * * * /data/projects/atacseq/claude/gbm-enhancer-atlas/infra/watchdog.sh

set -uo pipefail
LOG=/data/projects/atacseq/logs/watchdog.log
RETRY_FILE=/data/projects/atacseq/logs/.watchdog_retries
MAX_RETRIES=3

ts() { date +'%Y-%m-%d %H:%M:%S'; }
echo "" >>"$LOG"
echo "=== [$(ts)] watchdog tick ===" >>"$LOG"

# 1. Status snapshot
echo "DISK: $(df -h /data | tail -1)" >>"$LOG"
SNAKE_PIDS=$(pgrep -af 'bin/snakemake' | head -5)
if [[ -n "$SNAKE_PIDS" ]]; then
    echo "SNAKE_ALIVE: yes" >>"$LOG"
    echo "$SNAKE_PIDS" >>"$LOG"
else
    echo "SNAKE_ALIVE: no" >>"$LOG"
fi

LATEST_LOG=$(ls -t /data/projects/atacseq/logs/overnight_*.log 2>/dev/null | head -1)
if [[ -n "$LATEST_LOG" ]]; then
    echo "LATEST_LOG: $LATEST_LOG  size=$(du -h $LATEST_LOG | awk '{print $1}')" >>"$LOG"
    echo "--- last 10 lines ---" >>"$LOG"
    tail -10 "$LATEST_LOG" >>"$LOG"
fi

# 2. Decide if we need to restart
if [[ -n "$SNAKE_PIDS" ]]; then
    echo "DECISION: snakemake running — nothing to do" >>"$LOG"
    exit 0
fi

# Is the latest overnight log a clean finish?
if [[ -n "$LATEST_LOG" ]] && tail -40 "$LATEST_LOG" | grep -qE '(^\s*all)|( N jobs done\.|Nothing to be done|Finished job 0\.)' ; then
    # If snakemake exited with "Finished job 0." or "N jobs done" the run is complete.
    if tail -5 "$LATEST_LOG" | grep -qE 'Nothing to be done|jobs done|Finished job 0' ; then
        echo "DECISION: pipeline completed cleanly — not restarting" >>"$LOG"
        exit 0
    fi
fi

# Pipeline died with an error. Check retry counter.
RETRIES=0
if [[ -f "$RETRY_FILE" ]]; then RETRIES=$(cat "$RETRY_FILE"); fi
if (( RETRIES >= MAX_RETRIES )); then
    echo "DECISION: $RETRIES retries exhausted — leaving dead for JP to investigate" >>"$LOG"
    exit 1
fi
echo $(( RETRIES + 1 )) > "$RETRY_FILE"

echo "DECISION: snakemake dead, retry #$(( RETRIES + 1 )) — relaunching in tmux work" >>"$LOG"

# Make sure tmux is alive
tmux has-session -t work 2>/dev/null || tmux new-session -d -s work

NEW_LOG="/data/projects/atacseq/logs/overnight_$(date +%Y%m%d_%H%M%S)_watchdog.log"
CMD="cd /data/projects/atacseq/claude/gbm-enhancer-atlas && \
     export PATH=/data/miniconda3/bin:\$PATH && \
     export MAMBA_ROOT_PREFIX=/data/miniconda3 && \
     export GBM_PHASE=${GBM_PHASE:-A} && \
     micromamba run -n atacseq snakemake --cores all --keep-going --restart-times 2 \
       --rerun-incomplete 2>&1 | tee $NEW_LOG"
tmux send-keys -t work:0 "$CMD" Enter
echo "RELAUNCHED: new log at $NEW_LOG" >>"$LOG"
