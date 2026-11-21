#!/usr/bin/env bash
# Shut the pipeline down cleanly.
#
# Three modes (set MODE=... or pass as $1):
#
#   stop      — stop the EC2, keep the EBS attached. Cheapest for short pauses
#               (hours to days). Idle cost = EBS only (~$40/mo for 500GB gp3).
#               Resume with `./aws_spinup.sh EBS_VOLUME_ID=<vol-...>`.
#
#   snapshot  — stop the EC2, snapshot the data EBS to S3-backed snapshot,
#               wait for the snapshot to complete, then DELETE the live
#               volume. Idle cost drops from $0.08/GB/mo to $0.05/GB/mo
#               (snapshots are incremental, often much cheaper than that).
#               Resume with `./aws_restore.sh <snap-...>`.
#               Use for pauses of a month+. ← this is the default.
#
#   full      — same as snapshot, plus delete the EC2 instance's root
#               volume mapping too. Only if you're done with the project
#               for the quarter.

set -euo pipefail

INSTANCE_ID="i-000dddf603f4bde12"
REGION="us-east-1"
PROJECT_TAG="gbm-enhancer-atlas"

MODE="${1:-${MODE:-snapshot}}"

step() { printf "\n[%s] >>> %s\n" "$(date +%H:%M:%S)" "$1"; }

step "Stopping EC2 $INSTANCE_ID (mode=$MODE)"
aws ec2 stop-instances --instance-ids "$INSTANCE_ID" --region "$REGION" >/dev/null
aws ec2 wait instance-stopped --instance-ids "$INSTANCE_ID" --region "$REGION"
echo "    stopped"

# Find data volumes attached to this instance (exclude the root /dev/sda1 or /dev/xvda)
VOLUMES=$(aws ec2 describe-volumes --region "$REGION" \
  --filters "Name=attachment.instance-id,Values=$INSTANCE_ID" \
  --query 'Volumes[?Attachments[0].Device!=`/dev/sda1` && Attachments[0].Device!=`/dev/xvda`].VolumeId' \
  --output text)

if [[ -z "$VOLUMES" ]]; then
    echo "No data volumes attached. Nothing further to do."
    exit 0
fi

echo "Data volumes to process: $VOLUMES"

case "$MODE" in

  stop)
    step "Mode=stop — keeping EBS attached"
    echo "To resume next time:"
    for v in $VOLUMES; do echo "    EBS_VOLUME_ID=$v  ./infra/aws_spinup.sh"; done
    ;;

  snapshot|full)
    for VOL in $VOLUMES; do
      TS=$(date -u +%Y%m%dT%H%M%SZ)
      DESC="gbm-enhancer-atlas pipeline snapshot ${TS} from ${VOL}"
      step "Snapshotting $VOL → S3-backed snapshot"
      SNAP=$(aws ec2 create-snapshot --region "$REGION" \
        --volume-id "$VOL" --description "$DESC" \
        --tag-specifications "ResourceType=snapshot,Tags=[{Key=project,Value=$PROJECT_TAG},{Key=source_volume,Value=$VOL},{Key=created,Value=$TS}]" \
        --query 'SnapshotId' --output text)
      echo "    snapshot: $SNAP (pending — waiting for completion)"

      aws ec2 wait snapshot-completed --snapshot-ids "$SNAP" --region "$REGION"
      SNAP_SIZE_GB=$(aws ec2 describe-snapshots --snapshot-ids "$SNAP" --region "$REGION" \
        --query 'Snapshots[0].VolumeSize' --output text)
      echo "    snapshot complete: $SNAP (${SNAP_SIZE_GB} GB logical)"

      step "Detaching + deleting live volume $VOL"
      # Detach first (instance is stopped, so this is immediate)
      aws ec2 detach-volume --volume-id "$VOL" --region "$REGION" >/dev/null 2>&1 || true
      aws ec2 wait volume-available --volume-ids "$VOL" --region "$REGION"
      aws ec2 delete-volume --volume-id "$VOL" --region "$REGION"
      echo "    deleted"

      step "Done — resume info"
      echo "    To restore this run exactly where we left off:"
      echo "        ./infra/aws_restore.sh $SNAP"
      echo ""
      echo "    Snapshot cost estimate: ~\$$(python3 -c "print(f'{0.05*$SNAP_SIZE_GB:.2f}')") /month "
      echo "    (actual is often lower — snapshots are incremental and compressed)"
    done
    ;;

  *)
    echo "Unknown mode '$MODE' — use: stop, snapshot, or full"
    exit 1
    ;;
esac

if [[ "$MODE" == "full" ]]; then
    step "Terminating EC2 instance"
    aws ec2 terminate-instances --instance-ids "$INSTANCE_ID" --region "$REGION" >/dev/null
    aws ec2 wait instance-terminated --instance-ids "$INSTANCE_ID" --region "$REGION"
    echo "    terminated. To resume, aws_restore.sh will provision a fresh EC2."
fi

echo ""
echo "SHUTDOWN COMPLETE"
