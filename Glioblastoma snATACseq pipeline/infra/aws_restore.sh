#!/usr/bin/env bash
# Restore a prior pipeline run from an EBS snapshot and resume where we left off.
#
# Usage:
#     ./infra/aws_restore.sh snap-01234abcd
#     ./infra/aws_restore.sh                  # auto-pick the most recent snapshot
#                                             # tagged project=gbm-enhancer-atlas
#
# The snapshot becomes a new volume (same AZ as the EC2), attached to the
# instance at /dev/sdf. aws_spinup.sh is then invoked with EBS_VOLUME_ID set.
# Takes ~2-5 minutes end-to-end (snapshot → volume creation is immediate; the
# volume is usable while blocks lazy-load in the background).

set -euo pipefail

INSTANCE_ID="i-000dddf603f4bde12"
REGION="us-east-1"
AZ="us-east-1a"
PROJECT_TAG="gbm-enhancer-atlas"

step() { printf "\n[%s] >>> %s\n" "$(date +%H:%M:%S)" "$1"; }

SNAP="${1:-}"
if [[ -z "$SNAP" ]]; then
    step "No snapshot id passed; picking most recent snapshot tagged project=$PROJECT_TAG"
    SNAP=$(aws ec2 describe-snapshots --region "$REGION" --owner-ids self \
        --filters "Name=tag:project,Values=$PROJECT_TAG" \
        --query 'reverse(sort_by(Snapshots,&StartTime))[0].SnapshotId' --output text)
    if [[ "$SNAP" == "None" || -z "$SNAP" ]]; then
        echo "    No snapshot found. Run infra/aws_spinup.sh for a fresh volume instead."
        exit 1
    fi
    echo "    using $SNAP"
fi

step "Creating volume from snapshot $SNAP in $AZ"
VOL=$(aws ec2 create-volume --region "$REGION" --availability-zone "$AZ" \
    --snapshot-id "$SNAP" --volume-type gp3 \
    --tag-specifications "ResourceType=volume,Tags=[{Key=Name,Value=atacseq-data},{Key=project,Value=$PROJECT_TAG},{Key=restored_from,Value=$SNAP}]" \
    --query 'VolumeId' --output text)
echo "    new volume: $VOL"
aws ec2 wait volume-available --volume-ids "$VOL" --region "$REGION"

step "Resuming EC2 with EBS_VOLUME_ID=$VOL"
EBS_VOLUME_ID="$VOL" bash "$(dirname "$0")/aws_spinup.sh"
