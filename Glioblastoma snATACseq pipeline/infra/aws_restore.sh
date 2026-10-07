#!/usr/bin/env bash
# Restore the GBM Enhancer Atlas EC2 environment from an EBS snapshot.
#
# Usage:
#   ./infra/aws_restore.sh <snapshot-id> [--instance-type r6i.4xlarge] [--region us-east-1]
#
# What this does:
#   1. Creates a new 500 GB gp3 volume from the snapshot in us-east-1a.
#   2. Launches a new r6i.4xlarge Ubuntu 22.04 instance (defaults).
#   3. Attaches the restored volume to /dev/sdf.
#   4. SSHes in, mounts to /data, verifies the conda env + matrix artifacts.
#   5. Prints the public IP for subsequent work.
#
# The snapshot was created 2026-10-07 after the first full pipeline run
# finished (1.52M cells, 51 patients, 8 cohorts, matrix v2 complete).
# Keeping the snapshot costs ~$20/month for 500GB; restoring takes ~15 min.

set -euo pipefail

SNAPSHOT_ID="${1:-}"
INSTANCE_TYPE="${2:-r6i.4xlarge}"
REGION="${3:-us-east-1}"
AZ="${REGION}a"
AMI_ID="${AMI_ID:-ami-0c7217cdde317cfec}"  # Ubuntu 22.04 LTS us-east-1
KEY_NAME="${KEY_NAME:-bioinfo-key}"
SEC_GROUP="${SEC_GROUP:-sg-0ab95934ea158c30a}"

if [ -z "$SNAPSHOT_ID" ]; then
    echo "Usage: $0 <snapshot-id> [instance-type] [region]"
    echo ""
    echo "Known snapshot(s):"
    aws ec2 describe-snapshots --owner-ids self --region "$REGION" \
        --filters "Name=tag:Project,Values=GBMEnhancerAtlas" \
        --query 'Snapshots[].[SnapshotId,StartTime,VolumeSize,State,Description]' \
        --output table
    exit 2
fi

echo "=== Restoring from $SNAPSHOT_ID in $AZ ($INSTANCE_TYPE) ==="

# 1. Create volume from snapshot
echo "Creating 500 GB gp3 volume from snapshot..."
VOL_ID=$(aws ec2 create-volume \
    --region "$REGION" \
    --availability-zone "$AZ" \
    --snapshot-id "$SNAPSHOT_ID" \
    --volume-type gp3 \
    --size 500 \
    --tag-specifications 'ResourceType=volume,Tags=[{Key=Project,Value=GBMEnhancerAtlas},{Key=Source,Value=restored-snapshot}]' \
    --query 'VolumeId' --output text)
echo "  volume: $VOL_ID"

# Wait for availability
echo "  waiting for volume to become available..."
aws ec2 wait volume-available --region "$REGION" --volume-ids "$VOL_ID"

# 2. Launch new instance
echo "Launching new $INSTANCE_TYPE instance..."
INSTANCE_ID=$(aws ec2 run-instances \
    --region "$REGION" \
    --image-id "$AMI_ID" \
    --instance-type "$INSTANCE_TYPE" \
    --key-name "$KEY_NAME" \
    --security-group-ids "$SEC_GROUP" \
    --placement "AvailabilityZone=$AZ" \
    --tag-specifications 'ResourceType=instance,Tags=[{Key=Project,Value=GBMEnhancerAtlas},{Key=Source,Value=restored-snapshot}]' \
    --query 'Instances[0].InstanceId' --output text)
echo "  instance: $INSTANCE_ID"

echo "  waiting for instance to be running..."
aws ec2 wait instance-running --region "$REGION" --instance-ids "$INSTANCE_ID"

# 3. Attach volume
echo "Attaching $VOL_ID to $INSTANCE_ID as /dev/sdf..."
aws ec2 attach-volume --region "$REGION" --volume-id "$VOL_ID" \
    --instance-id "$INSTANCE_ID" --device /dev/sdf > /dev/null

echo "  waiting for attach..."
aws ec2 wait volume-in-use --region "$REGION" --volume-ids "$VOL_ID"

# 4. SSH + mount + verify
PUBLIC_IP=$(aws ec2 describe-instances --region "$REGION" \
    --instance-ids "$INSTANCE_ID" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
echo "  public IP: $PUBLIC_IP"

echo ""
echo "Waiting 60s for SSH to come up..."
sleep 60

echo "Mounting /data and verifying..."
ssh -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null \
    -i ~/atacseq_project/bioinfo-key.pem "ubuntu@$PUBLIC_IP" bash <<'REMOTE'
sudo mkdir -p /data
# nvme1n1 is the typical name for a secondary EBS volume on nitro instances
if lsblk | grep -q nvme1n1; then
    sudo mount /dev/nvme1n1 /data
else
    sudo mount /dev/xvdf /data
fi
df -h /data
echo "--- conda env ---"
ls /data/miniconda3/envs/atacseq/bin/python 2>&1 | head
echo "--- matrix artifacts ---"
ls -lh /data/projects/atacseq/claude/gbm-enhancer-atlas/matrix/ 2>&1 | head
REMOTE

echo ""
echo "=== DONE ==="
echo "Instance:  $INSTANCE_ID  (public IP $PUBLIC_IP)"
echo "Volume:    $VOL_ID"
echo "SSH:       ssh -i ~/atacseq_project/bioinfo-key.pem ubuntu@$PUBLIC_IP"
echo ""
echo "To stop (preserve EBS):     aws ec2 stop-instances --instance-ids $INSTANCE_ID"
echo "To terminate (delete all):  aws ec2 terminate-instances --instance-ids $INSTANCE_ID"
