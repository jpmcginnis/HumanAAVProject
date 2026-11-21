#!/usr/bin/env bash
# AWS EC2 spin-up for the GBM enhancer atlas pipeline.
# One command, idempotent: starts instance, mounts /data, activates conda, attaches tmux.
# Follows JP's restart checklist from tools-and-workflow.md.

set -euo pipefail

INSTANCE_ID="i-000dddf603f4bde12"
REGION="us-east-1"
SECURITY_GROUP="sg-0ab95934ea158c30a"
SSH_KEY="$HOME/atacseq_project/bioinfo-key.pem"
TMUX_SESSION="work"
CONDA_ENV="atacseq"
S3_BUCKET="jpm-atacseq-archive-2026"

# EBS volume: provision if missing. JP's note: original vol-0e0af156ca5270d37 deleted 2026-05-04.
# Pass EBS_VOLUME_ID=vol-xxxx as env var if a persistent volume exists; otherwise provision.
EBS_VOLUME_ID="${EBS_VOLUME_ID:-}"
EBS_SIZE_GB="${EBS_SIZE_GB:-2000}"   # 2TB default; override for cheaper 500GB

step() { printf "\n[%s] >>> %s\n" "$(date +%H:%M:%S)" "$1"; }

# --- 1. Start instance ---
step "Starting EC2 instance $INSTANCE_ID"
aws ec2 start-instances --instance-ids "$INSTANCE_ID" --region "$REGION" >/dev/null
aws ec2 wait instance-running --instance-ids "$INSTANCE_ID" --region "$REGION"

PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" --region "$REGION" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' --output text)
step "Instance running at $PUBLIC_IP"

# --- 2. Add current home IP to security group ---
MY_IP=$(curl -sS https://checkip.amazonaws.com | tr -d '\n')
step "Authorizing SSH from $MY_IP"
aws ec2 authorize-security-group-ingress \
    --group-id "$SECURITY_GROUP" --region "$REGION" \
    --protocol tcp --port 22 --cidr "$MY_IP/32" 2>/dev/null \
    || echo "    (already authorized, continuing)"

# --- 3. EBS volume handling ---
if [[ -z "$EBS_VOLUME_ID" ]]; then
    step "No EBS_VOLUME_ID provided — provisioning fresh ${EBS_SIZE_GB}GB gp3 volume"
    EBS_VOLUME_ID=$(aws ec2 create-volume \
        --availability-zone us-east-1a --region "$REGION" \
        --size "$EBS_SIZE_GB" --volume-type gp3 \
        --tag-specifications 'ResourceType=volume,Tags=[{Key=Name,Value=atacseq-data},{Key=project,Value=gbm-enhancer-atlas}]' \
        --query 'VolumeId' --output text)
    echo "    New volume: $EBS_VOLUME_ID"
    echo "    SAVE THIS — pass as EBS_VOLUME_ID=$EBS_VOLUME_ID on next invocation"
    aws ec2 wait volume-available --volume-ids "$EBS_VOLUME_ID" --region "$REGION"
    FRESH_VOLUME=1
else
    step "Using existing EBS volume $EBS_VOLUME_ID"
    FRESH_VOLUME=0
fi

# Attach volume
step "Attaching volume to /dev/sdf"
aws ec2 attach-volume \
    --volume-id "$EBS_VOLUME_ID" --instance-id "$INSTANCE_ID" \
    --device /dev/sdf --region "$REGION" >/dev/null 2>&1 \
    || echo "    (already attached, continuing)"

sleep 5

# --- 4. Remote setup ---
step "Running remote setup on the instance"

REMOTE_SETUP=$(cat <<EOF
set -euo pipefail

# Mount volume
if ! mountpoint -q /data; then
    sudo mkdir -p /data
    DEVICE=\$(lsblk -no NAME,SIZE | awk '/${EBS_SIZE_GB}G/ || /2T/ {print "/dev/"\$1; exit}')
    if [ "$FRESH_VOLUME" = "1" ]; then
        echo "    Formatting fresh volume \$DEVICE"
        sudo mkfs.ext4 -F "\$DEVICE"
    fi
    sudo mount "\$DEVICE" /data
    UUID=\$(sudo blkid -s UUID -o value "\$DEVICE")
    echo "    Mounted \$DEVICE at /data (UUID=\$UUID)"
    # Update fstab for next boot
    sudo sed -i '/\\/data/d' /etc/fstab
    echo "UUID=\$UUID /data ext4 defaults,nofail 0 2" | sudo tee -a /etc/fstab >/dev/null
fi

sudo chown -R ubuntu:ubuntu /data

# Conda
if [ ! -d /data/miniconda3 ]; then
    echo "    Installing miniconda to /data/miniconda3"
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p /data/miniconda3
fi
source /data/miniconda3/etc/profile.d/conda.sh

# Env — use the in-repo pinned env.yml (snapatac2, infercnvpy, scglue, pymc,
# cellxgene-census, playwright). Fall back to the S3 copy only if the repo
# hasn't been cloned yet (first-ever bootstrap).
if ! conda env list | grep -q "^${CONDA_ENV} "; then
    echo "    Creating conda env ${CONDA_ENV}"
    ENV_YML="/data/projects/atacseq/claude/gbm-enhancer-atlas/infra/atacseq_env.yml"
    if [ ! -f "\$ENV_YML" ]; then
        echo "    Repo env.yml not found yet; using S3 bootstrap copy"
        aws s3 cp "s3://${S3_BUCKET}/env/atacseq_env.yml" /tmp/atacseq_env.yml && ENV_YML=/tmp/atacseq_env.yml
    fi
    conda env create -f "\$ENV_YML" -n ${CONDA_ENV}
fi

# Project directory
mkdir -p /data/projects/atacseq/{raw,processed,results,logs,scripts,ref}
cd /data/projects/atacseq

# Sync pipeline repo — pulled from S3 (not GitHub) so pre-publication work stays private
if [ ! -d gbm-enhancer-atlas ]; then
    echo "    Pulling project tarball from S3"
    mkdir -p /data/projects/atacseq/claude
    cd /data/projects/atacseq/claude
    aws s3 cp s3://${S3_BUCKET}/code/gbm-enhancer-atlas.tar.gz /tmp/gbm-enhancer-atlas.tar.gz
    tar -xzf /tmp/gbm-enhancer-atlas.tar.gz
    ls -la gbm-enhancer-atlas/
fi

# Attach or create tmux
tmux has-session -t ${TMUX_SESSION} 2>/dev/null || tmux new-session -d -s ${TMUX_SESSION}
echo ""
echo "Ready. To attach: ssh -i ${SSH_KEY} ubuntu@${PUBLIC_IP} -t 'tmux attach -t ${TMUX_SESSION}'"
EOF
)

ssh -o StrictHostKeyChecking=no -i "$SSH_KEY" "ubuntu@$PUBLIC_IP" bash -s <<< "$REMOTE_SETUP"

step "SPIN-UP COMPLETE"
echo ""
echo "    Connect:        ssh -i $SSH_KEY ubuntu@$PUBLIC_IP"
echo "    Attach tmux:    ssh -i $SSH_KEY ubuntu@$PUBLIC_IP -t 'tmux attach -t $TMUX_SESSION'"
echo "    EBS volume ID:  $EBS_VOLUME_ID  (save for next spin-up)"
echo ""
echo "    Next:           cd /data/projects/atacseq/claude/gbm-enhancer-atlas && snakemake --cores all"
