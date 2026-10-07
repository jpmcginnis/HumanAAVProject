# GBM Enhancer Atlas — EC2 restore cheat sheet

Full point-in-time backup of the pipeline environment lives in **AWS EBS snapshot `snap-04deb97b5e367bedf`** (us-east-1).
This captures everything on the `/data` 500 GB volume as of 2026-10-07 02:42 UTC:

- Conda env (`/data/miniconda3`) with snapatac2, pyMC, polars, etc.
- Reference files: CATLAS atlas (40 GB), GBmap (11 GB), CATLAS BED
- Raw fragments for all 8 cohorts (237 GB)
- Intermediate per-sample h5ads (QC + CNV + labels)
- Processed cohort `catlas_quantified.h5ad` × 8
- Final matrix (`enhancer_candidate_matrix.parquet` 3.5 GB)
- Pan-malignant outputs (`pan_malignant_matrix.parquet` + scores)
- Logs

## Monthly cost

- Snapshot storage: **~$20/month** for 500 GB (incremental to S3 standard).
- Instance + attached volume: deleted. $0 when idle.

## Spin back up (15 min)

### Prereqs (one-time)
- AWS CLI configured with credentials for account `287883095118`, us-east-1 default.
- SSH key: `~/atacseq_project/bioinfo-key.pem` (chmod 400).

### Restore
```bash
cd "~/Desktop/GBM enhancer atlas 10-4-26/code"
./infra/aws_restore.sh snap-04deb97b5e367bedf
```

This script:
1. Creates a new 500 GB gp3 volume from the snapshot in us-east-1a.
2. Launches a new r6i.4xlarge Ubuntu 22.04 instance.
3. Attaches the restored volume to `/dev/sdf`.
4. SSHes in, mounts to `/data`, verifies the conda env + matrix artifacts.
5. Prints the public IP for subsequent work.

Takes ~15 min end-to-end. Hourly cost while running: ~$1.07 ($1.008 r6i.4xlarge + ~$0.08 EBS).

### Shut down again
```bash
# Stop (keeps EBS at ~$30/mo, fast restart):
aws ec2 stop-instances --instance-ids <new-instance-id> --region us-east-1

# Or terminate + delete volume (back to just the snapshot at ~$20/mo):
aws ec2 terminate-instances --instance-ids <new-instance-id> --region us-east-1
aws ec2 delete-volume --volume-id <new-volume-id> --region us-east-1
```

## Alternative — no EC2 needed for most queries

The main artifacts are also on your Google Drive at:
`~/Library/CloudStorage/GoogleDrive-jpmcginnis1@gmail.com/My Drive/AAV Gene Therapy/Enhancer ATACseq projects/October 2026 WT GBM analysis/`

Everything there can be queried locally with polars/pandas. **You only need to spin up EC2 for:**
- Processing a NEW cohort from scratch (fetch → ingest → quantify)
- Rebuilding the matrix with different cell-type grouping (requires per-cell data + 128GB RAM)
- Running the heavy PyMC ADVI stage of scoring on top candidates

For everything else — top-N queries, nearest-gene annotations, custom filters, cross-cell-type pivots — stay local.

## Backups

- **Snapshot** (`snap-04deb97b5e367bedf`): full /data volume, 500 GB → ~$20/mo. Permanent.
- **S3** (`s3://jpm-atacseq-archive-2026/final_atlas_2026_10_07/`): matrix + reports + 8 h5ads, 79 GB → ~$1.80/mo.
- **Google Drive**: full matrix + reports + 8 h5ads (via s3 sync).
- **GitHub**: all source code at https://github.com/jpmcginnis/HumanAAVProject (`Glioblastoma snATACseq pipeline/`).

Any one of these could restore the project; snapshot is fastest.
