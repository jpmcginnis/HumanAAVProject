# GBM Atlas v4 — UCSC Track Hub Plan (post-EBS-restore follow-up)

**Status:** Queued. Executes when `snap-092225222f7f7ed63` finishes restoring (started 2026-10-17 19:23 UTC, standard tier, 24-72 hr).

## Goal

Host per-cell-type bigWig / bigBed tracks on a public S3 bucket and have the UCSC Genome Browser auto-load them via a track hub, so clicking any peak's UCSC link in the v4 HTML atlas opens with 14 cell-type-specific GBM accessibility tracks visible alongside the ENCODE / phastCons / JASPAR reference tracks already pre-enabled.

End state: a user clicking "View in UCSC Genome Browser" from v4 sees, at the top of the UCSC view, 14 colored signal tracks showing the actual atlas-derived Tn5 accessibility across the peak for astrocyte / microglia / oligodendrocyte / OPC / endothelial / neuron / GABA_neuron / T_cell / TAM / malignant (pooled) + 4 Neftel subtypes. This lets them pick the focal ~500 bp AAV insert within a wider peak based on where signal actually concentrates in their target cell type.

## Steps (sequential, from EBS restore completion)

### 1. Restore + pull per-cohort fragments (~1-2 hr once restored)
```
# Verify restore completed
aws ec2 describe-snapshot-tier-status --filters "Name=snapshot-id,Values=snap-092225222f7f7ed63" --region us-east-1

# Create volume from snapshot
aws ec2 create-volume --snapshot-id snap-092225222f7f7ed63 --availability-zone us-east-1a --volume-type gp3

# Launch r6i.4xlarge in us-east-1a; attach volume
# SSH in, mount volume

# Expected contents of restored volume:
#  /data/projects/atacseq/
#    cohort_fragments/<cohort>/<sample>.fragments.tsv.gz  (8 cohorts)
#    per_cell_h5ads/<cohort>/catlas_quantified.h5ad       (includes greenwald)
#    pilot_v5/scoring_scratch/                             (intermediate matrices)
```

### 2. Generate per-cell-type pooled bigWigs (~4-8 hr compute)
For each of 14 cell types:
- Collect all cell barcodes assigned to that cell type across all 8 cohorts (from the catlas_quantified.h5ad obs tables + the atlas's cell-type labels)
- Filter each cohort's fragments.tsv to just those barcodes
- Pool filtered fragments across cohorts into one bedGraph per cell type
- Convert to bigWig with `bedGraphToBigWig <bedGraph> hg38.chrom.sizes <output.bigWig>`
- Scale by total cells per cell type for CPM-equivalent y-axis

Expected output: 14 bigWigs, 100-500 MB each, total ~2-5 GB.

Script skeleton (`scripts/build_percelltype_bigwigs.py`):
```python
import pandas as pd, subprocess, pyBigWig
from pathlib import Path

CT_LIST = ['astrocyte','microglia','oligodendrocyte','OPC','endothelial',
           'neuron','GABA_neuron','T_cell','TAM','malignant',
           'malignant_AC','malignant_MES','malignant_NPC','malignant_OPC']

# For each cell type, build a BED of per-base Tn5 insertion sites, then convert to bigWig
for ct in CT_LIST:
    # Load cell barcodes assigned to this cell type (from all cohort obs tables)
    barcodes = load_barcodes_for_celltype(ct)
    # Filter all cohort fragments.tsv files, extract midpoints, write bedGraph
    bedgraph = f'/tmp/{ct}.bedGraph'
    with open(bedgraph, 'w') as out:
        for cohort in COHORTS:
            fragments_path = f'/data/projects/atacseq/cohort_fragments/{cohort}/all.fragments.tsv.gz'
            # stream, filter by barcode, emit (chrom, pos, pos+1, 1) per Tn5 cut
            ...
    # Sort + aggregate
    subprocess.run(['sort', '-k1,1', '-k2,2n', bedgraph, '-o', bedgraph], check=True)
    # Collapse adjacent same-position cuts
    collapsed = f'/tmp/{ct}.collapsed.bedGraph'
    # ...
    # Convert
    subprocess.run(['bedGraphToBigWig', collapsed, 'hg38.chrom.sizes', f'/data/output/{ct}.bigWig'], check=True)
```

### 3. Also build per-cell-type bigBed from the union peak set + strength values (~15 min)
For each cell type: filter union peaks to those with reasonable strength in that cell type, output as BED9 with `itemRgb` color by cell type. Lets users show "all atlas-called peaks for this cell type" as a track distinct from the signal bigWig.

### 4. Upload to public S3 bucket (~30 min)
```
aws s3 mb s3://jpm-gbm-atlas-trackhub --region us-east-1
aws s3 cp /data/output/ s3://jpm-gbm-atlas-trackhub/hg38/ --recursive --acl public-read
```

Also upload the hub config files:
```
s3://jpm-gbm-atlas-trackhub/
  hub.txt
  genomes.txt
  hg38/
    trackDb.txt
    astrocyte.bigWig
    microglia.bigWig
    ... (12 more)
    astrocyte_peaks.bigBed
    ... (13 more)
```

### 5. Write hub config files (~15 min)

**hub.txt:**
```
hub GBM_snATAC_atlas_v4
shortLabel GBM Atlas v4 (McGinnis)
longLabel GBM snATAC-seq enhancer atlas v4 — 48 patients, 8 cohorts, 1.96M union peaks (McGinnis 2026)
genomesFile genomes.txt
email jpmcginnis1@gmail.com
descriptionUrl description.html
```

**genomes.txt:**
```
genome hg38
trackDb hg38/trackDb.txt
```

**hg38/trackDb.txt:**
```
track GBMatlas_bigwigs
container multiWig
shortLabel GBM atlas per-celltype signal
longLabel GBM atlas v4 per-cell-type Tn5 accessibility (CPM-scaled, pooled across evaluable patients)
type bigWig 0 100
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
autoScale on
alwaysZero on
priority 1

    track malignant_signal
    parent GBMatlas_bigwigs on
    shortLabel malignant (pooled)
    longLabel GBM malignant (pooled Neftel subtypes) — pooled across 40 patients, 8 cohorts
    bigDataUrl malignant.bigWig
    type bigWig 0 100
    color 220,20,60

    track astrocyte_signal
    parent GBMatlas_bigwigs on
    shortLabel astrocyte
    bigDataUrl astrocyte.bigWig
    color 34,139,34
    ... (12 more subtracks with distinct colors)

track GBMatlas_peaks
container multiWig
shortLabel GBM atlas per-celltype called peaks
longLabel GBM atlas v4 per-cell-type peak calls (strength >= floor)
type bigBed 9
visibility pack
priority 2
    ... (14 subtracks of bigBed peaks)
```

### 6. Update v4 HTML to load the hub (~10 min, one-line change)
In `v4_build_html.py`, change the UCSC URL construction to include `&hubUrl=https://s3.amazonaws.com/jpm-gbm-atlas-trackhub/hub.txt`. The custom red block for the current peak is kept — it's complementary, highlighting WHICH peak the user clicked on. The 14 cell-type signal tracks now also load.

### 7. Re-emit v4 HTML, send to user, ship

## Expected final UCSC view

For a given peak, user sees from top to bottom:
1. GBM atlas peak (red block — the current clicked peak)
2. 14 GBM atlas per-cell-type signal tracks (bigWig, with the tab's cell type auto-highlighted)
3. 14 GBM atlas per-cell-type peak call tracks (bigBed)
4. GENCODE comprehensive (genes)
5. ENCODE cCREs
6. ENCODE H3K27ac
7. phastCons 100-way conservation
8. JASPAR 2024 TFBS
9. RepeatMasker

Picking the focal ~500 bp AAV insert then becomes: find where in the peak the GBM-cell-type signal track peaks + where phastCons is highest + where JASPAR TFBS cluster is + avoid repeats. All visible in one UCSC view.

## Cost estimate

- S3 Standard-IA storage for ~5 GB track hub: ~$0.06/month
- S3 egress (data transfer): $0.09/GB → trivial for interactive browsing (~1 MB per peak view)
- Total: ~$1-2/year

Essentially free. The infrastructure persists after upload with no maintenance.

## Related: v5 trigger

Switching from Poisson-approx penetrance to exact per-cell penetrance (the swap-in that originally motivated the EBS restore) can also happen in the same session: re-run the scoring pipeline with exact per-cell counts, verify rank order barely shifts, re-emit as v5-final if the shifts matter OR annotate v4 as the final numeric values if they don't.

## Also persist per-cell h5ads to S3 at this time

Per the memory [[feedback-persist-per-cell-h5ads-to-s3]], the restored per-cell union-peak h5ads MUST be uploaded to S3 at this same step so this restore-nightmare doesn't repeat:
```
aws s3 sync /data/projects/atacseq/per_cell_h5ads/ s3://jpm-gbm-atlas-intermediates/per_cell_union_h5ads/ --storage-class STANDARD_IA
```
Verify completion, THEN tear down the EC2 instance.
