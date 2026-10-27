#!/usr/bin/env python3
"""Build per-patient metadata CSV for GBM atlas, parallel format to NB atlas metadata.
Extracts (cohort, patient, cell_types_present) from matrix column names and joins
with cohort_attribution.csv for publication / accession / DOI / platform info.
"""
import polars as pl
import re
from pathlib import Path
from collections import defaultdict

GBM = Path('/Users/jpmcginnis1/Desktop/GBM WT enhancer atlas/GBM WT snATACseq atlas')
MATRIX = GBM / 'data/matrix/enhancer_candidate_matrix.parquet'
COH    = GBM / 'cohort_attribution.csv'
OUT_CSV = GBM / 'metadata_2026-10-27.csv'
OUT_CSV_FULL = GBM / 'metadata_2026-10-27_fullcols.csv'

cohort_df = pl.read_csv(COH)
cohort_lookup = {r['cohort_slug']: r for r in cohort_df.iter_rows(named=True)}

# Parse column names: <cohort>__<patient>__<cell_type>
m_cols = pl.read_parquet(MATRIX, n_rows=1).columns[1:]  # skip peak_id
pat = re.compile(r'^(.+?)__(.+?)__(.+)$')
per_pt = defaultdict(lambda: {'cell_types': set(), 'has_malignant': False, 'has_tam': False})
for c in m_cols:
    mt = pat.match(c)
    if not mt: continue
    cohort, patient, ct = mt.groups()
    key = (cohort, patient)
    per_pt[key]['cell_types'].add(ct)
    if 'malignant' in ct: per_pt[key]['has_malignant'] = True
    if ct == 'TAM': per_pt[key]['has_tam'] = True

print(f'unique (cohort, patient) pairs: {len(per_pt)}')

ACCESSION_URL = {
    'E-MTAB-17183': 'https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-17183',
    'GSE276177': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276177',
    'GSE165037': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165037',
    'GSE163654': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163654',
    'GSE139136': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE139136',
    'GSE138794': 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138794',
    'HRA004942': 'https://ngdc.cncb.ac.cn/gsa-human/browse/HRA004942',
    'GDC':       'https://portal.gdc.cancer.gov/legacy-archive/search/f?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Single%20Cell%20Analysis%22%5D%7D%7D%5D%7D',
}

def accession_to_url(acc):
    """Pick first known accession token in the raw accession string."""
    for key, url in ACCESSION_URL.items():
        if key in str(acc):
            return url
    if 'GDC' in str(acc) or 'TCGA' in str(acc):
        return ACCESSION_URL['GDC']
    return ''

rows = []
for (cohort, patient), info in sorted(per_pt.items()):
    ca = cohort_lookup.get(cohort, {})
    tier = 'PASS' if info['has_malignant'] else 'EXCLUDE_no_malignant'
    rows.append({
        'donor_id': patient,
        'study_id': cohort,
        'atlas_usability_tier': tier,
        'disease_status_final': 'glioblastoma IDH-wildtype',
        'tissue_final': 'brain',
        'brain_region_final': 'glioblastoma tumor',
        'sample_type_final': 'tumor',
        'modality_dmap': ca.get('platform', ''),
        'data_source': ca.get('accession', ''),
        'data_url': accession_to_url(ca.get('accession', '')),
        'paper_doi_or_pmid': ca.get('publication_cite', ''),
        'access_status': ca.get('data_access', ''),
        'n_cells_cohort_total_pre_qc': ca.get('n_cells_total'),
        'n_malignant_cells_cohort_total': ca.get('n_malignant_cells'),
        'n_patients_cohort': ca.get('n_patients_included'),
        'cell_types_present': ', '.join(sorted(info['cell_types'])),
        'has_malignant_track': info['has_malignant'],
        'has_tam_track': info['has_tam'],
        'genome_build': ca.get('genome_build', 'hg38'),
        'cohort_notes': ca.get('notes', ''),
    })

df = pl.DataFrame(rows)
print(f'per-patient rows: {df.shape[0]}')
print(df.head(5))
print(df.group_by('study_id').agg(pl.len().alias('n_patients')).sort('n_patients', descending=True))

# Trimmed CSV = all columns except cell_types_present (which is huge) + cohort_notes
trim = df.drop(['cell_types_present','cohort_notes'])
trim.write_csv(OUT_CSV)
df.write_csv(OUT_CSV_FULL)
print(f'\nwrote {OUT_CSV}')
print(f'wrote {OUT_CSV_FULL}')
