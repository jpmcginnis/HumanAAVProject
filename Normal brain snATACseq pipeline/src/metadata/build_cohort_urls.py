#!/usr/bin/env python3
"""Build per-cohort data_url + paper_doi_or_pmid mapping, append to v18 → v19 metadata."""
import polars as pl
import urllib.request, json, re, time
from pathlib import Path

AT = Path('/Users/jpmcginnis1/Desktop/Normal cortex snATACseq atlas/atlas_v2.0.0')
V18 = AT / 'audit_2026-10-25/metadata_FINAL_v18.parquet'
V19 = AT / 'audit_2026-10-25/metadata_FINAL_v19.parquet'
V19_CSV = AT / 'audit_2026-10-25/metadata_FINAL_v19.csv'
V19_CSV_FULL = AT / 'audit_2026-10-25/metadata_FINAL_v19_fullcols.csv'

# Hand-curated for non-GEO cohorts (from the catalog audit doc)
MANUAL = {
    '4dn_hippo':         ('https://data.4dnucleome.org/experiment-set-replicates/?experiment_type.display_title=snATAC-seq&tissue=hippocampus',
                          'CATlas / Li-Ren 2023, Science DOI 10.1126/science.adf7044'),
    'hmba-bg':           ('s3://hmba-bican-sharing-802451596237-us-west-2/BasalGanglia_pre-print_ATAC/human/',
                          'HMBA basal ganglia preprint 2026.02.12.705594 (bioRxiv)'),
    'hmba-bg-retry':     ('s3://hmba-bican-sharing-802451596237-us-west-2/BasalGanglia_pre-print_ATAC/human/',
                          'HMBA basal ganglia preprint 2026.02.12.705594 (bioRxiv)'),
    'menon_corces_pd':   ('https://cellxgene.cziscience.com/collections/16876983-d454-43db-a8ec-77fafca9bf38',
                          'DOI 10.64898/2026.03.05.709922 (Menon/Corces 2026 PD multi-region atlas)'),
    'cxg-takeuchi-als':  ('https://figshare.com/articles/dataset/29628653',
                          'Takeuchi/Nagano 2026 ALS multiome preprint (figshare DOI 10.6084/m9.figshare.29628653)'),
    'cxg-lattke-ds':     ('https://cellxgene.cziscience.com/collections/ceb895f4-ff7f-403c-8b8c-83c6a3b5f5a0',
                          'Lattke et al. Down syndrome fetal cortex multiome'),
    'cxg-zhu-roussos':   ('https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204682',
                          'Zhu/Roussos developing cerebral cortex multi-omics (GSE204682)'),
    'cxg1':              ('https://cellxgene.cziscience.com/collections/a66d41a7-f1de-4b70-b983-7a4b790a0d5d',
                          'Braun/Linnarsson first-trimester neurodevelopment'),
    'seaad':             ('https://portal.brain-map.org/explore/seattle-alzheimers-disease; CELLxGENE 1ca90a2d-2943-483d-b678-b809bf464c30; s3://sea-ad-single-cell-profiling/',
                          'SEA-AD consortium (Allen Institute / Lein / Keene); multiple papers'),
    'veiga':             ('https://cellxgene.cziscience.com/collections/',
                          'Veiga adolescent cortex (CELLxGENE)'),
}

m = pl.read_parquet(V18)
cohorts = sorted(m['cohort_tag'].unique().to_list())
print(f'cohorts: {len(cohorts)}')

data_url = {}
paper = {}
for c in cohorts:
    if c in MANUAL:
        data_url[c], paper[c] = MANUAL[c]
    elif c.startswith('GSE') or re.match(r'^geo-GSE', c) or re.match(r'^gse\d+', c, re.I):
        # Normalize to GSE accession
        acc = c.replace('geo-','').replace('_margin','').replace('_C1','').replace('_C2','').replace('_catlas','')
        acc = re.match(r'(GSE\d+)', acc, re.I)
        acc = acc.group(1).upper() if acc else c
        data_url[c] = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}'
        paper[c] = f'{acc} — PMID TBD (fetch via GEO)'
    else:
        data_url[c] = 'TBD'
        paper[c] = 'TBD'

# Try to fetch PMIDs for GSE accessions in bulk via NCBI eutils
# elink: dbfrom=gds (GEO) → db=pubmed
gse_accs = sorted({v for k,v in data_url.items() if 'acc=GSE' in str(v)})
print(f'GSE accessions to look up: {len(gse_accs)}')

def gse_to_pmid(acc):
    try:
        # Step 1: GSE → GDS UID
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={acc}[Accession]&retmode=json'
        with urllib.request.urlopen(url, timeout=10) as r:
            data = json.load(r)
        ids = data.get('esearchresult',{}).get('idlist',[])
        if not ids: return None
        uid = ids[0]
        # Step 2: GDS UID → PubMed via elink
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=gds&db=pubmed&id={uid}&retmode=json'
        with urllib.request.urlopen(url, timeout=10) as r:
            data = json.load(r)
        linksets = data.get('linksets',[])
        for ls in linksets:
            for db in ls.get('linksetdbs',[]):
                if db.get('dbto') == 'pubmed' and db.get('links'):
                    return db['links'][0]
        return None
    except Exception as e:
        return None

# For each unique GSE accession, lookup
pmid_cache = {}
for i, url in enumerate(gse_accs):
    acc_match = re.search(r'acc=(GSE\d+)', url)
    if not acc_match: continue
    acc = acc_match.group(1)
    if acc in pmid_cache: continue
    pmid = gse_to_pmid(acc)
    pmid_cache[acc] = pmid
    print(f'  [{i+1}/{len(gse_accs)}] {acc} -> PMID {pmid}', flush=True)
    time.sleep(0.4)  # eutils rate-limit (3 req/sec max)

# Apply back to paper dict
for c in cohorts:
    url = data_url[c]
    acc_match = re.search(r'acc=(GSE\d+)', url)
    if acc_match:
        acc = acc_match.group(1)
        pmid = pmid_cache.get(acc)
        if pmid:
            paper[c] = f'{acc} — PMID {pmid}'
        # else: keep TBD

# Build join DF and merge
url_df = pl.DataFrame({
    'cohort_tag': cohorts,
    'data_url': [data_url[c] for c in cohorts],
    'paper_doi_or_pmid': [paper[c] for c in cohorts],
})
print()
print(url_df)

m2 = m.join(url_df, on='cohort_tag', how='left')
m2.write_parquet(V19)
useful = ['donor_id','study_id','atlas_usability_tier',
          'disease_status_final','disease_specific_final',
          'tissue_final','brain_region_final','region',
          'sample_type_final','age_final','sex_final','pmi_final',
          'modality_dmap','n_barcodes_raw_cohort','n_cells_usable_post_nnz1k_qc',
          'data_source','data_url','paper_doi_or_pmid','access_status',
          'disease_status_final_source','tissue_final_source',
          'brain_region_final_source','sample_type_final_source',
          'exclusion_reasons']
useful = [c for c in useful if c in m2.columns]
m2.select(useful).write_csv(V19_CSV)
m2.write_csv(V19_CSV_FULL)
print(f'\nwrote {V19_CSV}')
