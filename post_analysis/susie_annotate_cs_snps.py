"""
Annotates SuSiE credible set SNPs using Ensembl VEP REST API.
Reads all _credible_sets.tsv from SuSiE_results_6wind, queries VEP in batches,
and outputs consequence, gene name, and biotype for each unique SNP.
"""

import os
import json
import time
import requests
import pandas as pd

SUSIE_RESULTS = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_results_6wind'
OUTPUT_FILE   = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_analysis/credible_set_annotations.tsv'

ENSEMBL_VEP_URL = 'https://grch37.rest.ensembl.org/vep/human/id'
BATCH_SIZE      = 200
MAX_RETRIES     = 3
RETRY_DELAY     = 10

HEADERS = {
    'Content-Type': 'application/json',
    'Accept':       'application/json',
}


def vep_batch(rsids):
    payload = json.dumps({'ids': rsids})
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            resp = requests.post(ENSEMBL_VEP_URL, headers=HEADERS, data=payload, timeout=60)
            if resp.status_code == 429:
                wait = int(resp.headers.get('Retry-After', RETRY_DELAY))
                print(f'  Rate limited, waiting {wait}s...')
                time.sleep(wait)
                continue
            resp.raise_for_status()
            return resp.json()
        except Exception as e:
            print(f'  Attempt {attempt}/{MAX_RETRIES} failed: {e}')
            if attempt < MAX_RETRIES:
                time.sleep(RETRY_DELAY)
    return []


def parse_vep(variant):
    most_severe  = variant.get('most_severe_consequence', '')
    transcripts  = variant.get('transcript_consequences', [])

    # Among transcripts matching the most severe consequence, prefer protein_coding
    best = None
    for tc in transcripts:
        if most_severe in tc.get('consequence_terms', []):
            if best is None or tc.get('biotype') == 'protein_coding':
                best = tc

    if best is None and transcripts:
        best = transcripts[0]

    return {
        'most_severe_consequence': most_severe,
        'gene_symbol':             best.get('gene_symbol', '') if best else '',
        'gene_id':                 best.get('gene_id',     '') if best else '',
        'biotype':                 best.get('biotype',     '') if best else '',
    }


# ── Collect credible set SNPs (keep all rows) ──────────────────────────────────
all_rows = []

for hit_dir in sorted(os.listdir(SUSIE_RESULTS)):
    cs_file = os.path.join(SUSIE_RESULTS, hit_dir, f'{hit_dir}_credible_sets.tsv')
    if not os.path.exists(cs_file):
        continue
    df = pd.read_csv(cs_file, sep='\t')
    df.insert(0, 'hit', hit_dir)
    all_rows.append(df)

cs_df = pd.concat(all_rows, ignore_index=True)
unique_rsids = cs_df['ID'].unique().tolist()
print(f'Total rows: {len(cs_df)}  |  Unique rsIDs to annotate: {len(unique_rsids)}')

# ── Query Ensembl VEP in batches ───────────────────────────────────────────────
annotations = {}
for i in range(0, len(unique_rsids), BATCH_SIZE):
    batch = unique_rsids[i:i + BATCH_SIZE]
    print(f'Querying VEP batch {i + 1}–{i + len(batch)}...')
    results = vep_batch(batch)
    for v in results:
        rsid = v.get('id', v.get('input', ''))
        if rsid:
            annotations[rsid] = parse_vep(v)
    time.sleep(1)

# ── Join annotations onto all credible set rows ────────────────────────────────
ann_df = pd.DataFrame([
    {'ID': rsid, **annotations.get(rsid, {
        'most_severe_consequence': 'not_found',
        'gene_symbol': '',
        'gene_id':     '',
        'biotype':     '',
    })}
    for rsid in unique_rsids
])

out_df = cs_df.merge(ann_df, on='ID', how='left')
out_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

print(f'\nSaved {len(out_df)} rows → {OUTPUT_FILE}')
print(out_df[['hit', 'CS', 'ID', 'PIP', 'most_severe_consequence', 'gene_symbol', 'biotype']].to_string(index=False))
