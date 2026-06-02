"""
Annotates unique fine-mapping candidate SNPs with Ensembl VEP (GRCh37)
and produces a bar plot of genomic consequence categories.
"""

import os
import json
import time
import requests
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CANDIDATES_FILE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
OUTPUT_DIR      = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results'
ANNOT_FILE      = os.path.join(OUTPUT_DIR, 'fine_mapping_candidates_annotated.tsv')
PLOT_FILE       = os.path.join(OUTPUT_DIR, 'fine_mapping_consequence_barplot.png')

ENSEMBL_VEP_URL = 'https://grch37.rest.ensembl.org/vep/human/id'
BATCH_SIZE      = 200
MAX_RETRIES     = 3
RETRY_DELAY     = 10
HEADERS = {'Content-Type': 'application/json', 'Accept': 'application/json'}

CATEGORY_MAP = {
    'missense_variant':          'Missense',
    'synonymous_variant':        'Synonymous',
    'stop_gained':               'Stop gained',
    'stop_lost':                 'Stop lost',
    'start_lost':                'Start lost',
    'intron_variant':            'Intronic',
    'intergenic_variant':        'Intergenic',
    '5_prime_UTR_variant':       "5' UTR",
    '3_prime_UTR_variant':       "3' UTR",
    'upstream_gene_variant':     'Upstream',
    'downstream_gene_variant':   'Downstream',
    'splice_region_variant':     'Splice region',
    'splice_donor_variant':      'Splice donor/acceptor',
    'splice_acceptor_variant':   'Splice donor/acceptor',
    'non_coding_transcript_exon_variant': 'Non-coding exon',
}

CATEGORY_COLORS = {
    'Missense':               '#e41a1c',
    'Synonymous':             '#984ea3',
    'Stop gained':            '#b10026',
    'Stop lost':              '#b10026',
    'Start lost':             '#b10026',
    'Intronic':               '#4daf4a',
    'Intergenic':             '#999999',
    "5' UTR":                 '#377eb8',
    "3' UTR":                 '#f781bf',
    'Upstream':               '#ff7f00',
    'Downstream':             '#a65628',
    'Splice region':          '#e6ab02',
    'Splice donor/acceptor':  '#d95f02',
    'Non-coding exon':        '#66c2a5',
    'Other':                  '#cccccc',
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
    most_severe = variant.get('most_severe_consequence', '')
    transcripts = variant.get('transcript_consequences', [])
    best = None
    for tc in transcripts:
        if most_severe in tc.get('consequence_terms', []):
            if best is None or tc.get('biotype') == 'protein_coding':
                best = tc
    if best is None and transcripts:
        best = transcripts[0]
    return {
        'most_severe_consequence': most_severe,
        'gene_symbol': best.get('gene_symbol', '') if best else '',
        'biotype':     best.get('biotype',     '') if best else '',
        'category':    CATEGORY_MAP.get(most_severe, 'Other'),
    }


# ── 1. Load unique SNPs ───────────────────────────────────────────────────────
cands = pd.read_csv(CANDIDATES_FILE, sep='\t')
snps = cands.drop_duplicates(subset='ID').copy()
print(f'Unique SNPs to annotate: {len(snps)}')

# ── 2. VEP annotation ─────────────────────────────────────────────────────────
rsids = snps['ID'].tolist()
annotations = {}
for i in range(0, len(rsids), BATCH_SIZE):
    batch = rsids[i:i + BATCH_SIZE]
    print(f'Querying VEP {i + 1}–{i + len(batch)}...')
    for v in vep_batch(batch):
        rid = v.get('id', v.get('input', ''))
        if rid:
            annotations[rid] = parse_vep(v)
    time.sleep(1)

ann_df = pd.DataFrame([
    {'ID': rid, **annotations.get(rid, {
        'most_severe_consequence': 'not_found',
        'gene_symbol': '',
        'biotype':     '',
        'category':    'Other',
    })}
    for rid in rsids
])

snps = snps.merge(ann_df, on='ID', how='left')
snps.to_csv(ANNOT_FILE, sep='\t', index=False)
print(f'\nSaved annotated SNPs → {ANNOT_FILE}')
print(snps[['ID', 'CHROM', 'POS', 'most_severe_consequence', 'category', 'gene_symbol', 'biotype']].to_string(index=False))

# ── 3. Bar plot ───────────────────────────────────────────────────────────────
counts = snps['category'].value_counts().sort_values(ascending=True)
colors = [CATEGORY_COLORS.get(c, CATEGORY_COLORS['Other']) for c in counts.index]

fig, ax = plt.subplots(figsize=(8, max(4, len(counts) * 0.55)))
bars = ax.barh(counts.index, counts.values, color=colors, edgecolor='k', linewidth=0.5)

for bar, val in zip(bars, counts.values):
    ax.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height() / 2,
            str(val), va='center', ha='left', fontsize=10)

ax.set_xlabel('Number of candidate SNPs', fontsize=12)
ax.set_title('Genomic consequence of fine-mapping candidate SNPs\n'
             f'(n={len(snps)} unique SNPs, BY-corrected FDR < 0.05)', fontsize=11)
ax.set_xlim(0, counts.max() * 1.15)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
plt.savefig(PLOT_FILE, dpi=150)
plt.close()
print(f'\nPlot saved → {PLOT_FILE}')
