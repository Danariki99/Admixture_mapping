"""
GPN-MSA scores for fine-mapping candidate SNPs.
Input:  fine_mapping_candidates_annotated.tsv (REF/ALT already present)
Output: gpn_msa_finemapping/ — global TSV + per-hit TSV + scatter plots
"""

import os
import re
import certifi
os.environ['CURL_CA_BUNDLE'] = certifi.where()
os.environ['SSL_CERT_FILE']  = certifi.where()
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pyliftover import LiftOver
from huggingface_hub import hf_hub_download
import pysam

ANNOT_FILE  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_candidates_annotated.tsv'
OUTPUT_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/gpn_msa_finemapping'

GPN_BGZ_URL   = 'https://huggingface.co/datasets/songlab/gpn-msa-hg38-scores/resolve/main/scores.tsv.bgz'
GPN_THRESHOLD = -7.0

CONSEQUENCE_COLORS = {
    'missense_variant':         '#e41a1c',
    'synonymous_variant':       '#984ea3',
    'intron_variant':           '#4daf4a',
    'intergenic_variant':       '#999999',
    'upstream_gene_variant':    '#ff7f00',
    'downstream_gene_variant':  '#a65628',
    '3_prime_UTR_variant':      '#f781bf',
    '5_prime_UTR_variant':      '#377eb8',
    'splice_region_variant':    '#e6ab02',
    'splice_donor_variant':     '#d95f02',
    'splice_acceptor_variant':  '#d95f02',
    'stop_gained':              '#b10026',
    'other':                    '#cccccc',
}

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


os.makedirs(OUTPUT_BASE, exist_ok=True)

# ── 1. Load annotated candidates (already unique per SNP) ─────────────────────
snps = pd.read_csv(ANNOT_FILE, sep='\t')
print(f'SNPs to score: {len(snps)}')

missing_alleles = snps['REF'].isna().sum()
if missing_alleles > 0:
    print(f'WARNING: {missing_alleles} SNPs missing REF/ALT')

# ── 2. Liftover hg19 → hg38 ───────────────────────────────────────────────────
print('Lifting over hg19 → hg38...')
lo = LiftOver('hg19', 'hg38')

def liftover_pos(chrom, pos):
    result = lo.convert_coordinate(f'chr{chrom}', pos - 1)
    if result:
        new_chrom, new_pos, strand, _ = result[0]
        return new_chrom.replace('chr', ''), new_pos + 1, strand
    return None, None, None

snps[['CHROM_hg38', 'POS_hg38', 'strand_hg38']] = snps.apply(
    lambda r: pd.Series(liftover_pos(r['CHROM'], r['POS'])),
    axis=1
)

n_failed = snps['POS_hg38'].isna().sum()
print(f'  Liftover successful: {len(snps) - n_failed}/{len(snps)}')

# ── 3. Open GPN-MSA tabix file ────────────────────────────────────────────────
GPN_CACHE_DIR = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/gpn_msa_cache'
os.makedirs(GPN_CACHE_DIR, exist_ok=True)

print('\nOpening GPN-MSA tabix file...')
try:
    tbx = pysam.TabixFile(GPN_BGZ_URL)
    print('  Opened remotely via HTTP range requests')
except OSError:
    print('  Remote open failed — using local cache...')
    bgz_local = hf_hub_download('songlab/gpn-msa-hg38-scores', 'scores.tsv.bgz',
                                  repo_type='dataset', local_dir=GPN_CACHE_DIR)
    _          = hf_hub_download('songlab/gpn-msa-hg38-scores', 'scores.tsv.bgz.tbi',
                                  repo_type='dataset', local_dir=GPN_CACHE_DIR)
    tbx = pysam.TabixFile(bgz_local)
    print(f'  Opened local cache: {bgz_local}')

header = list(tbx.header)
print(f'  Header: {header}')

probe = snps.dropna(subset=['POS_hg38']).iloc[0]
probe_pos   = int(probe['POS_hg38'])
probe_chrom = str(int(float(probe['CHROM_hg38'])))

for prefix in [f'chr{probe_chrom}', probe_chrom]:
    try:
        recs = list(tbx.fetch(prefix, probe_pos - 1, probe_pos + 1))
        if recs:
            print(f'  Chrom format: "{prefix}"  |  sample: {recs[0]}')
            CHROM_PREFIX = 'chr' if prefix.startswith('chr') else ''
            break
    except (ValueError, KeyError):
        continue
else:
    print('  WARNING: could not probe tabix file')
    CHROM_PREFIX = 'chr'

# ── 4. Lookup scores ──────────────────────────────────────────────────────────
def lookup_score_tabix(row, tbx, chrom_prefix):
    chrom = str(int(float(row['CHROM_hg38']))) if pd.notna(row['CHROM_hg38']) else None
    pos   = row['POS_hg38']
    ref   = str(row['REF']) if pd.notna(row['REF']) else None
    alt   = str(row['ALT']) if pd.notna(row['ALT']) else None

    if chrom is None or pd.isna(pos) or ref is None or alt is None:
        return np.nan, 'missing_coords'

    pos = int(pos)
    region_chrom = f'{chrom_prefix}{chrom}'

    try:
        records = list(tbx.fetch(region_chrom, pos - 1, pos))
    except (ValueError, KeyError):
        return np.nan, 'chrom_not_in_tabix'

    for rec in records:
        fields = rec.split('\t')
        try:
            r_pos = int(fields[1])
            r_ref = fields[2]
            r_alt = fields[3]
            score = float(fields[4])
        except (IndexError, ValueError):
            continue
        if r_pos == pos:
            if r_ref == ref and r_alt == alt:
                return score, 'found'
            if r_ref == revcomp(ref) and r_alt == revcomp(alt):
                return score, 'found_complement'

    return np.nan, 'not_in_dataset'

print('\nLooking up GPN-MSA scores...')
results = snps.apply(
    lambda r: pd.Series(lookup_score_tabix(r, tbx, CHROM_PREFIX),
                        index=['gpn_msa_score', 'lookup_status']),
    axis=1
)
snps = pd.concat([snps, results], axis=1)

found = snps['lookup_status'].str.startswith('found').sum()
print(f'  Scores found: {found}/{len(snps)}')
print(snps.groupby('lookup_status')['ID'].count().to_string())

# ── 5. Save global TSV ────────────────────────────────────────────────────────
out_cols = ['ID', 'CHROM', 'POS', 'CHROM_hg38', 'POS_hg38', 'REF', 'ALT',
            'OR', 'beta', 'P', 'P_BY', 'hit',
            'gpn_msa_score', 'lookup_status',
            'most_severe_consequence', 'gene_symbol', 'biotype']
global_tsv = os.path.join(OUTPUT_BASE, 'gpn_msa_finemapping_all.tsv')
snps[out_cols].to_csv(global_tsv, sep='\t', index=False)
print(f'\nSaved global TSV → {global_tsv}')

# ── 6. Per-hit TSV + scatter plot (GPN score vs -log10 raw P) ─────────────────
def get_color(csq):
    return CONSEQUENCE_COLORS.get(csq, CONSEQUENCE_COLORS['other'])

for hit, hit_df in snps.groupby('hit'):
    hit_dir = os.path.join(OUTPUT_BASE, hit)
    os.makedirs(hit_dir, exist_ok=True)

    hit_tsv = os.path.join(hit_dir, f'{hit}_gpn_msa.tsv')
    hit_df.to_csv(hit_tsv, sep='\t', index=False)

    plot_df = hit_df.dropna(subset=['gpn_msa_score']).copy()
    if plot_df.empty:
        print(f'  {hit}: no scores, skipping plot')
        continue

    # y-axis: -log10(P), cap zeros at 1e-300
    plot_df['neglog10P'] = -np.log10(plot_df['P'].clip(lower=1e-300))
    colors = plot_df['most_severe_consequence'].apply(get_color)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(plot_df['gpn_msa_score'], plot_df['neglog10P'],
               c=colors, s=80, edgecolors='k', linewidths=0.4, zorder=3)

    ax.axvline(GPN_THRESHOLD, color='red', linestyle='--', linewidth=1, alpha=0.7)

    for _, row in plot_df.iterrows():
        if row['neglog10P'] > 10 or row['gpn_msa_score'] < GPN_THRESHOLD:
            ax.annotate(row['ID'], (row['gpn_msa_score'], row['neglog10P']),
                        fontsize=7, xytext=(4, 3), textcoords='offset points')

    seen_csq = plot_df['most_severe_consequence'].unique()
    handles = [
        plt.Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=CONSEQUENCE_COLORS.get(c, CONSEQUENCE_COLORS['other']),
                   markeredgecolor='k', markersize=8, label=c)
        for c in sorted(seen_csq)
    ]
    handles.append(plt.Line2D([0], [0], color='red', linestyle='--',
                               label=f'GPN threshold ({GPN_THRESHOLD})'))
    ax.legend(handles=handles, fontsize=8, loc='upper left', framealpha=0.8)

    ax.set_xlabel('GPN-MSA score (LLR)', fontsize=12)
    ax.set_ylabel('-log10(P)', fontsize=12)
    ax.set_title(f'{hit} — association strength vs evolutionary conservation', fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(hit_dir, f'{hit}_gpn_msa_scatter.png')
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f'  {hit}: saved → {hit_dir}/')
