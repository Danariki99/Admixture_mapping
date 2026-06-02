"""
Fetches GPN-MSA scores for SuSiE credible set SNPs.
Steps:
  1. Load credible_set_annotations.tsv, deduplicate by rsID (max PIP)
  2. Join REF/ALT from SuSiE inputs zscores files
  3. Liftover hg19 → hg38 with pyliftover
  4. Download GPN-MSA parquet files per chromosome from HuggingFace
  5. Lookup scores by (chrom, pos_hg38, ref, alt)
  6. Output TSV + scatter plot PIP vs GPN-MSA score
"""

import os
import re
import certifi
# Fix htslib SSL on conda — must be set before importing pysam
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

ANNOT_FILE   = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/post_analysis/credible_set_annotations.tsv'
SUSIE_INPUTS = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_inputs'
OUTPUT_BASE  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/gpn_msa_results'

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
    'other':                    '#cccccc',
}

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


os.makedirs(OUTPUT_BASE, exist_ok=True)

# ── 1. Load annotations ────────────────────────────────────────────────────────
annot = pd.read_csv(ANNOT_FILE, sep='\t')
annot = annot[~annot['ID'].str.startswith('Affx')].copy()

# Unique SNPs for liftover + score lookup (one row per rsID)
snps = (annot
        .groupby('ID', as_index=False)
        .agg(
            CHROM                   = ('CHROM', 'first'),
            POS                     = ('POS',   'first'),
            max_PIP                 = ('PIP',   'max'),
            most_severe_consequence = ('most_severe_consequence', 'first'),
            gene_symbol             = ('gene_symbol', 'first'),
            biotype                 = ('biotype',     'first'),
        ))
print(f'Unique SNPs to score: {len(snps)}')

# ── 2. Join REF/ALT from SuSiE inputs ─────────────────────────────────────────
ref_alt_rows = []
for hit_dir in sorted(os.listdir(SUSIE_INPUTS)):
    z_file = os.path.join(SUSIE_INPUTS, hit_dir, f'{hit_dir}_zscores.tsv')
    if not os.path.exists(z_file):
        continue
    df = pd.read_csv(z_file, sep='\t', usecols=['ID', 'REF', 'ALT'])
    ref_alt_rows.append(df)

ref_alt = (pd.concat(ref_alt_rows, ignore_index=True)
           .drop_duplicates(subset='ID')
           [['ID', 'REF', 'ALT']])

snps = snps.merge(ref_alt, on='ID', how='left')
missing_alleles = snps['REF'].isna().sum()
if missing_alleles > 0:
    print(f'WARNING: {missing_alleles} SNPs missing REF/ALT — will not get GPN score')

# ── 3. Liftover hg19 → hg38 ───────────────────────────────────────────────────
print('Lifting over hg19 → hg38...')
lo = LiftOver('hg19', 'hg38')

def liftover_pos(chrom, pos):
    result = lo.convert_coordinate(f'chr{chrom}', pos - 1)  # pyliftover is 0-based
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
if n_failed > 0:
    print(f'  Failed: {snps.loc[snps["POS_hg38"].isna(), "ID"].tolist()}')

# ── 4. Download GPN-MSA parquet files per chromosome ─────────────────────────
# ── 4. Open GPN-MSA tabix file (remote, HTTP range requests) ──────────────────
GPN_CACHE_DIR = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/gpn_msa_cache'
os.makedirs(GPN_CACHE_DIR, exist_ok=True)

print(f'\nOpening GPN-MSA tabix file...')
try:
    tbx = pysam.TabixFile(GPN_BGZ_URL)
    print(f'  Opened remotely via HTTP range requests')
except OSError:
    print(f'  Remote open failed (SSL) — downloading via hf_hub_download...')
    bgz_local = hf_hub_download('songlab/gpn-msa-hg38-scores', 'scores.tsv.bgz',
                                  repo_type='dataset', local_dir=GPN_CACHE_DIR)
    _          = hf_hub_download('songlab/gpn-msa-hg38-scores', 'scores.tsv.bgz.tbi',
                                  repo_type='dataset', local_dir=GPN_CACHE_DIR)
    tbx = pysam.TabixFile(bgz_local)
    print(f'  Opened local cache: {bgz_local}')

# Detect column format from header or first record
header = list(tbx.header)
print(f'  Header lines: {header}')

# Probe one position to detect chrom format (chr6 vs 6) and column order
probe_chrom6 = snps.dropna(subset=['POS_hg38']).iloc[0]
probe_pos    = int(probe_chrom6['POS_hg38'])
probe_chrom  = str(int(float(probe_chrom6['CHROM_hg38'])))

for prefix in [f'chr{probe_chrom}', probe_chrom]:
    try:
        probe_records = list(tbx.fetch(prefix, probe_pos - 1, probe_pos + 1))
        if probe_records:
            print(f'  Chrom format in file: "{prefix}"')
            print(f'  Sample record: {probe_records[0]}')
            CHROM_PREFIX = 'chr' if prefix.startswith('chr') else ''
            break
    except (ValueError, KeyError):
        continue
else:
    print('  WARNING: could not probe tabix file — check URL or connectivity')
    CHROM_PREFIX = 'chr'

# ── 5. Lookup scores via tabix ────────────────────────────────────────────────
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
        # Try to parse: expect chrom, pos, ref, alt, score (1-based pos)
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

print('\nLooking up GPN-MSA scores via tabix...')
results = snps.apply(
    lambda r: pd.Series(lookup_score_tabix(r, tbx, CHROM_PREFIX),
                        index=['gpn_msa_score', 'lookup_status']),
    axis=1
)
snps = pd.concat([snps, results], axis=1)

found = (snps['lookup_status'].str.startswith('found')).sum()
print(f'  Scores found: {found}/{len(snps)}')
print(snps.groupby('lookup_status')['ID'].count().to_string())

# ── 6. Join scores back onto full annotation table (per hit × CS × SNP) ───────
out_cols_snps = ['ID', 'CHROM', 'POS', 'CHROM_hg38', 'POS_hg38', 'REF', 'ALT',
                 'gpn_msa_score', 'lookup_status', 'max_PIP',
                 'most_severe_consequence', 'gene_symbol', 'biotype']

# Save global TSV (all unique SNPs)
global_tsv = os.path.join(OUTPUT_BASE, 'gpn_msa_scores_all.tsv')
snps[out_cols_snps].to_csv(global_tsv, sep='\t', index=False)
print(f'\nSaved global TSV → {global_tsv}')

# Merge scores back into per-hit annotation table
annot_scored = annot.merge(
    snps[['ID', 'CHROM_hg38', 'POS_hg38', 'REF', 'ALT', 'gpn_msa_score', 'lookup_status']],
    on='ID', how='left'
)

def get_color(csq):
    return CONSEQUENCE_COLORS.get(csq, CONSEQUENCE_COLORS['other'])

# ── 7. Per-hit TSV + scatter plot ─────────────────────────────────────────────
for hit, hit_df in annot_scored.groupby('hit'):
    hit_dir = os.path.join(OUTPUT_BASE, hit)
    os.makedirs(hit_dir, exist_ok=True)

    # TSV
    hit_tsv = os.path.join(hit_dir, f'{hit}_gpn_msa.tsv')
    hit_df.to_csv(hit_tsv, sep='\t', index=False)

    # Scatter plot
    plot_df = hit_df.dropna(subset=['gpn_msa_score']).copy()
    if plot_df.empty:
        print(f'  {hit}: no scores available, skipping plot')
        continue

    colors = plot_df['most_severe_consequence'].apply(get_color)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(plot_df['gpn_msa_score'], plot_df['PIP'],
               c=colors, s=80, edgecolors='k', linewidths=0.4, zorder=3)

    ax.axvline(GPN_THRESHOLD, color='red', linestyle='--', linewidth=1, alpha=0.7)

    for _, row in plot_df.iterrows():
        if row['PIP'] > 0.1 or row['gpn_msa_score'] < GPN_THRESHOLD:
            ax.annotate(row['ID'], (row['gpn_msa_score'], row['PIP']),
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
    ax.set_ylabel('PIP', fontsize=12)
    ax.set_title(f'{hit} — PIP vs GPN-MSA score', fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_path = os.path.join(hit_dir, f'{hit}_gpn_msa_scatter.png')
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f'  {hit}: saved TSV + plot → {hit_dir}/')
