"""
LocusZoom-like plots for the fine-mapped variants, one per significant hit
(only the regions we caught as significant).

For each hit we plot the fine-mapping ADD -log10(p) across the significant
windows, highlight the candidate SNPs and the lead (causal) SNP, and draw the
significance threshold = the p-value of the LAST significant SNP (the largest
raw p among that hit's candidates).

Figures are named 'Supplementary Figure N' continuing the numbering AFTER the
QQ plots (the last QQ figure number + 1).
"""

import os
import re
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

FINE_MAP_DIR  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
CAND_FILE     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv'
LEAD_FILE     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_causal_candidates.tsv'
QQ_DIR        = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/QQ_plots'
OUTPUT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/locuszoom_plots'
PHENO_TABLE   = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx'


def pheno_label(excel_df, pheno):
    row = excel_df.loc[excel_df['ID'] == pheno, 'ID2']
    name = row.iloc[0] if not row.empty else pheno
    return ' '.join(name.split('_')[1:]) if '_' in name else name


def load_hit_add(hit):
    """Concatenate ADD rows across all fine-mapping windows of a hit."""
    rows = []
    for f in sorted(glob.glob(os.path.join(FINE_MAP_DIR, hit, '*.glm.logistic.hybrid'))):
        df = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
        rows.append(df[df['TEST'] == 'ADD'])
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True)
    d['POS'] = pd.to_numeric(d['POS'], errors='coerce')
    d['P']   = pd.to_numeric(d['P'], errors='coerce')
    d = d.dropna(subset=['POS', 'P'])
    d = d[(d['P'] > 0) & np.isfinite(d['P'])].drop_duplicates(subset='ID')
    return d


def starting_number(qq_dir):
    """Next Supplementary-Figure number = (max QQ number) + 1."""
    nums = []
    for f in glob.glob(os.path.join(qq_dir, '*.png')):
        m = re.search(r'Supplementary[_ ]Figure[_ ](\d+)', os.path.basename(f))
        if m:
            nums.append(int(m.group(1)))
    return (max(nums) + 1) if nums else 1


def plot_locuszoom(d, lead_snp, cand_ids, threshold_p, chrom, title, out_png):
    d = d.copy()
    d['pos_mb']  = d['POS'] / 1e6
    d['mlog10p'] = -np.log10(d['P'])

    is_lead = d['ID'].astype(str) == str(lead_snp)
    is_cand = d['ID'].isin(cand_ids) & ~is_lead
    other   = ~is_lead & ~is_cand

    fig, ax = plt.subplots(figsize=(10.5, 4.2))

    ax.scatter(d.loc[other, 'pos_mb'], d.loc[other, 'mlog10p'],
               s=14, color='#bdbdbd', alpha=0.8, linewidths=0, label='Other SNPs')
    ax.scatter(d.loc[is_cand, 'pos_mb'], d.loc[is_cand, 'mlog10p'],
               s=18, color='#1f77b4', alpha=0.9, linewidths=0, label='Candidate SNPs')
    ax.scatter(d.loc[is_lead, 'pos_mb'], d.loc[is_lead, 'mlog10p'],
               s=42, color='#d62728', edgecolors='black', linewidths=0.8,
               zorder=5, label=f'Lead: {lead_snp}')

    if is_lead.any():
        x0 = float(d.loc[is_lead, 'pos_mb'].iloc[0])
        y0 = float(d.loc[is_lead, 'mlog10p'].iloc[0])
        ax.text(x0, y0, lead_snp, fontsize=9, ha='center', va='bottom')

    # significance threshold = last significant SNP (largest raw p among candidates)
    if threshold_p and threshold_p > 0:
        ax.axhline(-np.log10(threshold_p), color='red', ls='--', lw=1,
                   label='Significance threshold')

    ax.set_title(title)
    ax.set_xlabel(f'Position on {chrom} (Mb)')
    ax.set_ylabel(r'$-\log_{10}(p\mathrm{-value})$')
    ax.legend(frameon=True, fontsize=8, loc='upper left')
    ax.set_xlim(d['pos_mb'].min(), d['pos_mb'].max())
    ax.set_ylim(0, max(1.0, d['mlog10p'].max() * 1.05))
    fig.tight_layout()

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    fig.savefig(out_png, dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    excel_df = pd.read_excel(PHENO_TABLE, sheet_name='first_batch', usecols='B:C')
    cand = pd.read_csv(CAND_FILE, sep='\t')
    lead = pd.read_csv(LEAD_FILE, sep='\t')

    hits = sorted(cand['hit'].unique())          # the 7 fine-mapped (significant) hits
    counter = starting_number(QQ_DIR)
    print(f'Starting at Supplementary Figure {counter} (after {counter - 1} QQ plots)\n')

    for hit in hits:
        d = load_hit_add(hit)
        if d is None or d.empty:
            print(f'  {hit}: no ADD data, skipping')
            continue

        sub         = cand[cand['hit'] == hit]
        cand_ids    = set(sub['ID'])
        threshold_p = sub['P'].max()             # last significant SNP
        lead_row    = lead[lead['hit'] == hit]
        lead_snp    = (lead_row['candidate_snp'].iloc[0] if not lead_row.empty
                       else sub.loc[sub['P'].idxmin(), 'ID'])

        anc, pheno, chrom = hit.split('_')
        anc_disp = 'AMR' if anc == 'NAT' else anc
        pname    = pheno_label(excel_df, pheno)
        title    = f'{anc_disp} | {pname} | {chrom} | lead {lead_snp}'
        out_png  = os.path.join(OUTPUT_FOLDER, f'Supplementary Figure {counter}.png')

        plot_locuszoom(d, lead_snp, cand_ids, threshold_p, chrom, title, out_png)
        print(f'  [{counter}] {hit}  n_snps={len(d)}  candidates={len(cand_ids)}  '
              f'threshold_p={threshold_p:.2e}  lead={lead_snp}')
        counter += 1

    print(f'\nDone. Next free number = {counter}.')
