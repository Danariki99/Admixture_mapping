import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

DATASET    = 'ukbb'
TSV_FOLDER = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/{DATASET}'
ANCESTRIES = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']
META_COLS  = {'#CHROM', 'POS', 'end_POS', 'ABS_POS'}

total_hits = 0

for ancestry in ANCESTRIES:
    tsv = os.path.join(TSV_FOLDER, f'P_info_{ancestry}.tsv')
    if not os.path.exists(tsv):
        print(f'[{ancestry}] TSV not found, skipping')
        continue

    data = pd.read_csv(tsv, sep='\t')
    pheno_cols = [c for c in data.columns if c not in META_COLS]

    ancestry_hits = 0
    for pheno in pheno_cols:
        p = pd.to_numeric(data[pheno], errors='coerce')
        valid = p.notna() & (p > 0)
        if valid.sum() == 0:
            continue

        _, corrected, _, _ = multipletests(p[valid].values, alpha=0.05, method='fdr_by')
        sig_mask = corrected < 0.05

        if sig_mask.any():
            p_valid    = p[valid]
            p_sig      = p_valid[sig_mask]
            best_idx   = p_sig.idxmin()
            best_row   = data.loc[best_idx]
            best_raw   = p_sig.min()
            best_by    = corrected[sig_mask].min()
            n_sig      = sig_mask.sum()
            emp_thresh = p_valid[sig_mask].max()

            print(
                f'  HIT: {ancestry} {pheno} | '
                f'chr{int(best_row["#CHROM"])}:{int(best_row["POS"])} | '
                f'raw_p={best_raw:.2e} | BY_p={best_by:.2e} | '
                f'n_sig_windows={n_sig} | empirical_thresh={emp_thresh:.2e}'
            )
            ancestry_hits += 1

    if ancestry_hits == 0:
        print(f'[{ancestry}] No significant hits')
    else:
        print(f'[{ancestry}] {ancestry_hits} significant phenotype(s)')
        total_hits += ancestry_hits

print(f'\nTotal: {total_hits} significant (ancestry, pheno) pairs across all ancestries')
