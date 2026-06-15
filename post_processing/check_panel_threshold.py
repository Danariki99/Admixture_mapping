"""
Diagnostic check for the panel Manhattan threshold (fast, text-only — no plots).

For each ancestry / mode it compares two thresholds when ALL p-values of ALL
phenotypes are plotted:

  OLD (lenient)  = max p among FP1-significant windows across significant
                   phenotypes. This is the most permissive line and lets points
                   that are NOT FP1-significant (incl. from non-significant
                   phenotypes) appear above it  → apparent hits > real hits.

  NEW (adequate) = strictest line such that EVERY point above it is a genuine
                   FP1-significant window = largest significant p still below the
                   best (smallest) p among all non-significant points. We may drop
                   a few weak true hits, but never draw a false one.

This mirrors the logic now used inside create_panel_manhattan().
"""

import os
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

GENERAL_OUTPUT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/ukbb'
ANCESTRY_LIST = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']
META_COLS = {'#CHROM', 'POS', 'ABS_POS', 'end_POS'}
ALPHA = 0.20


def fp1_threshold_for(by_p, reject):
    by_p_sig = np.sort(by_p[reject])
    k_max = 0
    for k in range(1, len(by_p_sig) + 1):
        if by_p_sig[k - 1] * k <= 1:
            k_max = k
        else:
            break
    return by_p_sig[k_max - 1] if k_max > 0 else None


def check_ancestry(ancestry, suffix):
    plot_anc = 'AMR' if ancestry == 'NAT' else ancestry
    data_file = os.path.join(GENERAL_OUTPUT_FOLDER, f'P_info_{ancestry}.tsv')
    if not os.path.exists(data_file):
        return

    data = pd.read_csv(data_file, sep='\t')
    pheno_cols = [c for c in data.columns if c not in META_COLS]
    if not pheno_cols:
        return

    by_cols, sig_cols, lenient_contrib = {}, {}, {}
    for pheno in pheno_cols:
        valid_idx = data[pheno].notna()
        vals = data.loc[valid_idx, pheno].values
        if len(vals) == 0:
            continue
        reject, by_p, _, _ = multipletests(vals, alpha=ALPHA, method='fdr_by')
        by_series = pd.Series(np.nan, index=data.index)
        by_series.loc[valid_idx] = by_p
        by_cols[pheno] = by_series

        sig_series = pd.Series(False, index=data.index)
        fp1 = fp1_threshold_for(by_p, reject)
        if fp1 is not None:
            sig_series.loc[valid_idx] = by_p <= fp1
            if suffix == 'BY':
                lenient_contrib[pheno] = fp1
            else:
                m = sig_series & by_series.notna()
                lenient_contrib[pheno] = data.loc[m, pheno].max() if m.any() else None
        sig_cols[pheno] = sig_series

    cols = list(by_cols.keys())
    if suffix == 'BY':
        p_matrix = pd.DataFrame(by_cols)[cols]
    else:
        p_matrix = data[cols]
    sig_matrix = pd.DataFrame(sig_cols)[cols]

    pv   = p_matrix.values.flatten()
    sigv = sig_matrix.values.flatten()
    good = ~np.isnan(pv) & (pv > 0)
    pv, sigv = pv[good], sigv[good]
    sig_p, nonsig_p = pv[sigv], pv[~sigv]
    n_sig = int(sig_p.size)

    sig_phenos = {p: t for p, t in lenient_contrib.items() if t is not None}
    if not sig_phenos:
        print(f'\n[{plot_anc} | {suffix}]  no significant phenotypes — no line')
        return

    # OLD lenient threshold
    old_thresh   = max(sig_phenos.values())
    old_apparent = int((pv <= old_thresh).sum())

    # NEW adequate threshold
    new_thresh = None
    if n_sig > 0 and nonsig_p.size > 0:
        floor_nonsig = nonsig_p.min()
        eligible = sig_p[sig_p < floor_nonsig]
        new_thresh = eligible.max() if eligible.size > 0 else None
    elif n_sig > 0:
        new_thresh = sig_p.max()
    new_shown = int((sig_p <= new_thresh).sum()) if new_thresh is not None else 0

    print(f'\n[{plot_anc} | {suffix}]  significant phenotypes={len(sig_phenos)}  '
          f'real FP1 windows={n_sig}')
    print(f'  OLD lenient : p={old_thresh:.4g} (-log10={-np.log10(old_thresh):.2f})  '
          f'apparent hits above={old_apparent}  -> FALSE={old_apparent - n_sig}')
    if new_thresh is not None:
        print(f'  NEW adequate: p={new_thresh:.4g} (-log10={-np.log10(new_thresh):.2f})  '
              f'true hits shown={new_shown}  dropped={n_sig - new_shown}  FALSE=0')
    else:
        print(f'  NEW adequate: none — cannot separate any true hit from the '
              f'best non-significant point')


if __name__ == '__main__':
    for suffix in ('BY', 'raw'):
        print('=' * 72)
        print(f'PANEL THRESHOLD CHECK — suffix = {suffix}')
        print('=' * 72)
        for ancestry in ANCESTRY_LIST:
            check_ancestry(ancestry, suffix)
