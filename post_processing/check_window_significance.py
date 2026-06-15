"""
Check why a given (ancestry, phenotype, window) result was/wasn't called a hit
by the main post-processing pipeline (post_processing_functions.result_analysis).

Reproduces, for a single (ancestry, phenotype) pair:
  1. the lambda_GC QC filter (must be in [0.9, 1.1] to even be considered),
  2. the BY (fdr_by, alpha=0.20) correction over all genome-wide windows,
  3. the FP<=1 threshold on BY-significant windows,
  4. the >=5-contiguous-windows requirement,
and reports where the window of interest stands relative to all of these.

Usage:
  python check_window_significance.py --chrom 6 --pos 32383050 \
      --phenos HC1581 HC382 --ancestries EUR SAS
"""

import argparse
import os
import re
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

DATASET = 'ukbb'
BASE = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes'
GLM_TEMPLATE = BASE + f'/output/{DATASET}/output_ancestry_{{ancestry}}/{{pheno}}/output.{{pheno}}.glm.logistic.hybrid'
LOG_TEMPLATE = BASE + f'/output/{DATASET}/output_ancestry_{{ancestry}}/{{pheno}}/output.log'
WINDOW_POS_FILE = BASE + f'/post_processing_files/{DATASET}/positions.csv'

LAMBDA_MIN, LAMBDA_MAX = 0.9, 1.1
ALPHA = 0.20


def parse_lambda(log_path):
    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        for line in f:
            m = re.search(r'lambda \(based on median chisq\) = ([\d.]+?)\.?\s', line)
            if m:
                return float(m.group(1))
    return None


def max_contiguous_run(sig_df):
    s = sig_df.sort_values(['#CHROM', 'POS']).reset_index(drop=True)
    max_run = cur = 1
    for i in range(1, len(s)):
        if s.loc[i, '#CHROM'] == s.loc[i - 1, '#CHROM'] and s.loc[i, 'POS'] == s.loc[i - 1, 'end_POS']:
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 1
    return max_run


def analyze(ancestry, pheno, window_pos, target_chrom, target_pos):
    glm_file = GLM_TEMPLATE.format(ancestry=ancestry, pheno=pheno)
    log_file = LOG_TEMPLATE.format(ancestry=ancestry, pheno=pheno)

    print(f'\n=== {ancestry} / {pheno} ===')

    if not os.path.exists(glm_file):
        print(f'  file not found: {glm_file}')
        return

    lam = parse_lambda(log_file)
    print(f'  lambda_GC = {lam}')
    lam_ok = lam is not None and (LAMBDA_MIN <= lam <= LAMBDA_MAX)
    if not lam_ok:
        print(f'  -> FAILS lambda_GC QC ({LAMBDA_MIN}-{LAMBDA_MAX}): excluded from pipeline entirely')
        print(f'     (not present in P_info_{ancestry}.tsv, no BY correction ever computed)')

    df = pd.read_table(glm_file, sep='\t')
    df = df[['#CHROM', 'POS', 'P']]
    df = df[df['P'] != '.'].copy()
    df['P'] = pd.to_numeric(df['P'], errors='coerce')
    df = df.dropna(subset=['P'])
    df = pd.merge(df, window_pos, on=['#CHROM', 'POS'], how='inner')

    # BY correction across all genome-wide windows for this (ancestry, pheno)
    reject, by_p, _, _ = multipletests(df['P'].values, alpha=ALPHA, method='fdr_by')
    df['BY_P'] = by_p

    by_p_sig = np.sort(by_p[reject])
    k_max = 0
    for k in range(1, len(by_p_sig) + 1):
        if by_p_sig[k - 1] * k <= 1:
            k_max = k
        else:
            break

    if k_max > 0:
        fp1_threshold = by_p_sig[k_max - 1]
        sig = df[df['BY_P'] <= fp1_threshold].copy()
        run = max_contiguous_run(sig) if not sig.empty else 0
        print(f'  n windows total           = {len(df)}')
        print(f'  n BY-significant (<=0.20) = {reject.sum()}')
        print(f'  FP<=1 BY threshold        = {fp1_threshold:.4g}  '
              f'(n_sig={len(sig)}, expected FP={len(sig) * fp1_threshold:.2f})')
        print(f'  max contiguous run        = {run}  (need >=5 to be a "hit")')
        is_hit = run >= 5 and lam_ok
    else:
        fp1_threshold = None
        print('  no BY-significant windows at all (k_max=0)')
        is_hit = False

    # window of interest
    row = df[(df['#CHROM'] == target_chrom) & (df['POS'] == target_pos)]
    if row.empty:
        print(f'  window chr{target_chrom}:{target_pos} not found in this glm file')
        return
    row = row.iloc[0]
    raw_p, by_pv = row['P'], row['BY_P']
    print(f'  window chr{target_chrom}:{target_pos}-{int(row["end_POS"])}')
    print(f'    raw P = {raw_p:.4g}   BY P = {by_pv:.4g}')

    if not lam_ok:
        verdict = 'EXCLUDED upstream by lambda_GC filter (never reaches BY step)'
    elif fp1_threshold is None:
        verdict = 'no FP<=1 threshold exists -> not a hit'
    elif by_pv > fp1_threshold:
        verdict = f'BY P ({by_pv:.4g}) > FP<=1 threshold ({fp1_threshold:.4g}) -> not a hit'
    else:
        verdict = 'passes BY+FP<=1 threshold'
        if not is_hit:
            verdict += ', but contiguous-run-of-5 requirement fails -> not a hit'
        else:
            verdict += ' and contiguous-run requirement -> HIT'
    print(f'    => {verdict}')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--chrom', type=int, default=6)
    ap.add_argument('--pos', type=int, default=32383050)
    ap.add_argument('--phenos', nargs='+', default=['HC1581', 'HC382'])
    ap.add_argument('--ancestries', nargs='+', default=['EUR', 'SAS'])
    args = ap.parse_args()

    window_pos = pd.read_csv(WINDOW_POS_FILE, sep='\t')

    for pheno in args.phenos:
        for ancestry in args.ancestries:
            analyze(ancestry, pheno, window_pos, args.chrom, args.pos)


if __name__ == '__main__':
    main()
