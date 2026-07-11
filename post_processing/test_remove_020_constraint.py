"""
READ-ONLY test: does removing the alpha=0.20 (BY-FDR) filter change the
admixture-mapping hits?

For every ancestry x phenotype it reproduces the significance pipeline
(BY correction -> FP=1 -> >=5 contiguous windows) in two variants:
  A) WITH the 0.20 filter  (current code: FP=1 on by_p[reject], reject = by_p<=0.20)
  B) WITHOUT the filter     (FP=1 on ALL by_p)
and reports any difference in the significant-window sets and in the final hits.

It ALSO cross-checks variant A against the saved significant_positions.tsv to
confirm the reproduction is faithful.

Writes NOTHING. Reads only P_info_*.tsv and significant_positions.tsv.
"""

import os
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

GEN = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/ukbb'
ANCESTRIES = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']
META = {'#CHROM', 'POS', 'ABS_POS', 'end_POS'}
ALPHA = 0.20


def fp1_threshold(by_p, mode):
    """FP=1 threshold under different candidate-selection rules.
    mode: 'reject020' (current code, by_p<=0.20), 'all' (no filter),
          'remove1' (drop saturated by_p==1 windows only)."""
    if mode == 'reject020':
        vals = np.sort(by_p[by_p <= ALPHA])
    elif mode == 'remove1':
        vals = np.sort(by_p[by_p < 1.0])
    else:  # 'all'
        vals = np.sort(by_p)
    k_max = 0
    for k in range(1, len(vals) + 1):
        if vals[k - 1] * k <= 1:
            k_max = k
        else:
            break
    return vals[k_max - 1] if k_max > 0 else None


def max_contiguous_run(sig):
    """Replicates _max_contiguous_run: adjacency via POS == previous end_POS.
    Vectorised with numpy (fast even for large window sets)."""
    if len(sig) == 0:
        return 0
    s = sig.sort_values(['#CHROM', 'POS'])
    chrom = s['#CHROM'].to_numpy()
    pos   = s['POS'].to_numpy()
    end   = s['end_POS'].to_numpy()
    if len(s) == 1:
        return 1
    contig = (chrom[1:] == chrom[:-1]) & (pos[1:] == end[:-1])
    max_run = cur = 1
    for c in contig:
        cur = cur + 1 if c else 1
        if cur > max_run:
            max_run = cur
    return max_run


MODES = {'A_reject020': 'reject020', 'B_all': 'all', 'C_remove1': 'remove1'}
hits = {name: set() for name in MODES}          # (ancestry, pheno) passing >=5 contiguous

for anc in ANCESTRIES:
    f = os.path.join(GEN, f'P_info_{anc}.tsv')
    if not os.path.exists(f):
        continue
    data = pd.read_csv(f, sep='\t')
    phenos = [c for c in data.columns if c not in META]

    for pheno in phenos:
        valid = data[pheno].notna() & (data[pheno] > 0)
        sub = data.loc[valid, ['#CHROM', 'POS', 'end_POS', pheno]].copy()
        if sub.empty:
            continue
        _, by_p, _, _ = multipletests(sub[pheno].values, alpha=ALPHA, method='fdr_by')
        sub['by_p'] = by_p

        for name, mode in MODES.items():
            thr = fp1_threshold(by_p, mode)
            sig = sub[sub['by_p'] <= thr] if thr is not None else sub.iloc[0:0]
            if max_contiguous_run(sig) >= 5:
                hits[name].add((anc, pheno))

print('=' * 70)
print('TEST: effect of the candidate-selection rule on AM hits')
print('=' * 70)
print('  A_reject020 : current code (by_p <= 0.20 gate)')
print('  B_all       : no gate (FP=1 on all by_p)')
print('  C_remove1   : drop only saturated by_p == 1 windows\n')

for name in MODES:
    print(f'  {name:<12} final hits (>=5 contiguous): {len(hits[name])}')

A = hits['A_reject020']
print(f'\nvs current (A):')
for name in ('B_all', 'C_remove1'):
    extra = sorted(hits[name] - A)
    missing = sorted(A - hits[name])
    print(f'  {name}: identical to A = {hits[name] == A}   (+{len(extra)} extra, -{len(missing)} missing)')
    if 0 < len(extra) <= 25:
        print(f'      extra: {extra}')
hits_A = A  # for the cross-check block below

# ── cross-check variant A against the saved production output ─────────────────
sig_file = os.path.join(GEN, 'significant_positions.tsv')
if os.path.exists(sig_file):
    prod = pd.read_csv(sig_file, sep='\t')
    prod_hits = set(zip(prod['Ancestry'], prod['Phenotype']))
    print(f'\nCross-check vs significant_positions.tsv:')
    print(f'  production hits: {len(prod_hits)}   my variant-A hits: {len(hits_A)}')
    print(f'  A matches production: {hits_A == prod_hits}')
    if hits_A != prod_hits:
        print(f'    only in production: {sorted(prod_hits - hits_A)}')
        print(f'    only in my A:       {sorted(hits_A - prod_hits)}')
