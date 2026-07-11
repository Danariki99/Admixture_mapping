"""
EXAMPLE plot of the best Borzoi SED results (illustrative, not the full
production figure set).

Reads the per-(SNP, gene) SED summary produced by borzoi_sed_post_processing.py
(sed_gene_summary.csv: columns snp, gene, gene_name, chr, pos, max_abs_SED,
top_track, mean_SED_alltracks) and draws the top-N SNP-gene pairs by |SED| as a
horizontal bar chart:
    red  bar  -> ALT allele INCREASES predicted expression (SED > 0)
    blue bar  -> ALT allele DECREASES predicted expression (SED < 0)

Usage:
    python borzoi_sed_example_plot.py <result_folder | path_to_sed_gene_summary.csv>

If the SED summary is not found (Borzoi not run yet in this test), it prints a
message and exits 0 so the test pipeline does not fail.
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TOP_N = 20


def find_summary(arg):
    """Locate sed_gene_summary.csv from a direct path or by searching a folder."""
    if arg and os.path.isfile(arg):
        return arg
    roots = [arg] if arg else []
    roots.append(os.getcwd())
    for root in roots:
        if root and os.path.isdir(root):
            hits = glob.glob(os.path.join(root, '**', 'sed_gene_summary.csv'), recursive=True)
            if hits:
                return sorted(hits)[0]
    return None


def main():
    arg = sys.argv[1] if len(sys.argv) > 1 else None
    summary = find_summary(arg)
    if summary is None:
        print('[borzoi_sed_example_plot] sed_gene_summary.csv not found — run the '
              'Borzoi SED steps first (borzoi_sed_*.py). Skipping example plot.')
        return

    df = pd.read_csv(summary)
    if df.empty:
        print(f'[borzoi_sed_example_plot] {summary} is empty. Skipping.')
        return

    df['abs_SED'] = df['max_abs_SED'].abs()
    top = df.sort_values('abs_SED', ascending=False).head(TOP_N).iloc[::-1]

    labels = [f"{r.gene_name}  ({r.snp})" for r in top.itertuples()]
    vals   = top['max_abs_SED'].to_numpy()
    colors = np.where(vals >= 0, '#d62728', '#1f77b4')

    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.32 * len(top))))
    ax.barh(np.arange(len(top)), vals, color=colors, height=0.75)
    ax.axvline(0, color='black', lw=0.6)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel('Borzoi SED (ALT − REF predicted expression, strongest track)', fontsize=9)
    ax.set_title(f'Top {len(top)} Borzoi SED effects — example', fontsize=10, fontweight='bold')
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(facecolor='#d62728', label='ALT increases expression'),
                       Patch(facecolor='#1f77b4', label='ALT decreases expression')],
              fontsize=7, loc='lower right', frameon=False)

    out = os.path.join(os.path.dirname(summary), 'sed_example_top_effects.png')
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'[borzoi_sed_example_plot] saved -> {out}')


if __name__ == '__main__':
    main()
