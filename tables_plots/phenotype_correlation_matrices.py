"""
Phenotype correlation matrices.

Builds a sample x phenotype matrix from the .phe files (case/control coded 2/1,
recoded to 1/0), then draws two figures:

  1. correlation_matrix_all_phenos.png  — correlation between ALL phenotypes used
     (every .phe in the folder), hierarchically clustered.
  2. correlation_matrix_hit_phenos.png  — correlation only between the phenotypes
     that appear in our admixture-mapping hits, annotated.

Correlation is the Pearson correlation between the binary phenotypes (equivalent
to the phi coefficient), computed pairwise on the samples measured for both.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

PHENO_DIR   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb'
# analysis cohort: samples kept after the KING kinship cutoff (unrelated set)
KCUTOFF     = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id'
OUT_DIR     = os.path.dirname(os.path.abspath(__file__))

# phenotype dropped from the paper (TTE_asthma), excluded from the matrix -> 108 phenos
EXCLUDED_PHENOS = {'HC1036'}

# phenotypes that appear in the admixture-mapping hits (one per hit trait)
HIT_PHENOS = ['HC1574', 'HC219', 'HC326', 'HC643', 'HC1158', 'HC964',
              'HC937', 'HC1036', 'HC1581', 'HC382', 'HC164']
HIT_PHENOS = [p for p in HIT_PHENOS if p not in EXCLUDED_PHENOS]

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})


# ── build the sample x phenotype matrix ────────────────────────────────────────
def load_matrix():
    cols = {}
    for f in sorted(glob.glob(os.path.join(PHENO_DIR, '*.phe'))):
        pid = os.path.basename(f)[:-4]
        if pid in EXCLUDED_PHENOS:                              # TTE_asthma dropped
            continue
        s = pd.read_csv(f, sep=r'\s+', dtype={0: str})
        s = s.set_index(s.columns[0]).iloc[:, 0]                # index=IID, values=1/2
        s = pd.to_numeric(s, errors='coerce')
        cols[pid] = s.map({2: 1.0, 1: 0.0})                     # case=1, control=0, else NaN
    mat = pd.DataFrame(cols)
    print(f'Loaded {mat.shape[1]} phenotypes x {mat.shape[0]:,} samples (all)')

    # restrict to the analysis cohort (KING-cutoff, unrelated samples)
    keep = set(pd.read_csv(KCUTOFF, dtype=str).iloc[:, 0].str.strip())
    mat = mat.loc[mat.index.astype(str).isin(keep)]
    print(f'After KING-cutoff filter: {mat.shape[0]:,} samples')
    return mat


def main():
    mat = load_matrix()
    corr = mat.corr(method='pearson', min_periods=1000)        # phi coefficient, pairwise
    hits = [p for p in HIT_PHENOS if p in corr.columns]

    # ── export the matrix + summary statistics (mean correlation etc.) ────────
    corr.to_csv(os.path.join(OUT_DIR, 'phenotype_correlation_matrix.csv'))

    def offdiag(c):
        """Off-diagonal values of a correlation matrix, upper triangle only."""
        v = c.to_numpy(dtype=float)
        iu = np.triu_indices_from(v, k=1)
        x = v[iu]
        return x[np.isfinite(x)]

    rows = []
    for name, c in [(f'all {corr.shape[0]} phenotypes', corr),
                    (f'hit phenotypes ({len(hits)})', corr.loc[hits, hits])]:
        x = offdiag(c)
        rows.append({
            'set': name, 'n_pairs': len(x),
            'mean_r': x.mean(), 'mean_abs_r': np.abs(x).mean(),
            'median_r': np.median(x), 'sd_r': x.std(),
            'min_r': x.min(), 'max_r': x.max(),
            'pct_|r|>0.1': 100 * (np.abs(x) > 0.1).mean(),
            'pct_|r|>0.5': 100 * (np.abs(x) > 0.5).mean(),
            'pct_|r|>0.8': 100 * (np.abs(x) > 0.8).mean(),
        })
    stats = pd.DataFrame(rows)
    stats.to_csv(os.path.join(OUT_DIR, 'phenotype_correlation_summary.csv'), index=False)
    print('\n=== Correlation summary (off-diagonal pairs) ===')
    print(stats.round(4).to_string(index=False))

    # mean correlation excluding the redundant re-codings of the same trait
    x_hit = offdiag(corr.loc[hits, hits])
    print(f"\nHit phenotypes: mean r = {x_hit.mean():.4f} | excluding duplicate codings "
          f"(|r|>0.8): mean r = {x_hit[np.abs(x_hit) <= 0.8].mean():.4f}")
    print()

    # ── 1) all phenotypes — plain heatmap (NO clustering, just the matrix) ─────
    lab_all = list(corr.columns)          # HC#### codes only
    fig, ax = plt.subplots(figsize=(24, 22))
    sns.heatmap(corr.fillna(0.0), cmap='RdBu_r', center=0, vmin=-1, vmax=1, square=True,
                xticklabels=lab_all, yticklabels=lab_all, linewidths=0,
                cbar_kws={'label': 'Pearson r (phi)', 'shrink': 0.4,
                          'ticks': [-1, -0.5, 0, 0.5, 1]}, ax=ax)
    ax.tick_params(labelsize=6)
    ax.set_title(f'Phenotype correlation — all {corr.shape[0]} phenotypes  (KING-cutoff cohort)',
                 fontsize=16, fontweight='bold')
    out1 = os.path.join(OUT_DIR, 'correlation_matrix_all_phenos.png')
    fig.tight_layout()
    fig.savefig(out1, dpi=300, bbox_inches='tight')
    fig.savefig(out1.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f'Saved -> {out1}')

    # ── 2) hit phenotypes only — annotated heatmap ────────────────────────────
    hits = [p for p in HIT_PHENOS if p in corr.columns]
    sub = corr.loc[hits, hits]
    lab_hit = list(hits)                  # HC#### codes only
    fig, ax = plt.subplots(figsize=(9, 7.5))
    sns.heatmap(sub, annot=True, fmt='.2f', cmap='RdBu_r', center=0, vmin=-1, vmax=1,
                square=True, linewidths=0.5, linecolor='white',
                xticklabels=lab_hit, yticklabels=lab_hit,
                cbar_kws={'label': 'Pearson r (phi)'}, annot_kws={'size': 7}, ax=ax)
    ax.tick_params(labelsize=7)
    ax.set_title('Phenotype correlation — admixture-mapping hit phenotypes',
                 fontsize=11, fontweight='bold')
    fig.tight_layout()
    out2 = os.path.join(OUT_DIR, 'correlation_matrix_hit_phenos.png')
    fig.savefig(out2, dpi=300, bbox_inches='tight')
    fig.savefig(out2.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close(fig)
    print(f'Saved -> {out2}')


if __name__ == '__main__':
    main()
