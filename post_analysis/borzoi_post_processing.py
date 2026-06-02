"""
Borzoi post-processing: reads 4-fold sad.h5, exports CSVs, and produces
standard summary plots.

Plots generated:
  1. snp_ranking.png         — SNPs ranked by max |SAD| across all tracks
  2. tracktype_heatmap.png   — Mean |SAD| per SNP × track type (CAGE/RNA/DNASE/ATAC/CHIP)
  3. cage_heatmap.png        — Top 30 CAGE tracks × all SNPs (most interpretable for expression)
  4. top_snps_tracks.png     — Top 10 tracks for the 10 most impacted SNPs

CSVs:
  borzoi_results/f3c{0..3}/sad.csv
  borzoi_results/sad_mean.csv
"""

import os
import re
import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BORZOI_RESULTS = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results'
FOLDS = ['f3c0', 'f3c1', 'f3c2', 'f3c3']
PLOT_DIR = os.path.join(BORZOI_RESULTS, 'plots')
os.makedirs(PLOT_DIR, exist_ok=True)

# ── 1. Load all folds ─────────────────────────────────────────────────────────
sad_arrays   = []
snp_df       = None
track_labels = None
track_ids    = None

for fold in FOLDS:
    h5_path = os.path.join(BORZOI_RESULTS, fold, 'sad.h5')
    print(f'Reading {h5_path} ...')
    with h5py.File(h5_path, 'r') as f:
        sad      = f['SAD'][()].astype(np.float32)
        snp_ids  = [s.decode() for s in f['snp'][()]]
        chroms   = [s.decode() for s in f['chr'][()]]
        pos      = f['pos'][()].tolist()
        refs     = [s.decode() for s in f['ref_allele'][()]]
        alts     = [s.decode() for s in f['alt_allele'][()]]
        t_labels = [s.decode() for s in f['target_labels'][()]]
        t_ids    = [s.decode() for s in f['target_ids'][()]]

    if snp_df is None:
        snp_df       = pd.DataFrame({'snp': snp_ids, 'chr': chroms, 'pos': pos,
                                     'ref': refs, 'alt': alts})
        track_labels = t_labels
        track_ids    = t_ids

    sad_arrays.append(sad)

    # Per-fold CSV
    df_fold = pd.DataFrame(sad, columns=track_labels)
    df_fold = pd.concat([snp_df.reset_index(drop=True), df_fold], axis=1)
    out_csv = os.path.join(BORZOI_RESULTS, fold, 'sad.csv')
    df_fold.to_csv(out_csv, index=False)
    print(f'  → {out_csv}')

# ── 2. Mean across folds ──────────────────────────────────────────────────────
print('\nComputing mean across folds ...')
sad_mean = np.mean(np.stack(sad_arrays, axis=0), axis=0)   # (n_snps, n_tracks)

df_mean = pd.DataFrame(sad_mean, columns=track_labels)
df_mean = pd.concat([snp_df.reset_index(drop=True), df_mean], axis=1)
out_mean = os.path.join(BORZOI_RESULTS, 'sad_mean.csv')
df_mean.to_csv(out_mean, index=False)
print(f'  → {out_mean}')

snp_labels = snp_df['snp'].tolist()
n_snps, n_tracks = sad_mean.shape

# ── 3. Parse track types ──────────────────────────────────────────────────────
def track_type(label):
    return label.split(':')[0] if ':' in label else 'OTHER'

types = np.array([track_type(l) for l in track_labels])
type_order = ['CAGE', 'RNA', 'DNASE', 'ATAC', 'CHIP']

# ── Plot 1: SNP ranking by max |SAD| ─────────────────────────────────────────
print('\nPlot 1: SNP ranking ...')
max_abs_sad = np.max(np.abs(sad_mean), axis=1)
order       = np.argsort(max_abs_sad)[::-1]

fig, ax = plt.subplots(figsize=(8, max(6, n_snps * 0.22)))
colors = ['#d73027' if max_abs_sad[i] > np.percentile(max_abs_sad, 75) else '#4575b4'
          for i in order]
ax.barh(range(len(order)), max_abs_sad[order], color=colors)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([snp_labels[i] for i in order], fontsize=7)
ax.invert_yaxis()
ax.set_xlabel('Max |SAD| across all tracks (mean of 4 folds)')
ax.set_title('SNP ranking by maximum predicted effect')
ax.axvline(np.median(max_abs_sad), color='gray', linestyle='--', lw=1, label='median')
ax.legend(fontsize=8)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'snp_ranking.png'), dpi=150)
plt.close()
print('  → snp_ranking.png')

# ── Plot 2: Track-type heatmap (mean |SAD| per SNP × type) ───────────────────
print('Plot 2: Track-type heatmap ...')
type_matrix = np.zeros((n_snps, len(type_order)))
for j, tt in enumerate(type_order):
    mask = types == tt
    if mask.sum() > 0:
        type_matrix[:, j] = np.mean(np.abs(sad_mean[:, mask]), axis=1)

df_type = pd.DataFrame(type_matrix, index=snp_labels, columns=type_order)
# Order SNPs by total effect
snp_order = df_type.sum(axis=1).sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(7, max(6, n_snps * 0.22)))
sns.heatmap(df_type.loc[snp_order], cmap='YlOrRd', ax=ax,
            linewidths=0.3, cbar_kws={'label': 'Mean |SAD|'})
ax.set_title('Mean |SAD| per SNP by track type\n(mean of 4 folds)')
ax.set_xlabel('Track type')
ax.set_ylabel('SNP')
ax.tick_params(axis='y', labelsize=7)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'tracktype_heatmap.png'), dpi=150)
plt.close()
print('  → tracktype_heatmap.png')

# ── Plot 3: Top CAGE tracks heatmap ──────────────────────────────────────────
print('Plot 3: CAGE heatmap ...')
cage_mask = types == 'CAGE'
sad_cage  = sad_mean[:, cage_mask]
cage_lbls = [track_labels[i] for i in range(n_tracks) if types[i] == 'CAGE']
cage_lbls_short = [re.sub(r'^CAGE:', '', l).split(',')[0][:40] for l in cage_lbls]

# Select top 30 CAGE tracks by variance across SNPs
top30_idx = np.argsort(np.var(sad_cage, axis=0))[::-1][:30]
sad_cage_top = sad_cage[:, top30_idx]
cage_top_lbls = [cage_lbls_short[i] for i in top30_idx]

# Order SNPs by max |SAD| in CAGE tracks
snp_cage_order = np.argsort(np.max(np.abs(sad_cage_top), axis=1))[::-1]

fig, ax = plt.subplots(figsize=(14, max(6, n_snps * 0.22)))
sns.heatmap(sad_cage_top[snp_cage_order, :],
            xticklabels=cage_top_lbls,
            yticklabels=[snp_labels[i] for i in snp_cage_order],
            cmap='RdBu_r', center=0, ax=ax,
            linewidths=0.2, cbar_kws={'label': 'SAD (mean 4 folds)'})
ax.set_title('CAGE tracks: top 30 by variance across SNPs\n(red = ALT increases expression, blue = ALT decreases)')
ax.tick_params(axis='x', labelsize=7, rotation=45)
ax.tick_params(axis='y', labelsize=7)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'cage_heatmap.png'), dpi=150)
plt.close()
print('  → cage_heatmap.png')

# ── Plot 4: Top 10 tracks for the 10 most impacted SNPs ──────────────────────
print('Plot 4: Top tracks for top SNPs ...')
top10_snp_idx = np.argsort(max_abs_sad)[::-1][:10]

fig, axes = plt.subplots(2, 5, figsize=(22, 8))
axes = axes.flatten()

for ax_i, snp_idx in enumerate(top10_snp_idx):
    sad_snp   = sad_mean[snp_idx, :]
    top10_tr  = np.argsort(np.abs(sad_snp))[::-1][:10]
    vals      = sad_snp[top10_tr]
    lbls      = [re.sub(r'^(CAGE|RNA|DNASE|ATAC|CHIP):', '', track_labels[i])[:35]
                 for i in top10_tr]
    colors_bar = ['#d73027' if v > 0 else '#4575b4' for v in vals]

    axes[ax_i].barh(range(10), vals[::-1], color=colors_bar[::-1])
    axes[ax_i].set_yticks(range(10))
    axes[ax_i].set_yticklabels(lbls[::-1], fontsize=7)
    axes[ax_i].axvline(0, color='black', lw=0.5)
    axes[ax_i].set_title(snp_labels[snp_idx], fontsize=9)
    axes[ax_i].set_xlabel('SAD', fontsize=8)

plt.suptitle('Top 10 tracks for the 10 most impacted SNPs\n(red = ALT up, blue = ALT down)',
             fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'top_snps_tracks.png'), dpi=150)
plt.close()
print('  → top_snps_tracks.png')

print(f'\nAll outputs in: {BORZOI_RESULTS}')
print(f'Plots in:       {PLOT_DIR}')
