"""
Borzoi post-processing restricted to the fine-mapping causal candidate SNPs
(one per hit, extracted by fine_mapping_causal_candidate_extraction.py).

Same plot types as borzoi_post_processing.py, but computed only on the
candidate SNPs instead of the full candidate set.

Plots generated (in borzoi_results/plots_causal_candidates/):
  1. snp_ranking.png         — causal SNPs ranked by max |SUM| across all tracks
  2. tracktype_heatmap.png   — Mean |SUM| per SNP × track type (CAGE/RNA/DNASE/ATAC/CHIP)
  3. {type}_heatmap.png       — Top 30 tracks per type × causal SNPs
  4. {type}_tissue_heatmap.png — Tissue aggregation heatmaps (CAGE/RNA/DNASE/ATAC)
  5. snp_tracks.png          — Top 10 tracks for each causal SNP
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

BORZOI_RESULTS  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results'
CANDIDATES_FILE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_conditional_results/fine_mapping_causal_candidates.tsv'
PLOT_DIR        = os.path.join(BORZOI_RESULTS, 'plots_causal_candidates')
os.makedirs(PLOT_DIR, exist_ok=True)

# ── 1. Load causal candidates + SAD mean ──────────────────────────────────────
causal = pd.read_csv(CANDIDATES_FILE, sep='\t')
causal_snps = causal['candidate_snp'].unique().tolist()
print(f'Causal candidate SNPs ({len(causal_snps)}): {causal_snps}')

df_mean = pd.read_csv(os.path.join(BORZOI_RESULTS, 'sad_mean.csv'))
df_mean = df_mean[df_mean['snp'].isin(causal_snps)].reset_index(drop=True)

meta_cols    = ['snp', 'chr', 'pos', 'ref', 'alt']
track_labels = [c for c in df_mean.columns if c not in meta_cols]
sad_mean     = df_mean[track_labels].to_numpy()
snp_labels   = df_mean['snp'].tolist()
n_snps, n_tracks = sad_mean.shape
print(f'Loaded SAD matrix: {n_snps} SNPs x {n_tracks} tracks')

# ── 2. Parse track types ──────────────────────────────────────────────────────
def track_type(label):
    return label.split(':')[0] if ':' in label else 'OTHER'

types = np.array([track_type(l) for l in track_labels])
type_order = ['CAGE', 'RNA', 'DNASE', 'ATAC', 'CHIP']

# ── Plot 1: SNP ranking by max |SUM| ─────────────────────────────────────────
print('\nPlot 1: SNP ranking ...')
max_abs_sad = np.max(np.abs(sad_mean), axis=1)
order       = np.argsort(max_abs_sad)[::-1]

fig, ax = plt.subplots(figsize=(9, max(4, n_snps * 0.5)))
colors = ['#d73027' if max_abs_sad[i] > np.percentile(max_abs_sad, 75) else '#4575b4'
          for i in order]
ax.barh(range(len(order)), max_abs_sad[order], color=colors)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([snp_labels[i] for i in order], fontsize=9)
ax.invert_yaxis()
ax.set_xlabel('Max |SUM| across all tracks (mean of 4 folds)', fontsize=10)
ax.set_title('Causal candidate SNPs ranked by maximum predicted effect', fontsize=11)
ax.axvline(np.median(max_abs_sad), color='gray', linestyle='--', lw=1, label='median')
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'snp_ranking.png'), dpi=150, bbox_inches='tight')
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
snp_order = df_type.sum(axis=1).sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(7, max(4, n_snps * 0.5)))
sns.heatmap(df_type.loc[snp_order], cmap='YlOrRd', ax=ax,
            linewidths=0.3, cbar_kws={'label': 'Mean |SUM|'})
ax.set_title('Mean |SUM| per causal SNP by track type\n(mean of 4 folds)')
ax.set_xlabel('Track type')
ax.set_ylabel('SNP')
ax.tick_params(axis='y', labelsize=9)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'tracktype_heatmap.png'), dpi=150)
plt.close()
print('  → tracktype_heatmap.png')

# ── Plot 3: Per-track-type heatmaps (top 30 tracks each) ─────────────────────
print('Plot 3: Per-track-type heatmaps ...')

prefix_strip = {
    'CAGE':  lambda l: f"[CAGE] {re.sub(r'^CAGE:', '', l).split(',')[0].strip()[:38]}",
    'RNA':   lambda l: f"[RNA] {re.sub(r'^RNA:', '', l)[:38]}",
    'DNASE': lambda l: f"[DNASE] {re.sub(r'^DNASE:', '', l)[:36]}",
    'ATAC':  lambda l: f"[ATAC] {re.sub(r'^ATAC:', '', l)[:37]}",
    'CHIP':  lambda l: f"[CHIP] {re.sub(r'^CHIP:', '', l)[:37]}",
}

for tt in type_order:
    mask = types == tt
    if mask.sum() == 0:
        continue

    sad_tt = sad_mean[:, mask]
    lbl_tt = [prefix_strip[tt](track_labels[i]) for i in range(n_tracks) if types[i] == tt]

    n_top   = min(30, sad_tt.shape[1])
    top_idx = np.argsort(np.var(sad_tt, axis=0))[::-1][:n_top]
    sad_top = sad_tt[:, top_idx]
    top_lbls = [lbl_tt[i] for i in top_idx]

    snp_order_tt = np.argsort(np.max(np.abs(sad_top), axis=1))[::-1]

    fig, ax = plt.subplots(figsize=(18, max(4, n_snps * 0.5)))
    sns.heatmap(sad_top[snp_order_tt, :],
                xticklabels=top_lbls,
                yticklabels=[snp_labels[i] for i in snp_order_tt],
                cmap='RdBu_r', center=0, ax=ax,
                linewidths=0.3, cbar_kws={'label': 'SUM (mean 4 folds)'})
    ax.set_title(f'{tt} tracks: top {n_top} by variance across causal SNPs\n(red = ALT up, blue = ALT down)',
                 fontsize=11)
    ax.tick_params(axis='x', labelsize=7, rotation=90)
    ax.tick_params(axis='y', labelsize=9)
    plt.tight_layout()
    fname = f'{tt.lower()}_heatmap.png'
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {fname}')

# ── Plot 4: Tissue aggregation heatmaps (CAGE/RNA/DNASE/ATAC) ────────────────
print('Plot 4: Tissue heatmaps ...')

def tissue_heatmap_generic(sad_tt, lbls_tt, extractor, title, fname,
                            exclude_patterns=None, n_top=40):
    if exclude_patterns:
        keep = [i for i, l in enumerate(lbls_tt)
                if not any(re.search(p, l.lower()) for p in exclude_patterns)]
        sad_tt = sad_tt[:, keep]
        lbls_tt = [lbls_tt[i] for i in keep]

    if sad_tt.shape[1] == 0:
        print(f'  → {fname} skipped (no tracks after filtering)')
        return

    tissues = [extractor(l) for l in lbls_tt]
    unique  = sorted(set(t for t in tissues if t))

    mat = np.zeros((n_snps, len(unique)))
    for j, tissue in enumerate(unique):
        idx = [i for i, t in enumerate(tissues) if t == tissue]
        mat[:, j] = sad_tt[:, idx].mean(axis=1)

    n_show  = min(n_top, len(unique))
    top_idx = np.argsort(np.max(np.abs(mat), axis=0))[::-1][:n_show]
    mat_top = mat[:, top_idx]
    top_lbls = [unique[i].title() for i in top_idx]
    snp_ord = np.argsort(np.max(np.abs(mat_top), axis=1))[::-1]

    fig, ax = plt.subplots(figsize=(18, max(4, n_snps * 0.5)))
    sns.heatmap(mat_top[snp_ord, :],
                xticklabels=top_lbls,
                yticklabels=[snp_labels[i] for i in snp_ord],
                cmap='RdBu_r', center=0, ax=ax,
                linewidths=0.3, cbar_kws={'label': 'Mean SUM (mean 4 folds)'})
    ax.set_title(f'{title}\n(red = ALT up, blue = ALT down)', fontsize=11)
    ax.tick_params(axis='x', labelsize=7, rotation=90)
    ax.tick_params(axis='y', labelsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {fname}  ({len(unique)} tissues found, showing top {n_show})')

# CAGE
CAGE_EXCLUDE = [
    r'cell line:', r'universal rna', r'reference total rna', r'encode',
    r'mock treated', r'treated with', r'infection', r'pool\d',
    r'^cd\d+[+-]', r'pluriselect', r'expanded$',
    r'carcinoma', r'leukemia', r'lymphoma', r'sarcoma', r'melanoma',
    r'adenocarcinoma', r'blastoma', r'teratoma', r'mesothelioma',
    r'glioblastoma', r'glioma', r'tumor', r'cancer cell',
]

def cage_tissue(label):
    body = re.sub(r'^CAGE:', '', label)
    tissue = body.split(',')[0].strip()
    tissue = re.sub(r'\s+(adult|fetal|newborn|child|embryo|pool\d*).*$', '', tissue, flags=re.IGNORECASE)
    tissue = re.sub(r'\s*-\s*$', '', tissue)
    return tissue.strip().lower()

sad_cage  = sad_mean[:, types == 'CAGE']
lbls_cage = [track_labels[i] for i in range(n_tracks) if types[i] == 'CAGE']
tissue_heatmap_generic(sad_cage, lbls_cage, cage_tissue,
                       'CAGE tissue aggregation: top 40 real tissues by max |SUM|',
                       'cage_tissue_heatmap.png', exclude_patterns=CAGE_EXCLUDE)

# RNA
def rna_tissue(label):
    body = re.sub(r'^RNA:', '', label)
    body = re.sub(r'\s+(male|female)\s+.*$', '', body, flags=re.IGNORECASE)
    body = re.sub(r'\s+(donor|patient|individual).*$', '', body, flags=re.IGNORECASE)
    return body.strip().lower()

RNA_EXCLUDE = [r'cell line:', r'encode', r'mock', r'treated with', r'genetically modified',
               r'carcinoma', r'leukemia', r'lymphoma', r'sarcoma', r'blastoma', r'tumor',
               r'cancer cell', r'adenocarcinoma']

sad_rna  = sad_mean[:, types == 'RNA']
lbls_rna = [track_labels[i] for i in range(n_tracks) if types[i] == 'RNA']
tissue_heatmap_generic(sad_rna, lbls_rna, rna_tissue,
                       'RNA-seq tissue aggregation: top 40 cell/tissue types',
                       'rna_tissue_heatmap.png', exclude_patterns=RNA_EXCLUDE)

# DNASE
def dnase_tissue(label):
    body = re.sub(r'^DNASE:', '', label)
    body = re.sub(r'\s+(male|female)\s+.*$', '', body, flags=re.IGNORECASE)
    body = re.sub(r'\s+(donor|patient).*$', '', body, flags=re.IGNORECASE)
    return body.strip().lower()

DNASE_EXCLUDE = [r'cell line:', r'gm\d+', r'encode', r'mock', r'treated with',
                 r'genetically modified', r'infection', r'carcinoma', r'leukemia',
                 r'lymphoma', r'sarcoma', r'blastoma', r'tumor', r'cancer cell']

sad_dnase  = sad_mean[:, types == 'DNASE']
lbls_dnase = [track_labels[i] for i in range(n_tracks) if types[i] == 'DNASE']
tissue_heatmap_generic(sad_dnase, lbls_dnase, dnase_tissue,
                       'DNASE-seq tissue aggregation: top 40 cell/tissue types',
                       'dnase_tissue_heatmap.png', exclude_patterns=DNASE_EXCLUDE)

# ATAC
def atac_tissue(label):
    body = re.sub(r'^ATAC:', '', label)
    parts = [p.strip() for p in body.split('/')]
    if len(parts) >= 3:
        return f"{parts[1]} – {parts[2]}"
    elif len(parts) == 2:
        return parts[1]
    return body.strip().lower()

sad_atac  = sad_mean[:, types == 'ATAC']
lbls_atac = [track_labels[i] for i in range(n_tracks) if types[i] == 'ATAC']
tissue_heatmap_generic(sad_atac, lbls_atac, atac_tissue,
                       'ATAC-seq tissue/cell-type aggregation: top 40',
                       'atac_tissue_heatmap.png')

# CHIP: format is "CHIP:target:cell_line" — aggregate by target (TF / histone mark)
def chip_target(label):
    body = re.sub(r'^CHIP:', '', label)
    parts = body.split(':')
    return parts[0].strip().lower() if parts else body.strip().lower()

sad_chip  = sad_mean[:, types == 'CHIP']
lbls_chip = [track_labels[i] for i in range(n_tracks) if types[i] == 'CHIP']
tissue_heatmap_generic(sad_chip, lbls_chip, chip_target,
                       'ChIP-seq target aggregation (TF / histone mark): top 40',
                       'chip_target_heatmap.png')

# ── Plot 5: Top 10 tracks for each causal SNP ────────────────────────────────
print('Plot 5: Top tracks for each causal SNP ...')

def format_track_label(label, maxlen=38):
    parts = label.split(':', 1)
    if len(parts) == 2:
        return f'[{parts[0]}] {parts[1][:maxlen]}'
    return label[:maxlen + 8]

n_cols = min(5, n_snps)
n_rows = int(np.ceil(n_snps / n_cols))
fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.4 * n_cols, 8 * n_rows), squeeze=False)
axes = axes.flatten()

for ax_i, snp_idx in enumerate(order):  # order = SNPs sorted by max |SUM| desc
    sad_snp  = sad_mean[snp_idx, :]
    top10_tr = np.argsort(np.abs(sad_snp))[::-1][:10]
    vals     = sad_snp[top10_tr]
    lbls     = [format_track_label(track_labels[i]) for i in top10_tr]
    colors_bar = ['#d73027' if v > 0 else '#4575b4' for v in vals]

    axes[ax_i].barh(range(len(vals)), vals[::-1], color=colors_bar[::-1])
    axes[ax_i].set_yticks(range(len(vals)))
    axes[ax_i].set_yticklabels(lbls[::-1], fontsize=7)
    axes[ax_i].axvline(0, color='black', lw=0.5)
    axes[ax_i].set_title(snp_labels[snp_idx], fontsize=9)
    axes[ax_i].set_xlabel('SUM', fontsize=8)

for ax_i in range(n_snps, len(axes)):
    axes[ax_i].axis('off')

plt.suptitle('Top 10 tracks for each causal candidate SNP\n(red = ALT up, blue = ALT down)',
             fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'snp_tracks.png'), dpi=150)
plt.close()
print('  → snp_tracks.png')

print(f'\nAll outputs in: {PLOT_DIR}')
