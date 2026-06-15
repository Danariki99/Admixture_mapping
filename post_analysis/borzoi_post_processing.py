"""
Borzoi post-processing: reads 4-fold sad.h5, exports CSVs, and produces
standard summary plots.

Plots generated:
  1. snp_ranking.png         — SNPs ranked by max |SUM| across all tracks
  2. tracktype_heatmap.png   — Mean |SUM| per SNP × track type (CAGE/RNA/DNASE/ATAC/CHIP)
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

# ── Plot 1: SNP ranking by max |SUM| ─────────────────────────────────────────
print('\nPlot 1: SNP ranking ...')
max_abs_sad = np.max(np.abs(sad_mean), axis=1)
order       = np.argsort(max_abs_sad)[::-1]

fig, ax = plt.subplots(figsize=(9, max(7, n_snps * 0.28)))
colors = ['#d73027' if max_abs_sad[i] > np.percentile(max_abs_sad, 75) else '#4575b4'
          for i in order]
ax.barh(range(len(order)), max_abs_sad[order], color=colors)
ax.set_yticks(range(len(order)))
ax.set_yticklabels([snp_labels[i] for i in order], fontsize=8)
ax.invert_yaxis()
ax.set_xlabel('Max |SUM| across all tracks (mean of 4 folds)', fontsize=10)
ax.set_title('SNP ranking by maximum predicted effect', fontsize=11)
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
# Order SNPs by total effect
snp_order = df_type.sum(axis=1).sort_values(ascending=False).index

fig, ax = plt.subplots(figsize=(7, max(6, n_snps * 0.22)))
sns.heatmap(df_type.loc[snp_order], cmap='YlOrRd', ax=ax,
            linewidths=0.3, cbar_kws={'label': 'Mean |SUM|'})
ax.set_title('Mean |SUM| per SNP by track type\n(mean of 4 folds)')
ax.set_xlabel('Track type')
ax.set_ylabel('SNP')
ax.tick_params(axis='y', labelsize=7)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'tracktype_heatmap.png'), dpi=150)
plt.close()
print('  → tracktype_heatmap.png')

# ── Plot 3: Per-track-type heatmaps (top 30 tracks each) ─────────────────────
print('Plot 3: Per-track-type heatmaps ...')

# Label formatters per type — keep type prefix for readability
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

    sad_tt   = sad_mean[:, mask]
    lbl_tt   = [prefix_strip[tt](track_labels[i]) for i in range(n_tracks) if types[i] == tt]

    n_top    = min(30, sad_tt.shape[1])
    top_idx  = np.argsort(np.var(sad_tt, axis=0))[::-1][:n_top]
    sad_top  = sad_tt[:, top_idx]
    top_lbls = [lbl_tt[i] for i in top_idx]

    snp_order_tt = np.argsort(np.max(np.abs(sad_top), axis=1))[::-1]

    fig, ax = plt.subplots(figsize=(18, max(8, n_snps * 0.28)))
    sns.heatmap(sad_top[snp_order_tt, :],
                xticklabels=top_lbls,
                yticklabels=[snp_labels[i] for i in snp_order_tt],
                cmap='RdBu_r', center=0, ax=ax,
                linewidths=0.3, cbar_kws={'label': 'SUM (mean 4 folds)'})
    ax.set_title(f'{tt} tracks: top {n_top} by variance across SNPs\n(red = ALT up, blue = ALT down)',
                 fontsize=11)
    ax.tick_params(axis='x', labelsize=7, rotation=90)
    ax.tick_params(axis='y', labelsize=8)
    plt.tight_layout()
    fname = f'{tt.lower()}_heatmap.png'
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {fname}')

# ── Plot 4: CAGE tissue heatmap (real tissues only) ──────────────────────────
print('Plot 4: CAGE tissue heatmap ...')
cage_mask = types == 'CAGE'
sad_cage  = sad_mean[:, cage_mask]
cage_lbls = [track_labels[i] for i in range(n_tracks) if types[i] == 'CAGE']

# Exclude cell lines, technical references, and highly specific immune subtypes
EXCLUDE_PATTERNS = [
    r'cell line:', r'universal rna', r'reference total rna', r'encode',
    r'mock treated', r'treated with', r'infection', r'pool\d',
    r'^cd\d+[+-]', r'pluriselect', r'expanded$',
    r'carcinoma', r'leukemia', r'lymphoma', r'sarcoma', r'melanoma',
    r'adenocarcinoma', r'blastoma', r'teratoma', r'mesothelioma',
    r'glioblastoma', r'glioma', r'tumor', r'cancer cell',
]

def is_real_tissue(label):
    body = re.sub(r'^CAGE:', '', label).lower()
    return not any(re.search(p, body) for p in EXCLUDE_PATTERNS)

def cage_tissue(label):
    body = re.sub(r'^CAGE:', '', label)
    tissue = body.split(',')[0].strip()
    tissue = re.sub(r'\s+(adult|fetal|newborn|child|embryo|pool\d*).*$', '', tissue, flags=re.IGNORECASE)
    tissue = re.sub(r'\s*-\s*$', '', tissue)  # trailing dash
    return tissue.strip().lower()

# Filter to real tissues only
real_idx    = [i for i, l in enumerate(cage_lbls) if is_real_tissue(l)]
sad_real    = sad_cage[:, real_idx]
real_tissues = [cage_tissue(cage_lbls[i]) for i in real_idx]
unique_tissues = sorted(set(real_tissues))

# Aggregate: mean SAD per tissue per SNP
tissue_matrix = np.zeros((n_snps, len(unique_tissues)))
for j, tissue in enumerate(unique_tissues):
    idx = [i for i, t in enumerate(real_tissues) if t == tissue]
    tissue_matrix[:, j] = sad_real[:, idx].mean(axis=1)

# Top 40 tissues by max |SUM|
tissue_signal  = np.max(np.abs(tissue_matrix), axis=0)
top_tissue_idx = np.argsort(tissue_signal)[::-1][:40]
tissue_top     = tissue_matrix[:, top_tissue_idx]
top_names      = [unique_tissues[i].title() for i in top_tissue_idx]

snp_tissue_order = np.argsort(np.max(np.abs(tissue_top), axis=1))[::-1]

fig, ax = plt.subplots(figsize=(18, max(8, n_snps * 0.25)))
sns.heatmap(tissue_top[snp_tissue_order, :],
            xticklabels=top_names,
            yticklabels=[snp_labels[i] for i in snp_tissue_order],
            cmap='RdBu_r', center=0, ax=ax,
            linewidths=0.3, cbar_kws={'label': 'Mean CAGE SUM (mean 4 folds)'})
ax.set_title('CAGE tissue aggregation: top 40 real tissues by max |SUM|\n(red = ALT ↑ transcription, blue = ALT ↓ transcription)',
             fontsize=11)
ax.tick_params(axis='x', labelsize=7, rotation=90)
ax.tick_params(axis='y', labelsize=7)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'cage_tissue_heatmap.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f'  → cage_tissue_heatmap.png  ({len(unique_tissues)} real tissues found, showing top 40)')

# ── Plot 4b/4c/4d: Tissue heatmaps for RNA, DNASE, ATAC ──────────────────────
def tissue_heatmap_generic(sad_tt, lbls_tt, extractor, title, fname,
                            exclude_patterns=None, n_top=40):
    """Generic tissue aggregation heatmap for any track type."""
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

    n_show   = min(n_top, len(unique))
    top_idx  = np.argsort(np.max(np.abs(mat), axis=0))[::-1][:n_show]
    mat_top  = mat[:, top_idx]
    top_lbls = [unique[i].title() for i in top_idx]
    snp_ord  = np.argsort(np.max(np.abs(mat_top), axis=1))[::-1]

    fig, ax = plt.subplots(figsize=(18, max(8, n_snps * 0.25)))
    sns.heatmap(mat_top[snp_ord, :],
                xticklabels=top_lbls,
                yticklabels=[snp_labels[i] for i in snp_ord],
                cmap='RdBu_r', center=0, ax=ax,
                linewidths=0.3, cbar_kws={'label': 'Mean SUM (mean 4 folds)'})
    ax.set_title(f'{title}\n(red = ALT up, blue = ALT down)', fontsize=11)
    ax.tick_params(axis='x', labelsize=7, rotation=90)
    ax.tick_params(axis='y', labelsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_DIR, fname), dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  → {fname}  ({len(unique)} tissues found, showing top {n_show})')

# RNA tissue extractor: strip gender/age suffix
def rna_tissue(label):
    body = re.sub(r'^RNA:', '', label)
    body = re.sub(r'\s+(male|female)\s+.*$', '', body, flags=re.IGNORECASE)
    body = re.sub(r'\s+(donor|patient|individual).*$', '', body, flags=re.IGNORECASE)
    return body.strip().lower()

RNA_EXCLUDE = [r'cell line:', r'encode', r'mock', r'treated with', r'genetically modified',
               r'carcinoma', r'leukemia', r'lymphoma', r'sarcoma', r'blastoma', r'tumor',
               r'cancer cell', r'adenocarcinoma']

print('Plot 4b: RNA tissue heatmap ...')
sad_rna  = sad_mean[:, types == 'RNA']
lbls_rna = [track_labels[i] for i in range(n_tracks) if types[i] == 'RNA']
tissue_heatmap_generic(sad_rna, lbls_rna, rna_tissue,
                       'RNA-seq tissue aggregation: top 40 cell/tissue types',
                       'rna_tissue_heatmap.png', exclude_patterns=RNA_EXCLUDE)

# DNASE tissue extractor: same pattern as RNA
def dnase_tissue(label):
    body = re.sub(r'^DNASE:', '', label)
    body = re.sub(r'\s+(male|female)\s+.*$', '', body, flags=re.IGNORECASE)
    body = re.sub(r'\s+(donor|patient).*$', '', body, flags=re.IGNORECASE)
    return body.strip().lower()

DNASE_EXCLUDE = [r'cell line:', r'gm\d+', r'encode', r'mock', r'treated with',
                 r'genetically modified', r'infection', r'carcinoma', r'leukemia',
                 r'lymphoma', r'sarcoma', r'blastoma', r'tumor', r'cancer cell']

print('Plot 4c: DNASE tissue heatmap ...')
sad_dnase  = sad_mean[:, types == 'DNASE']
lbls_dnase = [track_labels[i] for i in range(n_tracks) if types[i] == 'DNASE']
tissue_heatmap_generic(sad_dnase, lbls_dnase, dnase_tissue,
                       'DNASE-seq tissue aggregation: top 40 cell/tissue types',
                       'dnase_tissue_heatmap.png', exclude_patterns=DNASE_EXCLUDE)

# ATAC: format is "ATAC:sampleID / tissue / cell_subtype" — use tissue + cell_subtype
def atac_tissue(label):
    body = re.sub(r'^ATAC:', '', label)
    parts = [p.strip() for p in body.split('/')]
    if len(parts) >= 3:
        return f"{parts[1]} – {parts[2]}"
    elif len(parts) == 2:
        return parts[1]
    return body.strip().lower()

print('Plot 4d: ATAC tissue heatmap ...')
sad_atac  = sad_mean[:, types == 'ATAC']
lbls_atac = [track_labels[i] for i in range(n_tracks) if types[i] == 'ATAC']
tissue_heatmap_generic(sad_atac, lbls_atac, atac_tissue,
                       'ATAC-seq tissue/cell-type aggregation: top 40',
                       'atac_tissue_heatmap.png')

# ── Plot 5: Top 10 tracks for the 10 most impacted SNPs ──────────────────────
print('Plot 5: Top tracks for top SNPs ...')
top10_snp_idx = np.argsort(max_abs_sad)[::-1][:10]

fig, axes = plt.subplots(2, 5, figsize=(22, 8))
axes = axes.flatten()

def format_track_label(label, maxlen=38):
    """[TYPE] description — keeps experiment type visible."""
    parts = label.split(':', 1)
    if len(parts) == 2:
        return f'[{parts[0]}] {parts[1][:maxlen]}'
    return label[:maxlen + 8]

for ax_i, snp_idx in enumerate(top10_snp_idx):
    sad_snp   = sad_mean[snp_idx, :]
    top10_tr  = np.argsort(np.abs(sad_snp))[::-1][:10]
    vals      = sad_snp[top10_tr]
    lbls      = [format_track_label(track_labels[i]) for i in top10_tr]
    colors_bar = ['#d73027' if v > 0 else '#4575b4' for v in vals]

    axes[ax_i].barh(range(10), vals[::-1], color=colors_bar[::-1])
    axes[ax_i].set_yticks(range(10))
    axes[ax_i].set_yticklabels(lbls[::-1], fontsize=7)
    axes[ax_i].axvline(0, color='black', lw=0.5)
    axes[ax_i].set_title(snp_labels[snp_idx], fontsize=9)
    axes[ax_i].set_xlabel('SUM', fontsize=8)

plt.suptitle('Top 10 tracks for the 10 most impacted SNPs\n(red = ALT up, blue = ALT down)',
             fontsize=12)
plt.tight_layout()
plt.savefig(os.path.join(PLOT_DIR, 'top_snps_tracks.png'), dpi=150)
plt.close()
print('  → top_snps_tracks.png')

print(f'\nAll outputs in: {BORZOI_RESULTS}')
print(f'Plots in:       {PLOT_DIR}')
