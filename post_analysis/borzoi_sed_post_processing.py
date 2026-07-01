"""
Borzoi SED post-processing: reads the 4-fold sed.h5 (gene-level SNP Expression
Difference), averages across folds, and exports per-(SNP, gene) tables.

For each SNP, every gene whose exons fall in the prediction window is scored:
a positive SED = ALT increases predicted expression, negative = decreases.

Outputs (in borzoi_results/sed_leads/):
  sed_mean_long.csv     — one row per (SNP, gene, track) with mean SED
  sed_gene_summary.csv  — per (SNP, gene): max |SED| track + which tissue
"""

import os
import re
import h5py
import numpy as np
import pandas as pd

SED_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/sed_leads'
FOLDS = ['f3c0', 'f3c1', 'f3c2', 'f3c3']

# ── Load + average folds ──────────────────────────────────────────────────────
sed_arrays = []
genes = snps = chrom = pos = ref = alt = track_labels = None

for fold in FOLDS:
    h5_path = os.path.join(SED_BASE, fold, 'sed.h5')
    print(f'Reading {h5_path} ...')
    with h5py.File(h5_path, 'r') as f:
        sed = f['SED'][()].astype(np.float32)          # (n_pairs, n_tracks)
        si  = f['si'][()]                              # snp index per (snp,gene) pair
        g   = [x.decode() for x in f['gene'][()]]
        s   = [x.decode() for x in f['snp'][()]]
        ch  = [x.decode() for x in f['chr'][()]]
        po  = f['pos'][()].tolist()
        rf  = [x.decode() for x in f['ref_allele'][()]]
        al  = [x.decode() for x in f['alt_allele'][()]]
        tl  = [x.decode() for x in f['target_labels'][()]]

    if genes is None:
        genes, snp_idx = g, si
        snps = [s[i] for i in si]
        chrom = [ch[i] for i in si]
        pos   = [po[i] for i in si]
        ref   = [rf[i] for i in si]
        alt   = [al[i] for i in si]
        track_labels = tl
    sed_arrays.append(sed)

sed_mean = np.mean(np.stack(sed_arrays, axis=0), axis=0)   # (n_pairs, n_tracks)
n_pairs, n_tracks = sed_mean.shape
print(f'\n{n_pairs} (SNP, gene) pairs x {n_tracks} tracks')

# ── ENSG → gene symbol map (from the same GENCODE GTF Borzoi used) ─────────────
GTF_PATH = '/private/home/rsmerigl/codes/cleaned_codes/borzoi/examples/hg38/genes/gencode41/gencode41_basic_nort.gtf'
import re
def build_gene_name_map(gtf_path, wanted_base_ids):
    name_map = {}
    with open(gtf_path) as fh:
        for line in fh:
            if '\tgene\t' not in line:
                continue
            gid = re.search(r'gene_id "([^"]+)"', line)
            gnm = re.search(r'gene_name "([^"]+)"', line)
            if gid and gnm:
                base = gid.group(1).split('.')[0]
                if base in wanted_base_ids:
                    name_map[base] = gnm.group(1)
    return name_map

wanted = {g.split('.')[0] for g in genes}
gene_name_map = build_gene_name_map(GTF_PATH, wanted)
gene_names = [gene_name_map.get(g.split('.')[0], g.split('.')[0]) for g in genes]
print(f'Mapped {len(gene_name_map)}/{len(wanted)} genes to symbols')

meta = pd.DataFrame({'snp': snps, 'gene': genes, 'gene_name': gene_names,
                     'chr': chrom, 'pos': pos, 'ref': ref, 'alt': alt})

# ── Long table: (SNP, gene, track) ────────────────────────────────────────────
long_rows = []
for j, tl in enumerate(track_labels):
    df = meta.copy()
    df['track'] = tl
    df['SED'] = sed_mean[:, j]
    long_rows.append(df)
long_df = pd.concat(long_rows, ignore_index=True)
long_out = os.path.join(SED_BASE, 'sed_mean_long.csv')
long_df.to_csv(long_out, index=False)
print(f'  → {long_out}')

# ── Per-(SNP, gene) summary: strongest-effect track ───────────────────────────
summ = []
for i in range(n_pairs):
    row = sed_mean[i, :]
    j = int(np.argmax(np.abs(row)))
    summ.append({
        'snp':       snps[i],
        'gene':      genes[i],
        'gene_name': gene_names[i],
        'chr':       chrom[i],
        'pos':       pos[i],
        'max_abs_SED':       row[j],
        'top_track':         track_labels[j],
        'mean_SED_alltracks': float(np.mean(row)),
    })
summ_df = pd.DataFrame(summ).sort_values(['snp', 'max_abs_SED'],
                                         key=lambda c: c.abs() if c.name == 'max_abs_SED' else c,
                                         ascending=[True, False])
summ_out = os.path.join(SED_BASE, 'sed_gene_summary.csv')
summ_df.to_csv(summ_out, index=False)
print(f'  → {summ_out}\n')

# quick console view: genes per SNP, ranked by |SED|
for snp in dict.fromkeys(snps):
    sub = summ_df[summ_df['snp'] == snp].reindex(
        summ_df[summ_df['snp'] == snp]['max_abs_SED'].abs().sort_values(ascending=False).index)
    print(f'{snp}  ({len(sub)} genes in region):')
    for _, r in sub.iterrows():
        print(f'    {r["gene_name"]:<14} {r["gene"]:<20} maxSED={r["max_abs_SED"]:+.3f}  [{r["top_track"]}]')
    print()
