import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

BIN_SIZE = 1_000_000  # 1 Mb bins

OUTPUT_FOLDER = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots'
INPUT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/output'

chromosome_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78',
    '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7',
    '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
]

# Collect unique (CHROM, POS) across all GWAS output files
snp_set = set()
for fold_name in sorted(os.listdir(INPUT_FOLDER)):
    fold_path = os.path.join(INPUT_FOLDER, fold_name)
    if not os.path.isdir(fold_path):
        continue
    for fname in os.listdir(fold_path):
        if not fname.endswith('.glm.logistic.hybrid') or fname.endswith('.adjusted'):
            continue
        fpath = os.path.join(fold_path, fname)
        df = pd.read_csv(fpath, sep='\t', usecols=['#CHROM', 'POS'])
        df = df[pd.to_numeric(df['#CHROM'], errors='coerce').notna()]
        df['#CHROM'] = df['#CHROM'].astype(int)
        # keep autosomes only (1-22)
        df = df[df['#CHROM'].between(1, 22)]
        snp_set.update(zip(df['#CHROM'], df['POS'].astype(int)))
    print(f'  loaded {fold_name}')

print(f'Total unique SNPs: {len(snp_set):,}')

snp_df = pd.DataFrame(list(snp_set), columns=['CHROM', 'POS'])
snp_df = snp_df.sort_values(['CHROM', 'POS'])

snp_df['BIN'] = (snp_df['POS'] // BIN_SIZE) * BIN_SIZE

bin_counts = snp_df.groupby(['CHROM', 'BIN']).size().reset_index(name='COUNT')

# Chromosome offsets: lay chromosomes end-to-end with a small gap
chrom_max_pos = snp_df.groupby('CHROM')['POS'].max()
chromosomes = list(range(1, 23))
chrom_offsets = {}
offset = 0
for chrom in chromosomes:
    chrom_offsets[chrom] = offset
    offset += int(chrom_max_pos.get(chrom, 250_000_000)) + 10_000_000  # 10 Mb gap

bin_counts['ABS_POS'] = bin_counts.apply(
    lambda r: r['BIN'] + chrom_offsets.get(r['CHROM'], 0), axis=1
)

fig, ax = plt.subplots(figsize=(14, 5))

for i, chrom in enumerate(chromosomes):
    chrom_data = bin_counts[bin_counts['CHROM'] == chrom]
    if chrom_data.empty:
        continue
    color = chromosome_colors[i % len(chromosome_colors)]
    ax.bar(
        chrom_data['ABS_POS'],
        chrom_data['COUNT'],
        width=BIN_SIZE * 0.85,
        color=color,
        linewidth=0,
    )

# Chromosome tick labels at midpoint of each chromosome
chrom_ticks = []
chrom_labels = []
for chrom in chromosomes:
    chrom_data = bin_counts[bin_counts['CHROM'] == chrom]
    if chrom_data.empty:
        continue
    mid = (chrom_data['ABS_POS'].min() + chrom_data['ABS_POS'].max()) / 2
    chrom_ticks.append(mid)
    chrom_labels.append(str(chrom))

ax.set_xticks(chrom_ticks)
ax.set_xticklabels(chrom_labels, fontsize=8)
ax.set_xlabel('Chromosome', fontsize=11)
ax.set_ylabel(f'SNPs per {BIN_SIZE // 1_000_000} Mb bin', fontsize=11)
ax.set_title('Genomic SNP Coverage — UKBB', fontsize=13)
ax.set_xlim(0, offset)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

out_png = os.path.join(OUTPUT_FOLDER, 'snp_coverage_ukbb.png')
out_pdf = os.path.join(OUTPUT_FOLDER, 'snp_coverage_ukbb.pdf')
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)
plt.close()

print(f'Saved: {out_png}')
print(f'Saved: {out_pdf}')
