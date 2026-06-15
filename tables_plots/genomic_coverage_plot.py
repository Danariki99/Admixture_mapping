"""
Genome-wide SNP density plot: number of UKBB genotyped SNPs per RFMix LAI window.

For each LAI window defined in the MSP files (all chromosomes), counts how many
UKBB genotyped SNPs fall within it. Alternating colours per chromosome.

Output: genomic_coverage_plot.png and .pdf (same folder as this script)
"""

import os
import glob
import subprocess
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── Paths ─────────────────────────────────────────────────────────────────────
VCF_FILE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz'
MSP_DIR  = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb'
OUT_DIR  = os.path.dirname(os.path.abspath(__file__))

# ── 1. Extract SNP positions from VCF ─────────────────────────────────────────
print('Extracting SNP positions from VCF (bcftools query)...')
result = subprocess.run(
    ['bcftools', 'query', '-f', '%CHROM\t%POS\n', VCF_FILE],
    capture_output=True, text=True, check=True
)
snps = pd.read_csv(
    pd.io.common.StringIO(result.stdout),
    sep='\t', header=None, names=['chrom', 'pos'],
    dtype={'chrom': int, 'pos': int}
)
print(f'  Total SNPs: {len(snps):,}')

snps_by_chrom = {c: np.sort(grp['pos'].values) for c, grp in snps.groupby('chrom')}

# ── 2. Read MSP windows (first 6 columns only) ────────────────────────────────
print('Reading MSP windows...')
windows = []
for path in sorted(glob.glob(os.path.join(MSP_DIR, '*.msp.tsv'))):
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.split('\t', 7)
            windows.append((int(parts[0]), int(parts[1]), int(parts[2])))

wins = pd.DataFrame(windows, columns=['chrom', 'spos', 'epos'])
wins = wins.sort_values(['chrom', 'spos']).reset_index(drop=True)
print(f'  Total MSP windows: {len(wins):,}')

# ── 3. Count VCF SNPs per window ──────────────────────────────────────────────
print('Counting SNPs per window...')
def count_snps(chrom, spos, epos):
    arr = snps_by_chrom.get(chrom)
    if arr is None:
        return 0
    return int(np.searchsorted(arr, epos, side='right') -
               np.searchsorted(arr, spos, side='left'))

wins['n_snps'] = wins.apply(lambda r: count_snps(r['chrom'], r['spos'], r['epos']), axis=1)

# ── 4. Genome-wide x-axis ─────────────────────────────────────────────────────
chrom_sizes = wins.groupby('chrom')['epos'].max().to_dict()
chroms = sorted(chrom_sizes.keys())

offsets, cum = {}, 0
for c in chroms:
    offsets[c] = cum
    cum += chrom_sizes[c]

wins['x_s'] = wins.apply(lambda r: offsets[r['chrom']] + r['spos'], axis=1)
wins['x_e'] = wins.apply(lambda r: offsets[r['chrom']] + r['epos'], axis=1)
chrom_centers = {c: offsets[c] + chrom_sizes[c] / 2 for c in chroms}

# ── 5. Alternating colours per chromosome ─────────────────────────────────────
palette = ['#3d6aad', '#82aed0']
colors = [palette[chroms.index(r['chrom']) % 2] for _, r in wins.iterrows()]

# ── 6. Plot ────────────────────────────────────────────────────────────────────
print('Plotting...')
MM_TO_INCH = 1 / 25.4

import matplotlib as mpl
mpl.rcParams.update({
    'font.size':          6,
    'axes.titlesize':     6,
    'axes.labelsize':     6,
    'xtick.labelsize':    5,
    'ytick.labelsize':    5,
    'pdf.fonttype':       42,
    'savefig.dpi':        600,
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
})

fig, ax = plt.subplots(figsize=(183 * MM_TO_INCH, 55 * MM_TO_INCH))

for i, (_, row) in enumerate(wins.iterrows()):
    ax.bar(row['x_s'], row['n_snps'],
           width=max(row['x_e'] - row['x_s'], 1),
           align='edge', color=colors[i], linewidth=0)

for c in chroms[1:]:
    ax.axvline(offsets[c], color='white', lw=0.6, zorder=2)

ax.set_xticks([chrom_centers[c] for c in chroms])
ax.set_xticklabels([str(c) for c in chroms])
ax.set_xlim(0, cum)
ax.set_title('UKBB genome-wide SNP coverage', fontsize=6)
ax.set_xlabel('Chromosome', fontsize=6)
ax.set_ylabel('SNPs per window', fontsize=6)
ax.tick_params(axis='both', labelsize=5)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()

out_png = os.path.join(OUT_DIR, 'genomic_coverage_plot.png')
out_pdf = os.path.join(OUT_DIR, 'genomic_coverage_plot.pdf')
plt.savefig(out_png, bbox_inches='tight', pad_inches=0.01)
plt.savefig(out_pdf, bbox_inches='tight', pad_inches=0.01, format='pdf')
print(f'  Saved → {out_png}')
print(f'  Saved → {out_pdf}')

plt.close()
print('Done.')
