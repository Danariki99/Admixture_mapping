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
VCF_FILE   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz'
MSP_DIR    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb'
OUT_DIR    = os.path.dirname(os.path.abspath(__file__))
SNP_CACHE  = os.path.join(OUT_DIR, 'genomic_coverage_snp_positions.tsv')

# ── 1. SNP positions: use cache if present, else extract from VCF (slow, ~30 min) ─
if os.path.exists(SNP_CACHE):
    print(f'Using cached SNP positions: {SNP_CACHE}')
    snps = pd.read_csv(SNP_CACHE, sep='\t', dtype={'chrom': int, 'pos': int})
else:
    print('Extracting SNP positions from VCF (bcftools query, ~30 min on 66GB VCF)...')
    result = subprocess.run(
        ['bcftools', 'query', '-f', '%CHROM\t%POS\n', VCF_FILE],
        capture_output=True, text=True, check=True
    )
    snps = pd.read_csv(
        pd.io.common.StringIO(result.stdout),
        sep='\t', header=None, names=['chrom', 'pos'],
        dtype={'chrom': int, 'pos': int}
    )
    snps.to_csv(SNP_CACHE, sep='\t', index=False)
    print(f'  Cached SNP positions → {SNP_CACHE}')
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

# ── 4. Distribution of SNPs per window ────────────────────────────────────────
print('Plotting distribution...')
MM_TO_INCH = 1 / 25.4

import matplotlib as mpl
from scipy.stats import gaussian_kde
mpl.rcParams.update({
    'font.size':          6,
    'axes.titlesize':     6,
    'axes.labelsize':     6,
    'xtick.labelsize':    5,
    'ytick.labelsize':    5,
    'pdf.fonttype':       42,
    'savefig.dpi':        600,
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
})

vals = wins['n_snps'].values.astype(float)
med  = np.median(vals)
print('  SNPs/window  median={:.0f}  IQR={:.0f}-{:.0f}  p99={:.0f}  max={:.0f}'.format(
    med, np.percentile(vals, 25), np.percentile(vals, 75),
    np.percentile(vals, 99), vals.max()))

# robust x-axis: focus on the bulk, the rare high-count windows stay in a thin tail
x_top = np.percentile(vals, 99) * 1.05
NBINS = 75         # middle ground between blocky (60) and too fine (100)

fig, ax = plt.subplots(figsize=(90 * MM_TO_INCH, 60 * MM_TO_INCH))

# histogram: light bars as a soft backdrop
counts, edges, _ = ax.hist(vals, bins=NBINS, range=(0, x_top),
                           color='#82aed0', edgecolor='white', linewidth=0.2, alpha=0.7)
bin_width = edges[1] - edges[0]

# smooth KDE overlay (mild extra smoothing so it still tracks the peak)
kde = gaussian_kde(vals)
kde.set_bandwidth(kde.factor * 1.15)
xs = np.linspace(0, x_top, 400)
ax.plot(xs, kde(xs) * len(vals) * bin_width, color='#3d6aad', lw=1.2)

# median reference
ax.axvline(med, color='#c0392b', ls='--', lw=0.7,
           label=f'median = {med:.0f} SNPs/window')

ax.set_xlim(0, x_top)
ax.set_title('UKBB SNP coverage per window', fontsize=6)
ax.set_xlabel('SNPs per window', fontsize=6)
ax.set_ylabel('Number of windows', fontsize=6)
ax.tick_params(axis='both', labelsize=5)
ax.legend(loc='upper right', fontsize=5, frameon=False)
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
