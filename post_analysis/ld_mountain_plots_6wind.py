"""
"Mountain" / TAD-style plot per hit — EXTENDED (6-window) version:
  - top   : LD matrix (r^2) as a 45-degree triangular heatmap
  - bottom: SNP association track (-log10 p, fine-mapping ADD) in the region

Uses the extended ±6-window LD (SuSiE_ld_6wind/{hit}/{hit}_ld_6wind.phased.vcor1).
p-values / positions come from BOTH fine_mapping_new (significant windows) and
fine_mapping_new_6wind (extension windows), matched by SNP ID.

Output: ld_mountain_plots_6wind/{hit}_mountain.png / .pdf
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

LD_DIR   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld_6wind'
FM_DIRS  = ['/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new',
            '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind']
OUT_DIR  = os.path.join(os.path.dirname(__file__), 'ld_mountain_plots_6wind')

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})


def load_fm_pvals(hit):
    """SNP ID -> (POS, OR) from ADD rows across the significant + extension windows."""
    rows = []
    for fm_dir in FM_DIRS:
        for f in glob.glob(os.path.join(fm_dir, hit, '*.glm.logistic.hybrid')):
            d = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
            d = d[d['TEST'] == 'ADD'][['ID', 'POS', 'OR']]
            rows.append(d)
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True).drop_duplicates('ID')
    d['POS'] = pd.to_numeric(d['POS'], errors='coerce')
    d['OR']  = pd.to_numeric(d['OR'], errors='coerce')
    return d.set_index('ID')


def make_plot(hit):
    ld_file  = os.path.join(LD_DIR, hit, f'{hit}_ld_6wind.phased.vcor1')
    var_file = ld_file + '.vars'
    if not os.path.exists(ld_file):
        print(f'  {hit}: no 6wind LD matrix, skip'); return

    snps = pd.read_csv(var_file, header=None)[0].tolist()
    N = len(snps)
    print(f'  {hit}: {N} SNPs — loading LD matrix ...')
    R = np.loadtxt(ld_file)                     # signed r, N x N
    ld2 = R ** 2                                 # r^2

    fm = load_fm_pvals(hit)
    pos = np.array([fm['POS'].get(s, np.nan) if fm is not None else np.nan for s in snps])
    ORv = np.array([fm['OR'].get(s, np.nan)  if fm is not None else np.nan for s in snps])

    fig, (axm, axp) = plt.subplots(
        2, 1, figsize=(9, 5.5), height_ratios=[3, 1.4],
        gridspec_kw={'hspace': 0.05}, sharex=True)

    ii, jj = np.meshgrid(np.arange(N + 1), np.arange(N + 1), indexing='ij')
    X = (ii + jj) / 2.0
    Y = (jj - ii) / 2.0
    C = ld2.copy().astype(float)
    C[np.tril_indices(N, -1)] = np.nan          # keep upper triangle only
    pc = axm.pcolormesh(X, Y, C, cmap='Reds', vmin=0, vmax=1,
                        shading='flat', rasterized=True)
    axm.set_ylim(0, N / 2.0)
    axm.set_xlim(0, N)
    axm.set_yticks([])
    for s in ('top', 'right', 'left'):
        axm.spines[s].set_visible(False)
    axm.set_title(hit + '  (±6 windows)', fontsize=9, fontweight='bold')
    # colorbar in an external inset so it does NOT shrink the triangle
    # (keeps top and bottom panels the same width -> aligned)
    cax = inset_axes(axm, width='2%', height='65%', loc='lower left',
                     bbox_to_anchor=(1.015, 0.0, 1, 1), bbox_transform=axm.transAxes,
                     borderpad=0)
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r'LD $r^2$', fontsize=7); cb.ax.tick_params(labelsize=6)

    # OR bar plot: bars emanate from the OR = 1 null line (red = risk, blue = protective)
    x = np.arange(N)
    valid = np.isfinite(ORv) & (ORv > 0)
    colors = np.where(ORv[valid] >= 1, '#d62728', '#1f77b4')
    axp.bar(x[valid], (ORv[valid] - 1), bottom=1.0, width=1.0, color=colors,
            linewidth=0, rasterized=True)
    axp.axhline(1.0, color='black', lw=0.5)
    axp.set_ylabel('OR', fontsize=8)
    # robust y-range so rare-variant outliers don't flatten the rest
    lo, hi = np.nanpercentile(ORv[valid], [1, 99])
    axp.set_ylim(min(0.9, lo), max(1.1, hi))
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)

    nt = 6
    idx = np.linspace(0, N - 1, nt).astype(int)
    labs = [f'{pos[i]/1e6:.2f}' if np.isfinite(pos[i]) else '' for i in idx]
    axp.set_xticks(idx); axp.set_xticklabels(labs, fontsize=6)
    chrom = hit.split('_')[-1]
    axp.set_xlabel(f'Position on {chrom} (Mb)', fontsize=8)
    axp.tick_params(axis='y', labelsize=6)

    os.makedirs(OUT_DIR, exist_ok=True)
    base = os.path.join(OUT_DIR, f'{hit}_mountain')
    fig.savefig(f'{base}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{base}.pdf', dpi=300, bbox_inches='tight', format='pdf')
    plt.close(fig)
    print(f'    saved -> {base}.png / .pdf')


if __name__ == '__main__':
    import sys
    hits = sys.argv[1:] or sorted(os.listdir(LD_DIR))
    for hit in hits:
        if os.path.isdir(os.path.join(LD_DIR, hit)):
            make_plot(hit)
