"""
"Mountain" / TAD-style plot per hit:
  - top   : LD matrix (r^2) as a 45-degree triangular heatmap
  - bottom: SNP association track (-log10 p, fine-mapping ADD test) in the region

LD from SuSiE_ld/{hit}/ (plink2 .phased.vcor1 = signed r, KING-cutoff cohort).
p-values / positions from fine_mapping_new/{hit}/ (ADD rows), matched by SNP ID.

Output: ld_mountain_plots/{hit}_mountain.png / .pdf
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

LD_DIR   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld'
FM_DIR   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
OUT_DIR  = os.path.join(os.path.dirname(__file__), 'ld_mountain_plots')

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})


def load_fm_pvals(hit):
    """SNP ID -> (POS, P) from the fine-mapping ADD rows."""
    rows = []
    for f in glob.glob(os.path.join(FM_DIR, hit, '*.glm.logistic.hybrid')):
        d = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
        d = d[d['TEST'] == 'ADD'][['ID', 'POS', 'P']]
        rows.append(d)
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True).drop_duplicates('ID')
    d['POS'] = pd.to_numeric(d['POS'], errors='coerce')
    d['P']   = pd.to_numeric(d['P'], errors='coerce')
    return d.set_index('ID')


def make_plot(hit):
    ld_file  = os.path.join(LD_DIR, hit, f'{hit}_ld.phased.vcor1')
    var_file = ld_file + '.vars'
    if not os.path.exists(ld_file):
        print(f'  {hit}: no LD matrix, skip'); return

    snps = pd.read_csv(var_file, header=None)[0].tolist()
    N = len(snps)
    print(f'  {hit}: {N} SNPs — loading LD matrix ...')
    R = np.loadtxt(ld_file)                     # signed r, N x N

    fm = load_fm_pvals(hit)
    pos = np.array([fm['POS'].get(s, np.nan) if fm is not None else np.nan for s in snps])
    P   = np.array([fm['P'].get(s, np.nan)   if fm is not None else np.nan for s in snps])
    mlog = -np.log10(np.where((P > 0) & np.isfinite(P), P, np.nan))

    # ── figure: heatmap (top) + p-value track (bottom), shared SNP-index x ──────
    fig, (axm, axp) = plt.subplots(
        2, 1, figsize=(9, 5.5), height_ratios=[3, 1.4],
        gridspec_kw={'hspace': 0.05}, sharex=True)

    # rotated triangular mesh: node (i,j) -> ((i+j)/2, (j-i)/2)
    ii, jj = np.meshgrid(np.arange(N + 1), np.arange(N + 1), indexing='ij')
    X = (ii + jj) / 2.0
    Y = (jj - ii) / 2.0
    C = R.copy().astype(float)                  # signed r (red = +, blue = -)
    C[np.tril_indices(N, -1)] = np.nan          # keep upper triangle only
    pc = axm.pcolormesh(X, Y, C, cmap='RdBu_r', vmin=-1, vmax=1,
                        shading='flat', rasterized=True)
    axm.set_ylim(0, N / 2.0)
    axm.set_xlim(0, N)
    axm.set_yticks([])
    for s in ('top', 'right', 'left'):
        axm.spines[s].set_visible(False)
    axm.set_title(hit, fontsize=9, fontweight='bold')
    cb = fig.colorbar(pc, ax=axm, fraction=0.03, pad=0.01)
    cb.set_label(r'LD $r$', fontsize=7); cb.ax.tick_params(labelsize=6)

    # p-value track
    axp.scatter(np.arange(N), mlog, s=6, color='#333333', linewidths=0, rasterized=True)
    axp.set_ylabel(r'$-\log_{10}(p)$', fontsize=8)
    axp.set_ylim(bottom=0)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)

    # x ticks: SNP index -> genomic position (Mb)
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
