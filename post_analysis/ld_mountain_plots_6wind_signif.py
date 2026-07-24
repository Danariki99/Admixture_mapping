"""
"Mountain" / TAD-style plot per hit — EXTENDED (6-window), FINE-MAP-coloured:
  - top   : LD matrix (r^2) as a 45-degree triangular heatmap
  - bottom: BY-corrected association track as a bar plot. EVERY SNP is drawn
            (height = -log10 of the BY-corrected p), all bars pointing UP.
            Only the SNPs that are BOTH (i) significant in the fine mapping
            (i.e. fine-mapping candidates, fine_mapping_all_candidates.tsv) AND
            (ii) CONCORDANT with the original admixture signal direction are
            coloured; everything else stays grey. Because concordant SNPs share
            the hit's direction, each plot uses at most TWO colours:
                grey  = all other SNPs
                red   = fine-mapped & concordant, risk hit       (OR >= 1)
                blue  = fine-mapped & concordant, protective hit (OR <  1)

Companion to ld_mountain_plots_6wind.py (which plots the OR itself as signed
bars). Uses the same extended +/-6-window LD and the same fine-mapping p-values.

Output: ld_mountain_plots_6wind_signif/{hit}_mountain.png / .pdf
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from statsmodels.stats.multitest import multipletests

LD_DIR   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_ld_6wind'
FM_DIRS  = ['/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new',
            '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new_6wind']
OUT_DIR  = os.path.join(os.path.dirname(__file__), 'ld_mountain_plots_6wind_signif')

# Fine-mapping candidates (already carry OR, admixture_OR and concordant_direction)
CAND_FILE = ('/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/'
             'ukbb/fine_mapping_conditional_results/fine_mapping_all_candidates.tsv')

BY_ALPHA = 0.05     # BY threshold, drawn only as a reference guide line

mpl.rcParams.update({'pdf.fonttype': 42, 'font.family': 'sans-serif',
                     'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans']})

# hit -> {SNP ID: OR} for the fine-mapped candidates that are concordant with
# the admixture signal (the only SNPs we colour).
_cand = pd.read_csv(CAND_FILE, sep='\t')
_cand = _cand[_cand['concordant_direction'] == True]                     # noqa: E712
CONCORDANT = {h: dict(zip(g['ID'], g['OR'])) for h, g in _cand.groupby('hit')}


def load_fm_pvals(hit):
    """SNP ID -> (POS, P, OR) from ADD rows across significant + extension windows."""
    rows = []
    for fm_dir in FM_DIRS:
        for f in glob.glob(os.path.join(fm_dir, hit, '*.glm.logistic.hybrid')):
            d = pd.read_csv(f, sep='\t', dtype=str, low_memory=False)
            d = d[d['TEST'] == 'ADD'][['ID', 'POS', 'P', 'OR']]
            rows.append(d)
    if not rows:
        return None
    d = pd.concat(rows, ignore_index=True).drop_duplicates('ID')
    for c in ('POS', 'P', 'OR'):
        d[c] = pd.to_numeric(d[c], errors='coerce')
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
    P   = np.array([fm['P'].get(s, np.nan)   if fm is not None else np.nan for s in snps])
    ORv = np.array([fm['OR'].get(s, np.nan)  if fm is not None else np.nan for s in snps])

    # ── BY correction over the SNPs that carry a valid p ────────────────────────
    finite = np.isfinite(P) & (P > 0)
    by_p = np.full(N, np.nan)
    if finite.sum() > 0:
        _, p_adj, _, _ = multipletests(P[finite], alpha=BY_ALPHA, method='fdr_by')
        by_p[finite] = p_adj
    height = -np.log10(np.where((by_p > 0) & np.isfinite(by_p), by_p, np.nan))

    # ── figure: heatmap (top) + BY bar track (bottom), shared SNP-index x ───────
    # gridspec with a dedicated thin colorbar column so the triangle and the bar
    # track keep exactly the same width (aligned) and it survives bbox='tight'.
    fig = plt.figure(figsize=(9, 5.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[3, 1.4], width_ratios=[45, 1],
                          hspace=0.05, wspace=0.02)
    axm = fig.add_subplot(gs[0, 0])
    axp = fig.add_subplot(gs[1, 0], sharex=axm)
    cax = fig.add_subplot(gs[0, 1])

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
    axm.tick_params(axis='x', labelbottom=False, bottom=False)   # no shared x labels on top panel
    for s in ('top', 'right', 'left'):
        axm.spines[s].set_visible(False)
    axm.set_title(hit + '  (±6 windows)', fontsize=9, fontweight='bold')
    cb = fig.colorbar(pc, cax=cax)
    cb.set_label(r'LD $r^2$', fontsize=7); cb.ax.tick_params(labelsize=6)

    # ── bar plot: every SNP drawn, only fine-mapped & concordant ones coloured ──
    x = np.arange(N)
    C_OTHER, C_RISK, C_PROT = '#bdbdbd', '#d62728', '#1f77b4'
    cand_or = CONCORDANT.get(hit, {})                       # {ID: OR} concordant candidates
    is_cand = np.array([s in cand_or for s in snps])
    colors = np.full(N, C_OTHER, dtype=object)
    colors[is_cand & (ORv >= 1)] = C_RISK                   # risk hit  -> red
    colors[is_cand & (ORv <  1)] = C_PROT                   # protective hit -> blue

    valid = np.isfinite(height)
    # draw grey first, then the coloured candidates on top so they are never hidden
    other = valid & ~is_cand
    axp.bar(x[other], height[other], width=1.0, color=C_OTHER,
            linewidth=0, rasterized=True)
    hi = valid & is_cand
    axp.bar(x[hi], height[hi], width=1.0, color=list(colors[hi]),
            linewidth=0, rasterized=True, zorder=3)
    axp.set_ylabel(r'$-\log_{10}(\mathrm{BY}\ p)$', fontsize=8)
    axp.set_ylim(bottom=0)
    for s in ('top', 'right'):
        axp.spines[s].set_visible(False)

    # legend: grey + whichever single highlight colour this hit actually uses
    handles = [Patch(facecolor=C_OTHER, label='Other SNPs')]
    if (hi & (ORv >= 1)).any():
        handles.append(Patch(facecolor=C_RISK, label='Fine-mapped · concordant (risk)'))
    if (hi & (ORv < 1)).any():
        handles.append(Patch(facecolor=C_PROT, label='Fine-mapped · concordant (protective)'))
    axp.legend(handles=handles, fontsize=6, loc='upper right', frameon=False,
               handlelength=1.0, borderaxespad=0.2)

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
