"""
Curated Borzoi figure: for each lead/causal SNP, show a hand-picked set of
tracks (not the automatic top-N), with the SAD values.

Bars coloured by direction: red = ALT increases signal (+), blue = ALT decreases
signal (-). One panel per SNP, 1x4 row at A4 width. Values are curated SAD scores.

Output: borzoi_results/plots_causal_candidates/selected_tracks.{png,pdf}
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

OUT_DIR    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/borzoi_results/plots_causal_candidates'
MM_TO_INCH = 1 / 25.4
A4_WIDTH_MM = 210

mpl.rcParams.update({
    'font.size':       6,
    'axes.titlesize':  6,
    'axes.labelsize':  5,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'pdf.fonttype':    42,
    'ps.fonttype':     42,
    'savefig.dpi':     600,
    'font.family':     'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
})

POS_COLOR = '#d73027'   # ALT up
NEG_COLOR = '#4575b4'   # ALT down

# Curated tracks per SNP: (compact label, SAD value)
DATA = {
    'rs3177928': {
        'title': 'rs3177928\nlung',
        'tracks': [
            ('lung, F 30y',          -83.30),
            ('lung, M 3y',           -82.96),
            ('L-lung upper, M 60y',  -51.99),
            ('L-lung upper, F 61y',  -50.21),
            ('L-lung lower, M 60y',  -48.93),
            ('R-lung lower, M 60y',  -46.04),
        ],
    },
    'rs28732226': {
        'title': 'rs28732226\nH3K4me3 / monocytes / B cells',
        'tracks': [
            ('H3K4me3 GM12878 (ChIP)', 28.05),
            ('CD14 monocyte (rep)',    21.49),
            ('IgD mem B cell',         17.39),
            ('CD14 monocyte (MS)',     16.64),
            ('naive B cell',           16.95),
            ('CD14 monocyte',          15.51),
        ],
    },
    'rs9274569': {
        'title': 'rs9274569\nlymphoid/myeloid (+); thyroid (−)',
        'tracks': [
            ('IgD mem B cell',      66.03),
            ('naive B cell',        64.88),
            ('CD14 monocyte',       59.47),
            ('immature NK cell',    53.76),
            ('naive CD4 T cell',    37.14),
            ('thyroid',             -0.85),
            ('thyroid (rep)',       -0.80),
        ],
    },
    'rs6906021': {
        'title': 'rs6906021\nlymphoid/myeloid (−); lung; thyroid',
        'tracks': [
            ('IgD mem B cell',      -173.30),
            ('naive B cell',        -153.88),
            ('immature NK cell',    -143.74),
            ('CD14 monocyte',       -118.57),
            ('CD8 mem T (MS)',       -77.89),
            ('lung, F 30y',         -27.09),
            ('L-lung lower, M 60y',  -9.93),
            ('thyroid, F embryo',    -2.13),
            ('thyroid, M 54y',       -1.22),
        ],
    },
}


def plot_snp(ax, info):
    # sort ascending by value → most positive at top, most negative at bottom
    tracks = sorted(info['tracks'], key=lambda t: t[1])
    labels = [t[0] for t in tracks]
    vals   = [t[1] for t in tracks]
    colors = [POS_COLOR if v >= 0 else NEG_COLOR for v in vals]

    y = np.arange(len(vals))
    ax.barh(y, vals, color=colors, height=0.72, linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=5)
    ax.axvline(0, color='black', lw=0.5)

    span = max(abs(min(vals)), abs(max(vals)))
    for yi, v in zip(y, vals):
        if v >= 0:
            ax.text(v + span * 0.03, yi, f'{v:+.1f}', va='center', ha='left', fontsize=4.5)
        else:
            ax.text(v - span * 0.03, yi, f'{v:+.1f}', va='center', ha='right', fontsize=4.5)

    ax.set_xlim(min(0, min(vals)) - span * 0.30, max(0, max(vals)) + span * 0.30)
    ax.set_title(info['title'], fontsize=6, fontweight='bold')
    ax.set_xlabel('SAD (ALT − REF)', fontsize=5)
    ax.tick_params(axis='both', width=0.5, length=2)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['bottom'].set_linewidth(0.5)


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    snps = list(DATA.keys())

    fig, axes = plt.subplots(
        1, 4, figsize=(A4_WIDTH_MM * MM_TO_INCH, 80 * MM_TO_INCH),
        gridspec_kw={'wspace': 0.9},
    )
    for ax, snp in zip(axes, snps):
        plot_snp(ax, DATA[snp])

    fig.suptitle('Borzoi predicted variant effects on selected tracks '
                 '(red = ALT ↑, blue = ALT ↓)', fontsize=7, y=1.06)
    fig.tight_layout(rect=[0, 0, 1, 0.9])

    base = os.path.join(OUT_DIR, 'selected_tracks')
    fig.savefig(f'{base}.png', bbox_inches='tight')
    fig.savefig(f'{base}.pdf', bbox_inches='tight', format='pdf')
    plt.close(fig)
    print(f'Saved → {base}.png / .pdf')


if __name__ == '__main__':
    main()
