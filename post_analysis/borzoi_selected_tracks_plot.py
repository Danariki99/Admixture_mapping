"""
Curated Borzoi SED figure: for each lead/causal SNP, the per-gene SNP Expression
Difference (SED) across a hand-picked set of relevant RNA tracks.

Values are taken from sed_mean_long.csv (RNA tracks of the relevant cell types;
GM12878 cell line, whole spleen/thymus and activated/stimulated samples excluded).
Bars coloured by direction: orange = ALT increases expression (+), light blue =
ALT decreases expression (-). One panel per SNP, 1x4 row at A4 width, no title.

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

POS_COLOR = '#f4a582'   # ALT up   (light orange)
NEG_COLOR = '#9ecae1'   # ALT down (light blue)

# Curated tracks per SNP: (compact label, SAD value)
DATA = {
    # rs3177928 — HLA-DRA, consistent decrease in lung
    'rs3177928': {
        'title': 'rs3177928\nHLA-DRA',
        'tracks': [
            ('lung, F 30y',        -71.11),
            ('lung, M 3y',         -63.39),
            ('lung, F 47y',        -35.39),
            ('lung, F 30y (rep)',  -34.18),
            ('lung, embryo',       -10.32),
            ('lung',                -7.35),
        ],
    },
    # rs28732226 — HLA-DRA (discordant, large +) and HLA-DRB9 (concordant, small -)
    'rs28732226': {
        'title': 'rs28732226\nHLA-DRA / HLA-DRB9',
        'tracks': [
            ('DRA · CD14 monocyte',     23.31),
            ('DRA · IgD mem B cell',    20.60),
            ('DRA · naive B cell',      19.79),
            ('DRA · immature NK cell',  16.18),
            ('DRA · CD8 mem T cell',    13.90),
            ('DRA · lung, M 3y',         5.28),
            ('DRA · lung, F 30y',        4.61),
            ('DRB9 · naive B cell',     -0.13),
            ('DRB9 · IgD mem B cell',   -0.08),
            ('DRB9 · CD14 monocyte',    -0.05),
        ],
    },
    # rs9274569 — HLA-DQB1, lymphoid/myeloid (+), thyroid (-)
    'rs9274569': {
        'title': 'rs9274569\nHLA-DQB1',
        'tracks': [
            ('IgD mem B cell',      22.30),
            ('naive B cell',        21.14),
            ('CD8 mem T cell',      19.49),
            ('CD4 mem T cell',      19.11),
            ('naive CD4 T cell',    17.66),
            ('CD14 monocyte',       12.38),
            ('thyroid, F 53y',      -1.39),
            ('thyroid, M 54y',      -1.21),
            ('thyroid',             -0.86),
        ],
    },
    # rs6906021 — HLA-DQB1, lymphoid/myeloid (-), lung (-), thyroid (-)
    'rs6906021': {
        'title': 'rs6906021\nHLA-DQB1',
        'tracks': [
            ('IgD mem B cell',     -100.06),
            ('naive B cell',        -88.19),
            ('immature NK cell',    -84.08),
            ('CD8 mem T cell',      -75.82),
            ('CD14 monocyte',       -65.32),
            ('lung, F 30y',         -14.20),
            ('lung, M 3y',          -11.63),
            ('thyroid, F 51y',       -0.90),
            ('thyroid, M 54y',       -0.78),
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
    fmt = lambda v: f'{v:+.2f}' if abs(v) < 1 else f'{v:+.1f}'
    for yi, v in zip(y, vals):
        if v >= 0:
            ax.text(v + span * 0.03, yi, fmt(v), va='center', ha='left', fontsize=4.5)
        else:
            ax.text(v - span * 0.03, yi, fmt(v), va='center', ha='right', fontsize=4.5)

    ax.set_xlim(min(0, min(vals)) - span * 0.30, max(0, max(vals)) + span * 0.30)
    ax.set_title(info['title'], fontsize=6, fontweight='bold')
    ax.set_xlabel('SED (ALT − REF)', fontsize=5)
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

    # no figure title (matches the Inkscape template)
    fig.tight_layout()

    base = os.path.join(OUT_DIR, 'selected_tracks')
    fig.savefig(f'{base}.png', bbox_inches='tight')
    fig.savefig(f'{base}.pdf', bbox_inches='tight', format='pdf')
    plt.close(fig)
    print(f'Saved → {base}.png / .pdf')


if __name__ == '__main__':
    main()
