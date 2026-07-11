"""
Multi-ancestry Manhattan panel where only the FP1-significant windows are
coloured — one colour per phenotype, shared GLOBALLY across ancestries so that
the same phenotype always gets the same colour.

For each ancestry that has at least one hit, a row is drawn containing every
phenotype for which that ancestry has at least one significant window. All
p-values are shown faint grey as backdrop; significant windows are coloured by
phenotype. No threshold line is drawn.

Outputs (in plots/ukbb/):
  manhattan_significant_colored_BY.pdf / .png
  manhattan_significant_colored_raw.pdf / .png
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from statsmodels.stats.multitest import multipletests

GENERAL_OUTPUT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/ukbb'
PLOT_OUTPUT_FOLDER    = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/plots/ukbb'
EXCEL_PATH            = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx'

ANCESTRY_LIST = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']
META_COLS     = {'#CHROM', 'POS', 'ABS_POS', 'end_POS'}
ALPHA         = 0.20

MM_TO_INCH   = 1 / 25.4
FONT_DEFAULT = 6
FONT_TICKS   = 5
FONT_TITLE   = 7
LW_AXES      = 0.5
LABEL_CHROMS = {1, 5, 9, 13, 17, 21}
CHR_COLORS   = ['#333333', '#aaaaaa']   # snputils-style alternating grey per chromosome

mpl.rcParams.update({
    'font.size':       FONT_DEFAULT,
    'pdf.fonttype':    42,
    'ps.fonttype':     42,
    'savefig.dpi':     600,
    'font.family':     'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
})


def fp1_threshold_for(by_p, reject):
    by_p_sig = np.sort(by_p[reject])
    k_max = 0
    for k in range(1, len(by_p_sig) + 1):
        if by_p_sig[k - 1] * k <= 1:
            k_max = k
        else:
            break
    return by_p_sig[k_max - 1] if k_max > 0 else None


def load_excel_labels(path):
    try:
        df = pd.read_excel(path, sheet_name='first_batch', usecols='B:C')
        return df
    except Exception as e:
        print(f'  (excel labels unavailable: {e})')
        return None


def pheno_label(excel_df, pheno):
    if excel_df is None:
        label = pheno
    else:
        row = excel_df.loc[excel_df['ID'] == pheno, 'ID2']
        label = row.iloc[0] if not row.empty else pheno
    # drop the leading code (split on '_', remove the first token)
    parts = label.split('_')
    if len(parts) > 1:
        label = '_'.join(parts[1:])
    label = label.replace('_', ' ')
    # manual short names
    if 'gastro-oesophageal reflux' in label.lower():
        label = 'gastric reflux'
    return label


def compute_ancestry(ancestry):
    """Return (data, by_cols, sig_cols, hit_phenos) or None if no data."""
    data_file = os.path.join(GENERAL_OUTPUT_FOLDER, f'P_info_{ancestry}.tsv')
    if not os.path.exists(data_file):
        return None
    data = pd.read_csv(data_file, sep='\t')
    pheno_cols = [c for c in data.columns if c not in META_COLS]
    if not pheno_cols:
        return None

    by_cols, sig_cols, hit_phenos = {}, {}, []
    for pheno in pheno_cols:
        valid_idx = data[pheno].notna()
        vals = data.loc[valid_idx, pheno].values
        if len(vals) == 0:
            continue
        reject, by_p, _, _ = multipletests(vals, alpha=ALPHA, method='fdr_by')
        by_series = pd.Series(np.nan, index=data.index)
        by_series.loc[valid_idx] = by_p
        by_cols[pheno] = by_series

        sig_series = pd.Series(False, index=data.index)
        fp1 = fp1_threshold_for(by_p, reject)
        if fp1 is not None:
            sig_series.loc[valid_idx] = by_p <= fp1
            if sig_series.any():
                hit_phenos.append(pheno)
        sig_cols[pheno] = sig_series

    return data, by_cols, sig_cols, hit_phenos


# Per-subplot size, identical across all panels (= original 2x4 panel element:
# 183/3 mm wide, 78/2 mm tall). Both split figures reuse the exact same element.
SUBPLOT_W_MM = 183 / 3
SUBPLOT_H_MM = 78 / 2


def make_panel(results, color_map, excel_df, suffix, ancestry_subset, out_name, n_cols=2):
    # Keep the given subset order; only ancestries that actually have hits.
    anc_with_hits = [a for a in ancestry_subset if results.get(a) and results[a][3]]
    if not anc_with_hits:
        print(f'[{out_name} | {suffix}] no ancestry with hits — nothing to plot')
        return

    N_COLS = n_cols
    N_ROWS = int(np.ceil(len(anc_with_hits) / N_COLS))
    fig, axes = plt.subplots(
        N_ROWS, N_COLS,
        figsize=(N_COLS * SUBPLOT_W_MM * MM_TO_INCH, N_ROWS * SUBPLOT_H_MM * MM_TO_INCH),
        gridspec_kw={'hspace': 0.35, 'wspace': 0.25},
    )
    axes = np.atleast_2d(axes)

    for i, ancestry in enumerate(anc_with_hits):
        row, col = i // N_COLS, i % N_COLS
        ax = axes[row, col]
        is_bottom = (i + N_COLS) >= len(anc_with_hits)

        data, by_cols, sig_cols, hit_phenos = results[ancestry]
        plot_anc = ancestry

        # plotted values
        p_lookup = by_cols if suffix == 'BY' else {p: data[p] for p in by_cols}

        # ── snputils-style backdrop: ALL p-values, alternating grey per chromosome ──
        for chrom, grp in sorted(data.groupby('#CHROM')):
            xs, ys = [], []
            for pheno in p_lookup:
                v = p_lookup[pheno].loc[grp.index].values
                xs.append(grp['ABS_POS'].values)
                ys.append(v)
            xs = np.concatenate(xs); ys = np.concatenate(ys)
            ok = ~np.isnan(ys) & (ys > 0)
            if ok.any():
                ax.scatter(xs[ok], -np.log10(ys[ok]),
                           color=CHR_COLORS[int(chrom) % 2], s=0.8,
                           linewidths=0, rasterized=True, zorder=1)

        # ── coloured significant windows, by phenotype ──
        for pheno in hit_phenos:
            sig_mask = sig_cols[pheno]
            idx = data.index[sig_mask.values]
            if len(idx) == 0:
                continue
            xs = data.loc[idx, 'ABS_POS'].values
            ys = p_lookup[pheno].loc[idx].values
            ok = ~np.isnan(ys) & (ys > 0)
            if ok.any():
                ax.scatter(xs[ok], -np.log10(ys[ok]),
                           color=color_map[pheno], s=4, linewidths=0,
                           rasterized=True, zorder=3)

        # ── chromosome ticks: mark all, label only LABEL_CHROMS ──
        all_ticks, all_labels = [], []
        for chrom, grp in sorted(data.groupby('#CHROM')):
            all_ticks.append(grp['ABS_POS'].mean())
            all_labels.append(str(int(chrom)) if int(chrom) in LABEL_CHROMS else '')
        ax.set_xticks(all_ticks)
        if is_bottom:
            ax.set_xticklabels(all_labels, fontsize=FONT_TICKS)
            ax.set_xlabel('Chromosome', fontsize=FONT_DEFAULT, labelpad=2)
        else:
            ax.set_xticklabels([])

        if col == 0:
            ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=FONT_DEFAULT, labelpad=2)
        ax.tick_params(axis='y', labelsize=FONT_TICKS, width=LW_AXES, length=2)
        ax.tick_params(axis='x', width=LW_AXES, length=2)
        ax.set_title(plot_anc, fontsize=FONT_TITLE, fontweight='bold', pad=2)
        ax.set_xlim(data['ABS_POS'].min(), data['ABS_POS'].max())
        ax.set_ylim(bottom=0)
        ax.margins(y=0.05)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
        ax.spines['bottom'].set_linewidth(LW_AXES)
        ax.spines['left'].set_linewidth(LW_AXES)

    # hide unused subplots
    for j in range(len(anc_with_hits), N_ROWS * N_COLS):
        axes[j // N_COLS, j % N_COLS].set_visible(False)

    # ── legend: only phenotypes present in this subset (global colours) ──
    legend_phenos = sorted({p for a in anc_with_hits for p in results[a][3]})
    handles = [Line2D([0], [0], marker='o', linestyle='', markersize=3,
                      markerfacecolor=color_map[p], markeredgewidth=0,
                      label=pheno_label(excel_df, p))
               for p in legend_phenos]
    fig.legend(handles=handles, loc='upper center', ncol=min(3, len(handles)),
               frameon=False, fontsize=FONT_TICKS, bbox_to_anchor=(0.5, 0.0),
               columnspacing=1.0, handletextpad=0.4, labelspacing=0.3)

    out_base = os.path.join(PLOT_OUTPUT_FOLDER, f'manhattan_{out_name}_{suffix}')
    plt.savefig(f'{out_base}.pdf', dpi=600, bbox_inches='tight')
    plt.savefig(f'{out_base}.png', dpi=600, bbox_inches='tight')
    plt.close()
    print(f'[{out_name} | {suffix}] saved → {out_base}.pdf / .png  '
          f'({len(anc_with_hits)} ancestries, {len(legend_phenos)} phenotypes)')


def main():
    os.makedirs(PLOT_OUTPUT_FOLDER, exist_ok=True)
    excel_df = load_excel_labels(EXCEL_PATH)

    # ── PASS 1: compute everything + global set of hit phenotypes ──
    results = {}
    global_hit_phenos = set()
    for ancestry in ANCESTRY_LIST:
        res = compute_ancestry(ancestry)
        results[ancestry] = res
        if res and res[3]:
            global_hit_phenos.update(res[3])
            print(f'{ancestry}: hits in {sorted(res[3])}')
        else:
            print(f'{ancestry}: no hits')

    # ── global colour map: same phenotype → same colour across ancestries ──
    hit_list = sorted(global_hit_phenos)
    cmap = plt.get_cmap('tab20')
    color_map = {p: cmap(i % 20) for i, p in enumerate(hit_list)}
    print(f'\nGlobal hit phenotypes ({len(hit_list)}): {hit_list}\n')

    # ── PASS 2: two split panels, same element size ──
    # novelty = NAT/AMR + EAS ; existing = AFR, EUR, SAS, WAS
    NOVELTY  = ['NAT', 'EAS', 'EUR']
    EXISTING = ['AFR', 'SAS', 'WAS']
    for suffix in ('BY', 'raw'):
        make_panel(results, color_map, excel_df, suffix, NOVELTY,  out_name='novelty',  n_cols=2)
        make_panel(results, color_map, excel_df, suffix, EXISTING, out_name='existing', n_cols=2)


if __name__ == '__main__':
    main()
