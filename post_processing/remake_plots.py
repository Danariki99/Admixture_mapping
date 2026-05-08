import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests
from adjustText import adjust_text
from post_processing_functions import fetch_cytoband

EXCEL_PATH = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx'

CHROMOSOME_COLORS = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78',
    '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7',
    '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
]

TIMEOUT_MAP = {
    'chr6:31346445-31377047':    '6p21.33',
    'chr8:124070432-124092625':  '8q24.13',
    'chr6:31905130-32007956':    '6p21.33',
    'chr6:32207393-32288190':    '6p21.32',
    'chr10:116036889-116139029': '10q25.3',
    'chr6:31428169-31435326':    '6p21.33',
    'chr9:85752837-85810910':    '9q21.32',
    'chr17:1820750-1925859':     '17p13.3',
}


def load_excel(path):
    return pd.read_excel(path, sheet_name='first_batch', usecols='B:C')


def get_pheno_label(excel_df, pheno):
    row = excel_df.loc[excel_df['ID'] == pheno, 'ID2']
    return row.iloc[0].replace('_', ' ') if not row.empty else pheno


def get_cytoband_label(max_row):
    chrom = int(max_row['#CHROM'])
    pos   = int(max_row['POS'])
    end   = int(max_row['end_POS']) if 'end_POS' in max_row.index else pos
    name  = fetch_cytoband(f'chr{chrom}', pos, end)
    if name == 'Timeout':
        coord = f'chr{chrom}:{pos}-{end}'
        name  = TIMEOUT_MAP.get(coord, coord)
    return name


def make_manhattan(data, pheno_cols, by_pvals, sig_dict,
                   chrom_positions, chrom_labels,
                   ancestry, title, ylabel, threshold_line,
                   output_path, excel_df):
    plt.figure(figsize=(12, 6))
    texts = []

    for pheno in pheno_cols:
        y_col = by_pvals[pheno] if by_pvals else pheno
        for chrom, cdata in data.groupby('#CHROM'):
            y = cdata[y_col] if isinstance(y_col, str) else y_col.loc[cdata.index]
            plt.scatter(cdata['ABS_POS'], -np.log10(y),
                        color=CHROMOSOME_COLORS[int(chrom) - 1], s=7)

        if pheno in sig_dict:
            max_row  = sig_dict[pheno]
            value    = 0.2 if ancestry == 'SAS' else 0.6
            offset_x = np.random.uniform(-1e9, 1e9)
            offset_y = np.random.uniform(-value, value)
            label    = get_cytoband_label(max_row)
            pheno_lbl = get_pheno_label(excel_df, pheno)
            by_y = data.loc[max_row.name, pheno] if pheno in data.columns else max_row[pheno]
            text = plt.annotate(
                f'{pheno_lbl}\n{label}',
                (max_row['ABS_POS'] + offset_x, -np.log10(by_y) + offset_y)
            )
            texts.append(text)

    adjust_text(texts)
    for y_val, color, label in threshold_line:
        plt.axhline(y=y_val, color=color, linestyle='--', linewidth=1, label=label)
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=len(threshold_line))
    plt.xticks(chrom_positions, chrom_labels)
    plt.xlabel('Chromosome')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', dpi=200)
    plt.close()
    print(f'Saved: {output_path}')


def process_ancestry(ancestry, tsv_path, plot_folder, excel_df):
    print(f'\n=== {ancestry} ===')
    data = pd.read_csv(tsv_path, sep='\t')

    # add end_POS from positions.csv if missing
    if 'end_POS' not in data.columns:
        pos_file = os.path.join(os.path.dirname(tsv_path), 'positions.csv')
        if os.path.exists(pos_file):
            window_pos = pd.read_csv(pos_file, sep='\t')[['#CHROM', 'POS', 'end_POS']]
            window_pos['#CHROM'] = window_pos['#CHROM'].astype(data['#CHROM'].dtype)
            window_pos['POS']    = window_pos['POS'].astype(data['POS'].dtype)
            data = pd.merge(data, window_pos, on=['#CHROM', 'POS'], how='left')

    all_meta   = {'#CHROM', 'POS', 'end_POS', 'ABS_POS'}
    meta_cols  = [c for c in ['#CHROM', 'POS', 'end_POS', 'ABS_POS'] if c in data.columns]
    pheno_cols = [c for c in data.columns if c not in all_meta]

    if not pheno_cols:
        print(f'  No phenotype columns found, skipping.')
        return False

    chrom_positions = data.groupby('#CHROM')['ABS_POS'].mean().tolist()
    chrom_labels    = sorted(data['#CHROM'].unique())

    # ── compute BY correction per pheno ──
    by_series = {}   # pheno -> Series of BY p-values (same index as data)
    sig_by    = {}   # pheno -> row with min raw p (for BY plot annotation)
    sig_raw   = {}   # pheno -> row with min raw p (for raw plot annotation)
    empirical_thresholds = {}  # pheno -> max raw p that passed BY

    for pheno in pheno_cols:
        p = data[pheno].copy()
        valid_mask = p.notna() & (p > 0)
        if valid_mask.sum() == 0:
            continue

        by_p_vals              = np.ones(len(data))
        _, corrected, _, _     = multipletests(p[valid_mask].values, alpha=0.05, method='fdr_by')
        by_p_vals[valid_mask]  = corrected
        by_series[pheno]       = pd.Series(by_p_vals, index=data.index)

        sig_mask = by_p_vals < 0.05
        if sig_mask.any():
            best_idx        = p[valid_mask][corrected < 0.05].idxmin()
            row             = data.loc[best_idx].copy()
            sig_by[pheno]   = row
            sig_raw[pheno]  = row
            empirical_thresholds[pheno] = p[valid_mask][corrected < 0.05].max()
            print(f'  Hit: {pheno} — best window chr{int(row["#CHROM"])}:{int(row["POS"])} '
                  f'(raw p={row[pheno]:.2e}, empirical threshold≈{empirical_thresholds[pheno]:.2e})')

    if not by_series:
        print('  No valid phenotypes after BY computation.')
        return False

    # ── build a merged DataFrame with BY columns for scatter ──
    data_by = data[meta_cols].copy()
    for pheno, s in by_series.items():
        data_by[pheno] = s

    # global empirical threshold (max across phenos, for display)
    if empirical_thresholds:
        global_empirical = max(empirical_thresholds.values())
    else:
        global_empirical = None

    n_windows    = len(data)
    bonf_thresh  = 0.05 / n_windows

    plot_anc = 'AMR' if ancestry == 'NAT' else ancestry

    valid_phenos   = [p for p in pheno_cols if p in by_series]
    threshold_by   = [(-np.log10(0.05), 'red', 'BY threshold (FDR 0.05)')]
    threshold_raw  = [(-np.log10(bonf_thresh), 'blue', f'Bonferroni (0.05/{n_windows})')]
    if global_empirical is not None:
        threshold_raw.append((-np.log10(global_empirical), 'red', 'Empirical BY equivalent'))

    # ── Individual plots ──
    make_manhattan(
        data=data_by, pheno_cols=valid_phenos, by_pvals=None, sig_dict=sig_by,
        chrom_positions=chrom_positions, chrom_labels=chrom_labels,
        ancestry=ancestry, title=f'Manhattan Plot {plot_anc} — BY corrected p',
        ylabel='-log10(BY corrected p)', threshold_line=threshold_by,
        output_path=os.path.join(plot_folder, f'manhattan_plot_{plot_anc}_BY.png'),
        excel_df=excel_df,
    )
    make_manhattan(
        data=data, pheno_cols=valid_phenos, by_pvals=None, sig_dict=sig_raw,
        chrom_positions=chrom_positions, chrom_labels=chrom_labels,
        ancestry=ancestry, title=f'Manhattan Plot {plot_anc} — raw p',
        ylabel='-log10(raw p)', threshold_line=threshold_raw,
        output_path=os.path.join(plot_folder, f'manhattan_plot_{plot_anc}_raw.png'),
        excel_df=excel_df,
    )
    return True



def make_combined(ancestries, plot_folder, suffix, output):
    ncols = 2
    nrows = (len(ancestries) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(24, 6 * nrows))
    axes = axes.flatten()
    for i, anc in enumerate(ancestries):
        plot_anc = 'AMR' if anc == 'NAT' else anc
        img_path = os.path.join(plot_folder, f'manhattan_plot_{plot_anc}_{suffix}.png')
        if os.path.exists(img_path):
            axes[i].imshow(plt.imread(img_path))
        axes[i].axis('off')
    for j in range(len(ancestries), len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    plt.savefig(output, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'Saved: {output}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dataset', help='Dataset name (e.g. ukbb)')
    args = parser.parse_args()

    dataset      = args.dataset
    tsv_folder   = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/post_processing_files/{dataset}'
    plot_folder  = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/plots/{dataset}'
    os.makedirs(plot_folder, exist_ok=True)

    excel_df     = load_excel(EXCEL_PATH)
    ancestries   = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']

    processed = []
    for anc in ancestries:
        tsv = os.path.join(tsv_folder, f'P_info_{anc}.tsv')
        if not os.path.exists(tsv):
            print(f'Skipping {anc}: {tsv} not found')
            continue
        if process_ancestry(anc, tsv, plot_folder, excel_df):
            processed.append(anc)

    # ── Combined BY ──
    make_combined(processed, plot_folder, suffix='BY',
                  output=os.path.join(plot_folder, 'manhattan_combined_BY.png'))

    # ── Combined raw ──
    make_combined(processed, plot_folder, suffix='raw',
                  output=os.path.join(plot_folder, 'manhattan_combined_raw.png'))


if __name__ == '__main__':
    main()
