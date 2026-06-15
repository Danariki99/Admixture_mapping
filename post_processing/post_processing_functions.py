import subprocess
import os
import time
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests
import requests
from collections import defaultdict
from snputils.visualization.manhattan_plot import manhattan_plot as su_manhattan


def parse_lambda(log_path):
    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        for line in f:
            m = re.search(r'lambda \(based on median chisq\) = ([\d.]+?)\.?\s', line)
            if m:
                return float(m.group(1))
    return None


def positions_extraction(input_file, output_folder):
    # Define the command
    command = f"awk -F'\t' '!/^##/ {{print $1\"\t\"$2}}' {input_file}"

    # Run the command
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    output, error = process.communicate()

    data = StringIO(output.decode())
    df = pd.read_csv(data, sep="\t")

    # Split the 'POS' column into 'start_pos' and 'end_pos'
    df[['POS', 'end_POS']] = df['POS'].str.split('_', expand=True)

    # Create output file path
    output_file = os.path.join(output_folder, 'positions.csv')
    df.to_csv(output_file, sep='\t', index=False)

    return output_file

def result_analysis(
    ancestry_list,
    phe_folder,
    general_file_ini,
    window_pos_file,
    general_output_file,
    plot_output_folder,
    general_output_folder,
    dataset_name=None,
    apply_fb_filter=False,
    fb_template=None,
    fb_root=None,
    fb_threshold=0.9,
    fb_chunksize=200000
):
    significance_threshold = 0.05
    lambda_min, lambda_max = 0.9, 1.1

    excel_df = pd.read_excel(
        '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx',
        sheet_name='first_batch',
        usecols="B:C"
    )

    significant_df = pd.DataFrame(columns=['#CHROM', 'POS', 'end_POS', 'ABS_POS', 'P', 'Phenotype', 'Ancestry'])

    chromosome_colors = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b',
        '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#aec7e8', '#ffbb78',
        '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#c7c7c7',
        '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31'
    ]

    timeout_map = {
        'chr6:31346445-31377047':   '6p21.33',
        'chr8:124070432-124092625': '8q24.13',
        'chr6:31905130-32007956':   '6p21.33',
        'chr6:32207393-32288190':   '6p21.32',
        'chr10:116036889-116139029':'10q25.3',
        'chr6:31428169-31435326':   '6p21.33',
        'chr9:85752837-85810910':   '9q21.32',
        'chr17:1820750-1925859':    '17p13.3',
    }

    os.makedirs(plot_output_folder, exist_ok=True)
    os.makedirs(general_output_folder, exist_ok=True)

    window_pos = pd.read_csv(window_pos_file, sep='\t')

    for ancestry in ancestry_list:
        print(ancestry)

        general_file = general_file_ini.replace('#', ancestry)
        output_file  = general_output_file.replace('#', ancestry)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        pheno_list = os.listdir(phe_folder)

        significant_dict = {}  # pheno -> sig DataFrame (only phenos with hits)
        valid_phenos = []
        data_all = None

        # ─────────────────────────────────────────────
        # LOOP OVER ALL PHENOTYPES (INDEPENDENT)
        # ─────────────────────────────────────────────
        for phe_file in pheno_list:
            pheno = phe_file.replace('.phe', '')
            current_file = general_file.replace('*', pheno)

            if not os.path.exists(current_file):
                print(f"[{ancestry}] Missing file: {pheno}")
                continue

            lam = parse_lambda(os.path.join(os.path.dirname(current_file), 'output.log'))

            if lam is None:
                print(f"[{ancestry}] Skipping {pheno}: λGC=None")
                continue

            if not (lambda_min <= lam <= lambda_max):
                print(f"[{ancestry}] Skipping {pheno}: λGC={lam}")
                continue

            print(f"[{ancestry}] Processing {pheno}: λGC={lam}")

            df = pd.read_table(current_file, sep="\t")
            df = df[['#CHROM', 'POS', 'P']]
            df = df[df['P'] != '.'].copy()
            df['P'] = pd.to_numeric(df['P'], errors='coerce')
            df = df.dropna(subset=['P'])
            df = df.rename(columns={'P': pheno})

            # ✔ ONLY GLOBAL WINDOWS (NO FIRST_PHENO BIAS)
            df = pd.merge(df, window_pos, on=['#CHROM', 'POS'], how='inner')

            # ABS position
            if data_all is None:
                max_pos = df.groupby('#CHROM')['POS'].max().cumsum().shift(fill_value=0)

            df['ABS_POS'] = df['POS'] + df['#CHROM'].map(max_pos)

            # BY correction
            reject, by_p, _, _ = multipletests(df[pheno].values, alpha=0.20, method='fdr_by')
            df[f'{pheno}_BY'] = by_p

            # FP=1 criterion on BY-significant items only (p_BY <= 0.20)
            # Filtering by reject first ensures threshold <= 0.20 by construction
            by_p_sig = np.sort(by_p[reject])
            k_max = 0
            for k in range(1, len(by_p_sig) + 1):
                if by_p_sig[k - 1] * k <= 1:
                    k_max = k
                else:
                    break
            if k_max > 0:
                fp1_threshold = by_p_sig[k_max - 1]
                sig = df[df[f'{pheno}_BY'] <= fp1_threshold].copy()
                print(f"  [{ancestry}] {pheno}: {len(sig)} significant windows (BY threshold={fp1_threshold:.4f}, FP attesi={len(sig)*fp1_threshold:.2f})")
            else:
                fp1_threshold = None
                sig = pd.DataFrame(columns=df.columns)

            if not sig.empty:
                sig = sig.copy()
                sig['Phenotype'] = pheno
                sig['Ancestry'] = ancestry
                significant_dict[pheno] = sig

            valid_phenos.append(pheno)

            # merge into global table
            keep_cols = ['#CHROM', 'POS', 'ABS_POS', pheno, f'{pheno}_BY']
            if data_all is None:
                data_all = df[keep_cols]
            else:
                data_all = pd.merge(
                    data_all,
                    df[keep_cols],
                    on=['#CHROM', 'POS', 'ABS_POS'],
                    how='outer'
                )

        if data_all is None:
            print(f"[{ancestry}] No valid phenotypes")
            continue

        # ─────────────────────────────
        # SAVE RAW DATA
        # ─────────────────────────────
        cols_raw = [c for c in data_all.columns if not c.endswith('_BY')]
        data_to_save = data_all[cols_raw].copy()
        if 'end_POS' not in data_to_save.columns:
            data_to_save = pd.merge(data_to_save, window_pos[['#CHROM', 'POS', 'end_POS']],
                                    on=['#CHROM', 'POS'], how='left')
        data_to_save.to_csv(output_file, index=False, sep='\t')

        # ─────────────────────────────
        # MANHATTAN PLOTS (BY + raw)
        # ─────────────────────────────
        plot_ancestry = 'AMR' if ancestry == 'NAT' else ancestry

        # max raw/BY p among accepted windows across all phenotypes — used as significance line
        max_fp1_raw = None
        max_fp1_by  = None
        for pheno, sig_df in significant_dict.items():
            t_raw = sig_df[pheno].max()
            t_by  = sig_df[f'{pheno}_BY'].max()
            if max_fp1_raw is None or t_raw > max_fp1_raw:
                max_fp1_raw = t_raw
            if max_fp1_by is None or t_by > max_fp1_by:
                max_fp1_by = t_by

        for plot_mode in ('BY', 'raw'):
            p_cols  = [f'{p}_BY' if plot_mode == 'BY' else p for p in valid_phenos]
            plot_df = data_all[['#CHROM', 'POS']].copy()
            plot_df['P'] = data_all[p_cols].min(axis=1, skipna=True)
            plot_df = plot_df.dropna(subset=['P']).copy()

            thresh = max_fp1_by if plot_mode == 'BY' else max_fp1_raw

            # snputils divides significance_threshold by len(df) internally,
            # so multiply back to place the line exactly at -log10(thresh)
            su_manhattan(
                plot_df,
                significance_threshold=thresh * len(plot_df) if thresh is not None else 1.0,
                line_color='r' if thresh is not None else 'none',
                figsize=(12, 6),
                title=f'Manhattan Plot {plot_ancestry} — {plot_mode}',
                save=False,
            )

            ax = plt.gca()
            texts = []

            for pheno in valid_phenos:
                sig_df = significant_dict.get(pheno, pd.DataFrame())
                if sig_df.empty:
                    continue

                y_col    = f'{pheno}_BY' if plot_mode == 'BY' else pheno
                max_row  = sig_df.loc[sig_df[pheno].idxmin()]
                value    = 0.2 if ancestry == 'SAS' else 0.6
                offset_x = np.random.uniform(-1e9, 1e9)
                offset_y = np.random.uniform(-value, value)

                end_pos = int(max_row['end_POS']) if 'end_POS' in max_row.index else int(max_row['POS'])
                name    = fetch_cytoband(f"chr{int(max_row['#CHROM'])}", int(max_row['POS']), end_pos)
                if name == 'Timeout':
                    coord = f"chr{int(max_row['#CHROM'])}:{int(max_row['POS'])}-{end_pos}"
                    name  = timeout_map.get(coord, coord)

                pheno_label     = excel_df.loc[excel_df['ID'] == pheno, 'ID2'].iloc[0].replace('_', ' ')
                annotation_text = f'{pheno_label}\n{name}'
                y_val = max_row[y_col]
                text  = ax.annotate(annotation_text,
                                    (max_row['ABS_POS'] + offset_x, -np.log10(y_val) + offset_y))
                texts.append(text)

                if plot_mode == 'BY':
                    sig_df_copy = sig_df.copy()
                    sig_df_copy['Phenotype'] = pheno
                    sig_df_copy['Ancestry']  = ancestry
                    significant_df = pd.concat([significant_df, sig_df_copy])

            adjust_text(texts, ax=ax)

            if plot_mode == 'BY':
                ax.set_ylabel('-log10(BY corrected p)')
            else:
                ax.set_ylabel('-log10(raw p)')

            plt.savefig(os.path.join(plot_output_folder, f'manhattan_plot_{plot_ancestry}_{plot_mode}.png'),
                        bbox_inches='tight')
            plt.close()

    significant_file = os.path.join(general_output_folder, 'significant_positions.tsv')

    if not significant_df.empty:
        # Require a contiguous run of ≥5 windows per hit (Ancestry × Phenotype)
        # (guaranteed by FDR≤0.2 + FP≤1 on count, this checks they are actually adjacent)
        def _max_contiguous_run(group):
            s = group.sort_values(['#CHROM', 'POS']).reset_index(drop=True)
            max_run = cur = 1
            for i in range(1, len(s)):
                if s.loc[i, '#CHROM'] == s.loc[i-1, '#CHROM'] and s.loc[i, 'POS'] == s.loc[i-1, 'end_POS']:
                    cur += 1
                    max_run = max(max_run, cur)
                else:
                    cur = 1
            return max_run

        valid_hits = []
        for (anc, phe), group in significant_df.groupby(['Ancestry', 'Phenotype']):
            run = _max_contiguous_run(group)
            if run < 5:
                print(f"  [{anc}] {phe}: skipped (max contiguous run = {run} < 5)")
            else:
                valid_hits.append({'Ancestry': anc, 'Phenotype': phe})

        if not valid_hits:
            print("\nNo hits with ≥5 contiguous windows after filtering.")
            return None

        significant_df = significant_df.merge(pd.DataFrame(valid_hits), on=['Ancestry', 'Phenotype'])

        significant_df.to_csv(significant_file, sep='\t', index=False)

        print(f"\n{'='*40}")
        print("Significant hits summary:")
        summary = (significant_df
                   .groupby(['Ancestry', 'Phenotype'])
                   .size()
                   .reset_index(name='n_windows')
                   .sort_values(['Ancestry', 'Phenotype']))
        for _, row in summary.iterrows():
            print(f"  {row['Ancestry']}  {row['Phenotype']}  {row['n_windows']} windows")
        print(f"{'='*40}\n")

        return significant_file

    print("\nNo significant hits found.")
    return None


def _normalize_chrom(value):
    value = str(value).strip()
    if value.lower().startswith('chr'):
        value = value[3:]
    return value


def _resolve_fb_path(fb_template, fb_root, dataset, ancestry, chrom):
    context = {
        'dataset': dataset or '',
        'ancestry': ancestry or '',
        'chrom': chrom
    }
    if fb_template:
        try:
            candidate = fb_template.format(**context)
        except KeyError:
            candidate = fb_template.format(chrom=chrom)
        if os.path.exists(candidate):
            return candidate
        gz_candidate = f"{candidate}.gz"
        if os.path.exists(gz_candidate):
            return gz_candidate
    if not fb_root:
        return None
    chrom_str = f"{chrom}"
    base_paths = {fb_root}
    if dataset:
        base_paths.add(os.path.join(fb_root, dataset))
    if ancestry:
        base_paths.add(os.path.join(fb_root, ancestry))
    if dataset and ancestry:
        base_paths.add(os.path.join(fb_root, dataset, ancestry))
        base_paths.add(os.path.join(fb_root, ancestry, dataset))
    patterns = []
    for base in base_paths:
        patterns.extend([
            os.path.join(base, f"chr{chrom_str}", "*.fb*"),
            os.path.join(base, f"chr_{chrom_str}", "*.fb*"),
            os.path.join(base, f"*chr{chrom_str}*.fb*"),
        ])
        if ancestry:
            patterns.append(os.path.join(base, f"*{ancestry}*chr{chrom_str}*.fb*"))
    for pattern in patterns:
        matches = sorted(glob.glob(pattern))
        if matches:
            return matches[0]
    return None


def _read_fb_confident_positions(fb_file, threshold, chunksize=200000):
    """Load FB probabilities with Polars, immediately convert to pandas for downstream logic."""
    del chunksize  # retained for signature compatibility
    confident = defaultdict(set)
    try:
        pdf = pd.read_csv(fb_file, sep='\t')
    except FileNotFoundError:
        return confident
    except pd.errors.EmptyDataError:
        return confident
    columns = pdf.columns
    lower_map = {col.lower(): col for col in columns}
    chrom_col = next((lower_map[candidate] for candidate in ['chromosome', '#chrom', 'chrom', 'chr'] if candidate in lower_map), None)
    pos_col = next((lower_map[candidate] for candidate in ['position', 'pos'] if candidate in lower_map), None)
    if chrom_col is None or pos_col is None:
        return confident
    meta_cols = {chrom_col, pos_col}
    meta_cols.update({col for col in columns if col and col.lower().startswith('genetic')})
    prob_cols = [col for col in columns if col not in meta_cols]
    if not prob_cols:
        return confident

    pdf[prob_cols] = pdf[prob_cols].apply(pd.to_numeric, errors='coerce')
    max_conf = pdf[prob_cols].max(axis=1)
    keep_mask = max_conf >= threshold
    if not keep_mask.any():
        return confident

    keep_df = pdf.loc[keep_mask, [chrom_col, pos_col]].copy()
    keep_df[pos_col] = pd.to_numeric(keep_df[pos_col], errors='coerce')
    keep_df = keep_df.dropna(subset=[pos_col])
    for chrom, pos in keep_df.itertuples(index=False):
        norm_chrom = _normalize_chrom(chrom)
        if norm_chrom:
            confident[norm_chrom].add(int(pos))
    return confident


def filter_windows_by_confidence(data_df, ancestry, dataset, fb_template=None, fb_root=None, threshold=0.9, chunksize=200000):
    summary = {
        'filtered': False,
        'threshold': threshold,
        'missing_chromosomes': [],
        'processed_chromosomes': [],
        'n_windows_before': len(data_df),
        'n_windows_after': len(data_df),
        'n_windows_removed': 0
    }
    if data_df.empty or (fb_template is None and fb_root is None):
        return data_df, summary
    keep_keys = set()
    for chrom in sorted(data_df['#CHROM'].unique()):
        fb_path = _resolve_fb_path(fb_template, fb_root, dataset, ancestry, chrom)
        if not fb_path or not os.path.exists(fb_path):
            summary['missing_chromosomes'].append(chrom)
            continue
        confident = _read_fb_confident_positions(fb_path, threshold, chunksize)
        chrom_key = _normalize_chrom(chrom)
        if chrom_key in confident:
            keep_keys.update({f"{chrom_key}:{pos}" for pos in confident[chrom_key]})
        summary['processed_chromosomes'].append(chrom)
    if not keep_keys:
        return data_df, summary
    chrom_series = data_df['#CHROM'].apply(_normalize_chrom)
    pos_series = pd.to_numeric(data_df['POS'], errors='coerce')
    valid_idx = chrom_series.notna() & pos_series.notna()
    composite = chrom_series[valid_idx].astype(str) + ':' + pos_series[valid_idx].astype(int).astype(str)
    final_mask = pd.Series(False, index=data_df.index)
    final_mask.loc[composite.index] = composite.isin(keep_keys)
    filtered_df = data_df.loc[final_mask].reset_index(drop=True)
    summary['filtered'] = True
    summary['n_windows_after'] = len(filtered_df)
    summary['n_windows_removed'] = summary['n_windows_before'] - summary['n_windows_after']
    return filtered_df, summary


def SNPs_extraction(input_file, output_dir, max_retries=5, retry_delay=120):
    r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'SNPs_gene_extraction.R')
    for attempt in range(1, max_retries + 1):
        if attempt > 1:
            print(f"  SNPs_extraction retry {attempt}/{max_retries} in {retry_delay}s...")
            time.sleep(retry_delay)
        result = subprocess.run(
            ['Rscript', r_script, '-i', input_file, '-o', output_dir],
            capture_output=True, text=True
        )
        output_lines = result.stdout.strip().split('\n')
        output_file_path = output_lines[-1] if output_lines and output_lines[-1] else ""
        if result.returncode == 0 and output_file_path and os.path.exists(output_file_path):
            return output_file_path
        print(f"  SNPs_extraction attempt {attempt} failed (code {result.returncode}):")
        print("  --- stderr ---")
        for line in result.stderr.strip().split('\n')[-30:]:
            print(f"    {line}")
        print("  --- stdout ---")
        for line in result.stdout.strip().split('\n')[-10:]:
            print(f"    {line}")
    raise RuntimeError(
        f"SNPs_extraction failed after {max_retries} attempts — Ensembl unavailable. Rerun when the server is up."
    )


# Create a function to find the closest SNP to the middle of a given window
def find_snps_in_window(window, snps_df):
    snps_in_window = snps_df[snps_df['pos'].between(window['POS'], window['end_POS'])].copy()
    if not snps_in_window.empty:
        snps_in_window['Phenotype'] = window['Phenotype']
        snps_in_window['Ancestry'] = window['Ancestry']
        snps_in_window['window_start'] = window['POS']
        snps_in_window['window_end'] = window['end_POS']
        snps_in_window['P'] = window['P']
        return snps_in_window[['CHR', 'pos', 'rfid', 'P', 'Phenotype', 'Ancestry', 'window_start', 'window_end']]
    return pd.DataFrame(columns=['CHR', 'pos', 'rfid', 'P', 'Phenotype', 'Ancestry', 'window_start', 'window_end'])

def associate_SNPs_to_windows(snps_file, window_file, output_folder):
    # Load the snps file
    snps_df = pd.read_csv(snps_file, sep='\t')
    snps_df.rename(columns={'chr_name': 'CHR'}, inplace=True)
    snps_df.rename(columns={'refsnp_id': 'rfid'}, inplace=True)
    snps_df.rename(columns={'chrom_start': 'pos'}, inplace=True)
    snps_df = snps_df.drop(columns=['allele'])

    window_df = pd.read_csv(window_file, sep='\t')

    # Convert the windows to intervals and calculate the mid point
    window_df['window'] = pd.IntervalIndex.from_arrays(window_df['POS'], window_df['end_POS'], closed='both')
    window_df['mid'] = window_df['window'].map(lambda x: x.mid)

    all_snps = pd.DataFrame()
    for _, window in window_df.iterrows():
        snps_in_window = find_snps_in_window(window, snps_df)
        all_snps = pd.concat([all_snps, snps_in_window])

    all_snps.columns = ['chr', 'pos', 'rfid', 'P', 'phenotype', 'ancestry', 'start', 'end']

    # Save the new file
    output_file = os.path.join(output_folder, "significant_SNPs_with_P_values.txt")
    all_snps.to_csv(output_file, sep='\t', index=False)
    return output_file

def FUMA_files_creation(snps_filename, output_folder):

    # create the necessary folders
    output_folder_snps = os.path.join(output_folder, 'snps')
    if not os.path.exists(output_folder_snps):
        os.makedirs(output_folder_snps)

    output_folder_wind = os.path.join(output_folder, 'wind')
    if not os.path.exists(output_folder_wind):
        os.makedirs(output_folder_wind)

    snps = pd.read_csv(snps_filename, sep='\t')

    ancestries = snps['ancestry'].unique()

    phenotypes = snps['phenotype'].unique()

    for ancestry in ancestries:

        for phenotype in phenotypes:

            snps_subset = snps[(snps['ancestry'] == ancestry) & (snps['phenotype'] == phenotype)]

            fuma_snps = snps_subset.drop(columns=['phenotype', 'ancestry', 'start', 'end'])
            fuma_snps.rename(columns={'chr' : 'CHR'}, inplace=True)
            fuma_wind = snps_subset.drop(columns=['pos', 'rfid', 'P', 'phenotype', 'ancestry'])
            fuma_wind = fuma_wind.drop_duplicates()

            output_file_snps = os.path.join(output_folder_snps, f'{ancestry}_{phenotype}_snps.txt')
            output_file_wind = os.path.join(output_folder_wind, f'{ancestry}_{phenotype}_wind.txt')
            if not snps_subset.empty:
                fuma_snps.to_csv(output_file_snps, sep='\t', index=False)
                fuma_wind.to_csv(output_file_wind, sep='\t', index=False)
            
    return output_folder_snps, output_folder_wind

def create_combined_manhattan(ancestry_list, plot_output_folder, suffix='BY'):
    ncols = 2
    nrows = (len(ancestry_list) + 1) // 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(24, 6 * nrows))
    for i, anc in enumerate(ancestry_list):
        row, col = i // ncols, i % ncols
        img_anc  = 'AMR' if anc == 'NAT' else anc
        img_path = os.path.join(plot_output_folder, f'manhattan_plot_{img_anc}_{suffix}.png')
        ax = axes[row, col]
        if os.path.exists(img_path):
            ax.imshow(plt.imread(img_path))
        ax.axis('off')
    for j in range(len(ancestry_list), nrows * ncols):
        axes[j // ncols, j % ncols].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_output_folder, f'manhattan_combined_{suffix}.png'), dpi=150, bbox_inches='tight')
    plt.close()


def create_panel_manhattan(ancestry_list, general_output_folder, plot_output_folder, suffix='BY'):
    MM_TO_INCH  = 1 / 25.4
    N_ROWS, N_COLS = 2, 4
    FONT_DEFAULT = 6
    FONT_TICKS   = 5
    FONT_TITLE   = 7
    LW_AXES      = 0.5
    LW_THRESH    = 0.5
    CHR_COLORS   = ['#333333', '#aaaaaa']

    with plt.rc_context({'pdf.fonttype': 42, 'ps.fonttype': 42}):
        fig, axes = plt.subplots(
            N_ROWS, N_COLS,
            figsize=(183 * MM_TO_INCH, 78 * MM_TO_INCH),
            gridspec_kw={'hspace': 0.35, 'wspace': 0.25},
        )

        for i, ancestry in enumerate(ancestry_list):
            row, col = i // N_COLS, i % N_COLS
            ax = axes[row, col]
            plot_anc = 'AMR' if ancestry == 'NAT' else ancestry

            data_file = os.path.join(general_output_folder, f'P_info_{ancestry}.tsv')
            if not os.path.exists(data_file):
                ax.set_visible(False)
                continue

            data = pd.read_csv(data_file, sep='\t')
            meta_cols = {'#CHROM', 'POS', 'ABS_POS', 'end_POS'}
            pheno_cols = [c for c in data.columns if c not in meta_cols]
            if not pheno_cols:
                ax.set_visible(False)
                continue

            # Recompute BY per phenotype + FP1 significance mask (per window)
            by_cols  = {}
            sig_cols = {}
            for pheno in pheno_cols:
                valid_idx = data[pheno].notna()
                vals = data.loc[valid_idx, pheno].values
                if len(vals) == 0:
                    continue
                reject, by_p, _, _ = multipletests(vals, alpha=0.20, method='fdr_by')
                by_series = pd.Series(np.nan, index=data.index)
                by_series.loc[valid_idx] = by_p
                by_cols[pheno] = by_series

                sig_series = pd.Series(False, index=data.index)
                by_p_sig = np.sort(by_p[reject])
                k_max = 0
                for k in range(1, len(by_p_sig) + 1):
                    if by_p_sig[k - 1] * k <= 1:
                        k_max = k
                    else:
                        break
                if k_max > 0:
                    fp1_threshold = by_p_sig[k_max - 1]
                    sig_series.loc[valid_idx] = by_p <= fp1_threshold
                sig_cols[pheno] = sig_series

            # Full p-value matrix: ALL phenotypes (one point per window x phenotype)
            cols = list(by_cols.keys())
            if suffix == 'BY':
                p_matrix = pd.DataFrame(by_cols)[cols]
            else:
                p_matrix = data[cols]
            sig_matrix = pd.DataFrame(sig_cols)[cols]

            # Adequate threshold: the strictest line such that EVERY point above it is
            # a genuine FP1-significant window. We may drop a few weak true hits, but
            # never draw a false one. = largest significant p still below the best
            # (smallest) p among all non-significant points.
            pv   = p_matrix.values.flatten()
            sigv = sig_matrix.values.flatten()
            good = ~np.isnan(pv) & (pv > 0)
            pv, sigv = pv[good], sigv[good]
            sig_p, nonsig_p = pv[sigv], pv[~sigv]

            max_thresh = None
            if sig_p.size > 0:
                if nonsig_p.size > 0:
                    floor_nonsig = nonsig_p.min()
                    eligible = sig_p[sig_p < floor_nonsig]
                    max_thresh = eligible.max() if eligible.size > 0 else None
                else:
                    max_thresh = sig_p.max()

            n_sig   = int(sig_p.size)
            n_shown = int((sig_p <= max_thresh).sum()) if max_thresh is not None else 0
            print(f'  [{plot_anc} | {suffix}] sig windows={n_sig}  shown above line={n_shown}  '
                  f'dropped={n_sig - n_shown}  '
                  f'threshold p={max_thresh if max_thresh is not None else "none"}')

            # Scatter all phenotype p-values per chromosome
            for chrom, grp in sorted(data.groupby('#CHROM')):
                abs_pos = grp['ABS_POS'].values
                sub     = p_matrix.loc[grp.index].values   # (n_windows, n_phenotypes)
                xs = np.repeat(abs_pos, sub.shape[1])
                ys = sub.flatten()
                valid = ~np.isnan(ys) & (ys > 0)
                if valid.any():
                    ax.scatter(
                        xs[valid],
                        -np.log10(ys[valid]),
                        color=CHR_COLORS[int(chrom) % 2],
                        s=0.8, linewidths=0, rasterized=True,
                    )

            # Significance line
            if max_thresh is not None:
                ax.axhline(y=-np.log10(max_thresh), color='r', linestyle='--',
                           linewidth=LW_THRESH, zorder=3)

            # Chromosome ticks — mark all, label only 1,5,9,13,17,21
            LABEL_CHROMS = {1, 5, 9, 13, 17, 21}
            all_ticks, all_labels = [], []
            for chrom, grp in sorted(data.groupby('#CHROM')):
                all_ticks.append(grp['ABS_POS'].mean())
                all_labels.append(str(int(chrom)) if int(chrom) in LABEL_CHROMS else '')
            ax.set_xticks(all_ticks)

            # Labels: ylabel only col==0, xlabel only last row
            if row == N_ROWS - 1:
                ax.set_xticklabels(all_labels, fontsize=FONT_TICKS)
                ax.set_xlabel('Chromosome', fontsize=FONT_DEFAULT, labelpad=2)
            else:
                ax.set_xticklabels([])
                ax.set_xlabel('')

            if col == 0:
                ax.set_ylabel(r'$-\log_{10}(p)$', fontsize=FONT_DEFAULT, labelpad=2)
            else:
                ax.set_ylabel('')
            ax.tick_params(axis='y', labelsize=FONT_TICKS, width=LW_AXES, length=2)
            # All chr ticks same length; unlabelled ones slightly shorter
            ax.tick_params(axis='x', width=LW_AXES, length=2)

            # Spines — restore bottom, hide top/right
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(True)
            ax.spines['bottom'].set_linewidth(LW_AXES)
            ax.spines['left'].set_linewidth(LW_AXES)

            ax.set_title(plot_anc, fontsize=FONT_TITLE, fontweight='bold', pad=2)
            ax.set_xlim(data['ABS_POS'].min(), data['ABS_POS'].max())
            # Remove auto-margin so points sit at y=0, same as snputils style
            ax.set_ylim(bottom=0)
            ax.margins(y=0.05)

        for j in range(len(ancestry_list), N_ROWS * N_COLS):
            axes[j // N_COLS, j % N_COLS].set_visible(False)

        out_base = os.path.join(plot_output_folder, f'manhattan_panel_{suffix}')
        plt.savefig(f'{out_base}.pdf', dpi=600, bbox_inches='tight')
        plt.savefig(f'{out_base}.png', dpi=600, bbox_inches='tight')
        plt.close()
        print(f'Panel saved → {out_base}.pdf / .png')


def fetch_cytoband(chromosome, start, end, genome="hg19", retries=3, retry_delay=5):
    import time
    url = f"https://api.genome.ucsc.edu/getData/track?genome={genome}&track=cytoBand&chrom={chromosome}&start={start}&end={end}"

    for attempt in range(1, retries + 1):
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            data = response.json()
            bands = data.get("cytoBand", [])
            if not bands:
                return "Unknown"
            chrom = bands[0]["chrom"].replace("chr", "")
            name  = bands[0]["name"]
            return f"{chrom}{name}"
        except requests.exceptions.ReadTimeout:
            print(f"Timeout on attempt {attempt}/{retries} for {chromosome}:{start}-{end}")
        except requests.exceptions.RequestException as e:
            print(f"Request error on attempt {attempt}/{retries}: {e}")
        except (KeyError, IndexError, ValueError) as e:
            print(f"Parse error on attempt {attempt}/{retries}: {e}")
        if attempt < retries:
            time.sleep(retry_delay)

    return "Timeout"

    
