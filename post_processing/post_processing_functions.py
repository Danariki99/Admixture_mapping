import subprocess
import os
import glob
import re
import pandas as pd
import polars as pl
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO
from adjustText import adjust_text
from statsmodels.stats.multitest import multipletests
import requests
from collections import defaultdict


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
            _, by_p, _, _ = multipletests(df[pheno].values, alpha=0.05, method='fdr_by')
            df[f'{pheno}_BY'] = by_p

            sig = df[df[f'{pheno}_BY'] < significance_threshold].copy()

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
        chrom_positions = data_all.groupby('#CHROM')['ABS_POS'].mean().tolist()
        chrom_labels    = sorted(data_all['#CHROM'].unique())
        n_windows       = len(data_all)
        bonf_thresh     = significance_threshold / n_windows
        plot_ancestry   = 'AMR' if ancestry == 'NAT' else ancestry

        # empirical BY-equivalent threshold (max raw p that passed BY)
        empirical_thresh = None
        for pheno, sig_df in significant_dict.items():
            t = sig_df[pheno].max()
            if empirical_thresh is None or t > empirical_thresh:
                empirical_thresh = t

        for plot_mode in ('BY', 'raw'):
            plt.figure(figsize=(12, 6))
            texts = []

            for pheno in valid_phenos:
                y_col = f'{pheno}_BY' if plot_mode == 'BY' else pheno
                for chrom, chrom_data in data_all.groupby('#CHROM'):
                    plt.scatter(chrom_data['ABS_POS'], -np.log10(chrom_data[y_col]),
                                color=chromosome_colors[chrom - 1], s=7)

                sig_df = significant_dict.get(pheno, pd.DataFrame())
                if sig_df.empty:
                    continue

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
                text  = plt.annotate(annotation_text,
                                     (max_row['ABS_POS'] + offset_x, -np.log10(y_val) + offset_y))
                texts.append(text)

                if plot_mode == 'BY':
                    sig_df = sig_df.copy()
                    sig_df['Phenotype'] = pheno
                    sig_df['Ancestry']  = ancestry
                    significant_df = pd.concat([significant_df, sig_df])

            adjust_text(texts)

            if plot_mode == 'BY':
                plt.axhline(y=-np.log10(significance_threshold), color='r', linestyle='--',
                            label='BY threshold (FDR 0.05)')
                plt.ylabel('-log10(BY corrected p)')
            else:
                plt.axhline(y=-np.log10(bonf_thresh), color='b', linestyle='--',
                            label=f'Bonferroni (0.05/{n_windows})')
                if empirical_thresh is not None:
                    plt.axhline(y=-np.log10(empirical_thresh), color='r', linestyle='--',
                                label='Empirical BY equivalent')
                plt.ylabel('-log10(raw p)')

            plt.xticks(chrom_positions, chrom_labels)
            plt.xlabel('Chromosome')
            plt.title(f'Manhattan Plot {plot_ancestry} — {plot_mode}')
            plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=2)
            plt.savefig(os.path.join(plot_output_folder, f'manhattan_plot_{plot_ancestry}_{plot_mode}.png'),
                        bbox_inches='tight')
            plt.close()

    significant_file = os.path.join(general_output_folder, 'significant_positions.tsv')

    if not significant_df.empty:
        significant_df.to_csv(significant_file, sep='\t', index=False)
        return significant_file

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
        pl_df = pl.read_csv(fb_file, separator='\t', has_header=True)
    except FileNotFoundError:
        return confident
    except pl.exceptions.NoDataError:
        return confident

    pdf = pl_df.to_pandas()
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


def SNPs_extraction(input_file, output_dir):
    # Call the R script using subprocess
    current_dir = os.getcwd()
    
    result = subprocess.run(['Rscript', os.path.join(current_dir, 'SNPs_gene_extraction.R'), '-i', input_file, '-o', output_dir], capture_output=True, text=True)

    # Split the stdout into lines and get the last line (the output file path)
    output_lines = result.stdout.strip().split('\n')
    output_file_path = output_lines[-1] if output_lines else ""

    return output_file_path


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

    
