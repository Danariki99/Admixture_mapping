import os
import re
import time
import requests
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import multipletests

DATASET      = 'ukbb'
HIT_FOLDER   = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{DATASET}/wind'
GLM_TEMPLATE = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/{DATASET}/output_ancestry_{{ancestry}}/{{pheno}}/output.{{pheno}}.glm.logistic.hybrid'
PHENO_TABLE  = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx'
OUTPUT_FILE  = 'table_hits_wind.xlsx'

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


def parse_lambda(log_path):
    if not os.path.exists(log_path):
        return None
    with open(log_path) as f:
        for line in f:
            m = re.search(r'lambda \(based on median chisq\) = ([\d.]+?)\.?\s', line)
            if m:
                return float(m.group(1))
    return None


def fetch_cytoband(chromosome, start, end, genome='hg19', retries=3, retry_delay=5):
    url = (f'https://api.genome.ucsc.edu/getData/track?genome={genome}'
           f'&track=cytoBand&chrom={chromosome}&start={start}&end={end}')
    for attempt in range(1, retries + 1):
        try:
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            bands = resp.json().get('cytoBand', [])
            if not bands:
                return 'Unknown'
            return f"{bands[0]['chrom'].replace('chr', '')}{bands[0]['name']}"
        except requests.exceptions.ReadTimeout:
            print(f'    Timeout attempt {attempt}/{retries}')
        except Exception as e:
            print(f'    Error attempt {attempt}/{retries}: {e}')
        if attempt < retries:
            time.sleep(retry_delay)
    return 'Timeout'


excel_df = pd.read_excel(PHENO_TABLE, sheet_name='first_batch', usecols='B:C')

rows = []
for hit_file in sorted(os.listdir(HIT_FOLDER)):
    parts    = hit_file.replace('_wind.txt', '').split('_')
    ancestry = parts[0]
    pheno    = parts[1]

    wind_df = pd.read_csv(os.path.join(HIT_FOLDER, hit_file), sep='\t')

    glm_file = GLM_TEMPLATE.format(ancestry=ancestry, pheno=pheno)
    if not os.path.exists(glm_file):
        print(f'[SKIP] {ancestry} {pheno}: GLM file not found')
        continue

    glm = pd.read_csv(glm_file, sep='\t')
    glm = glm[glm['P'] != '.'].copy()
    glm['P']      = pd.to_numeric(glm['P'],      errors='coerce')
    glm['#CHROM'] = pd.to_numeric(glm['#CHROM'], errors='coerce')
    glm['POS']    = pd.to_numeric(glm['POS'],    errors='coerce')
    glm = glm.dropna(subset=['P'])

    # Compute BY correction and FP=1 threshold genome-wide (same logic as post_processing)
    _, by_p, _, _ = multipletests(glm['P'].values, alpha=0.20, method='fdr_by')
    by_sorted = np.sort(by_p)
    k_max = 0
    for k in range(1, len(by_sorted) + 1):
        if by_sorted[k - 1] * k <= 1:
            k_max = k
        else:
            break
    fdr_threshold = by_sorted[k_max - 1] if (k_max > 0 and by_sorted[k_max - 1] < 1.0) else None
    glm['BY'] = by_p
    sig_by = glm.loc[glm['BY'] <= fdr_threshold, 'BY'] if fdr_threshold is not None else pd.Series(dtype=float)
    fdr_min = sig_by.min() if not sig_by.empty else None

    log_file  = os.path.join(os.path.dirname(glm_file), 'output.log')
    lambda_gc = parse_lambda(log_file)

    pheno_row  = excel_df.loc[excel_df['ID'] == pheno, 'ID2']
    pheno_name = pheno_row.iloc[0] if not pheno_row.empty else pheno

    # One row per chromosome (handles hits with significant windows on multiple chromosomes)
    for chrom, wind_chr_df in wind_df.groupby('chr'):
        merged = pd.merge(
            glm, wind_chr_df,
            left_on=['#CHROM', 'POS'], right_on=['chr', 'start'],
            how='inner'
        )

        if merged.empty:
            print(f'[WARN] {ancestry} {pheno} chr{chrom}: no overlap between GLM and wind file, using global min p on chr')
            glm_chr = glm[glm['#CHROM'] == chrom]
            best    = glm_chr.loc[glm_chr['P'].idxmin()]
            chr_val = int(best['#CHROM'])
            start   = int(wind_chr_df['start'].min())
            end     = int(wind_chr_df['end'].max())
        else:
            best    = merged.loc[merged['P'].idxmin()]
            chr_val = int(merged['chr'].iloc[0])
            start   = int(merged['start'].min())
            end     = int(merged['end'].max())

        p    = best['P']
        oddr = best['OR']
        l95  = best['L95']
        u95  = best['U95']

        print(f'  {ancestry} {pheno}: chr{chr_val}:{start}-{end}  fetching cytoband...')
        cytoband = fetch_cytoband(f'chr{chr_val}', start, end)
        if cytoband == 'Timeout':
            cytoband = TIMEOUT_MAP.get(f'chr{chr_val}:{start}-{end}', f'chr{chr_val}:{start}-{end}')

        rows.append({
            'Phenotype':     pheno_name,
            'CytoBand':      cytoband,
            'Ancestry':      ancestry,
            'OR (CI 95%)':   f'{oddr} ({l95}, {u95})',
            'p value':       p,
            'FDR_min':       fdr_min,
            'FDR_max':       fdr_threshold,
            'lambda_GC':     lambda_gc,
            'chr':           chr_val,
            'START':         start,
            'END':           end,
            'n_sig_windows': len(wind_chr_df),
        })
        lambda_str = f'{lambda_gc:.3f}' if lambda_gc is not None else 'N/A'
        print(f'    → {cytoband}  p={p:.2e}  λGC={lambda_str}  n_windows={len(wind_chr_df)}')

out = pd.DataFrame(rows)
out.to_excel(OUTPUT_FILE, index=False)
print(f'\nSaved {len(out)} rows to {OUTPUT_FILE}')
