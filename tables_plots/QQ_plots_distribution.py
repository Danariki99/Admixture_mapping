import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

DATASET      = 'ukbb'
HIT_FOLDER   = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/{DATASET}/wind'
GLM_TEMPLATE = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/{DATASET}/output_ancestry_{{ancestry}}/{{pheno}}/output.{{pheno}}.glm.logistic.hybrid'
OUTPUT_FOLDER = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/QQ_plots'
PHENO_TABLE  = '/private/home/rsmerigl/codes/cleaned_codes/Admixture_mapping/tables_plots/ukbb_v1.xlsx'


def parse_lambda(log_path):
    with open(log_path) as f:
        for line in f:
            m = re.search(r'lambda \(based on median chisq\) = ([\d.]+?)\.?\s', line)
            if m:
                return float(m.group(1))
    return None


if __name__ == '__main__':
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    excel_df = pd.read_excel(PHENO_TABLE, sheet_name='first_batch', usecols='B:C')

    hits = sorted(os.listdir(HIT_FOLDER))
    counter = 0

    for hit_file in hits:
        parts    = hit_file.replace('_wind.txt', '').split('_')
        ancestry = parts[0]
        pheno    = parts[1]

        glm_file = GLM_TEMPLATE.format(ancestry=ancestry, pheno=pheno)
        if not os.path.exists(glm_file):
            print(f'[SKIP] Missing GLM file: {ancestry} {pheno}')
            continue

        log_file  = os.path.join(os.path.dirname(glm_file), 'output.log')
        lambda_gc = parse_lambda(log_file) if os.path.exists(log_file) else None

        df = pd.read_csv(glm_file, sep='\t')
        p_values = pd.to_numeric(df['P'], errors='coerce').dropna()
        p_values = p_values[p_values > 0]
        if p_values.empty:
            print(f'[SKIP] No valid p-values: {ancestry} {pheno}')
            continue

        p_sorted  = np.clip(np.sort(p_values.values), 1e-300, 1)
        n         = len(p_sorted)
        expected  = -np.log10(np.arange(1, n + 1) / (n + 1))
        observed  = -np.log10(p_sorted)

        row        = excel_df.loc[excel_df['ID'] == pheno, 'ID2']
        pheno_name = row.iloc[0] if not row.empty else pheno
        short_name = '_'.join(pheno_name.split('_')[1:])
        if len(short_name) > 60:
            short_name = short_name[:57] + '...'

        lambda_str = f'{lambda_gc:.3f}' if lambda_gc is not None else 'N/A'
        counter += 1
        print(f'[{counter}] {ancestry} {pheno}: λGC={lambda_str}')

        plt.figure(figsize=(6, 6))
        plt.scatter(expected, observed, s=10, color='steelblue', edgecolor='none')
        plt.plot([expected.min(), expected.max()], [expected.min(), expected.max()],
                 color='firebrick', linestyle='--', linewidth=1)
        plt.xlabel('Expected -log10(p)')
        plt.ylabel('Observed -log10(p)')
        plt.title(f'{ancestry}  {short_name}\nλGC = {lambda_str}')
        plt.tight_layout()

        output_path = os.path.join(OUTPUT_FOLDER, f'Supplementary_Figure_{counter}_{ancestry}_{pheno}.png')
        plt.savefig(output_path, dpi=300)
        plt.close()

    print(f'\nDone. {counter} QQ plots saved to {OUTPUT_FOLDER}')
