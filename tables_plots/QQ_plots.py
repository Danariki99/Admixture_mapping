import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def load_pheno_table(path, sheet_name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Phenotype table not found: {path}")
    return pd.read_excel(path, sheet_name=sheet_name)


if __name__ == '__main__':
    #define all the data
    HIT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'
    OUTPUT_FOLDER ='/private/groups/ioannidislab/smeriglio/out_cleaned_codes/QQ_plots'
    PHENO_TABLE = 'ukbb_v1.xlsx'
    PHENO_SHEET = 'first_batch'

    df_first_batch = load_pheno_table(PHENO_TABLE, PHENO_SHEET)

    p_file_template = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_*/#/output.#.glm.logistic.hybrid'
    hits_list = os.listdir(HIT_FOLDER)
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    for hit in hits_list:
        ancestry = hit.split('_')[0]
        pheno = hit.split('_')[1]
        pheno_row = df_first_batch[df_first_batch['ID'] == pheno]['ID2']
        if not pheno_row.empty:
            pheno_name = pheno_row.iloc[0]
        else:
            pheno_name = 'Unknown'
        p_file = p_file_template.replace('*', ancestry).replace('#', pheno)
        if not os.path.exists(p_file):
            print(f"Missing file for {ancestry} {pheno}: {p_file}")
            continue

        df = pd.read_csv(p_file, sep='\t')
        p_values = pd.to_numeric(df['P'], errors='coerce').dropna()
        if p_values.empty:
            print(f"No valid p-values for {ancestry} {pheno}")
            continue

        p_values = np.sort(p_values)
        n = len(p_values)
        expected = -np.log10((np.arange(1, n + 1)) / (n + 1))
        observed = -np.log10(p_values)

        plt.figure(figsize=(6, 6))
        plt.scatter(expected, observed, s=10, color='steelblue', edgecolor='none')
        plt.plot([expected.min(), expected.max()], [expected.min(), expected.max()], color='firebrick', linestyle='--', linewidth=1)
        plt.xlabel('Expected -log10(p)')
        plt.ylabel('Observed -log10(p)')
        plt.title(f'{ancestry} {pheno_name}')
        plt.tight_layout()
        if '/' in pheno_name:
            pheno_name = pheno_name.replace('/', '_')
        output_path = os.path.join(OUTPUT_FOLDER, f'{ancestry}_{pheno_name}_QQ_plot.png')
        plt.savefig(output_path, dpi=300)
        plt.close()
