import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import chi2

def load_pheno_table(path, sheet_name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Phenotype table not found: {path}")
    return pd.read_excel(path, sheet_name=sheet_name)


def clean_pheno_name(name):
    """ID2 -> display name: drop the HC#### code and the TTE_/AD_ token."""
    parts = str(name).split('_')
    if len(parts) > 1:
        parts = parts[1:]                       # drop the leading HC#### code
    if parts and parts[0] in ('TTE', 'AD'):
        parts = parts[1:]                       # drop the TTE_/AD_ token
    return '_'.join(parts)


# hits dropped from the paper (TTE_asthma)
EXCLUDED_PHENOS = {'HC1036'}


if __name__ == '__main__':
    #define all the data
    # FINAL hits of the paper (the old fine_mapping_ancestries_PCA_verbose folder
    # holds a different, obsolete hit set)
    HIT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_new'
    OUTPUT_FOLDER ='/private/groups/ioannidislab/smeriglio/out_cleaned_codes/QQ_plots'
    PHENO_TABLE = 'ukbb_v1.xlsx'
    PHENO_SHEET = 'first_batch'

    df_first_batch = load_pheno_table(PHENO_TABLE, PHENO_SHEET)

    p_file_template = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_*/#/output.#.glm.logistic.hybrid'
    # one QQ plot per (ancestry, phenotype): EUR_HC937 has two hit folders
    # (chr6 and chr10) but a single admixture-mapping result, so de-duplicate.
    hits_list = sorted({'_'.join(h.split('_')[:2]) for h in os.listdir(HIT_FOLDER)
                        if os.path.isdir(os.path.join(HIT_FOLDER, h))})
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    # Supplementary Figure 1 is the phenotype correlation matrix, so the QQ plots
    # start at 2 (the locuszoom script then auto-continues after the last QQ).
    counter = 1

    for hit in hits_list:
        ancestry = hit.split('_')[0]
        pheno = hit.split('_')[1]
        if pheno in EXCLUDED_PHENOS:               # TTE_asthma hit removed from the paper
            print(f'Skipping excluded hit: {hit}')
            continue
        counter += 1
        pheno_row = df_first_batch[df_first_batch['ID'] == pheno]['ID2']
        if not pheno_row.empty:
            pheno_name = pheno_row.iloc[0]
        else:
            pheno_name = 'Unknown'

        if pheno_name == 'HC1007_TTE_acute_upper_respiratory_infections_of_multiple_and_unspecified_sites':
            pheno_name = 'HC1007_TTE_acute_upper_respiratory_infections'
        print(pheno_name)
        p_file = p_file_template.replace('*', ancestry).replace('#', pheno)
        if not os.path.exists(p_file):
            print(f"Missing file for {ancestry} {pheno}: {p_file}")
            continue

        df = pd.read_csv(p_file, sep='\t')
        p_values = pd.to_numeric(df['P'], errors='coerce').dropna()
        if p_values.empty:
            print(f"No valid p-values for {ancestry} {pheno}")
            continue

        # Avoid zeros that would break the lambda calculation
        p_values = np.clip(np.sort(p_values.values), 1e-300, 1)
        n = len(p_values)
        expected = -np.log10((np.arange(1, n + 1)) / (n + 1))
        observed = -np.log10(p_values)

        # Genomic inflation factor (lambda GC)
        chi2_stats = chi2.isf(p_values, df=1)
        lambda_gc = np.median(chi2_stats) / 0.454936423119572  # median of chi2(1)
        print(f"λGC {ancestry} {pheno}: {lambda_gc:.3f} (n={n})")

        plt.figure(figsize=(6, 6))
        plt.scatter(expected, observed, s=10, color='steelblue', edgecolor='none')
        plt.plot([expected.min(), expected.max()], [expected.min(), expected.max()], color='firebrick', linestyle='--', linewidth=1)
        plt.xlabel('Expected -log10(p)')
        plt.ylabel('Observed -log10(p)')
        plt.title(f'{ancestry} {clean_pheno_name(pheno_name)}\nλGC = {lambda_gc:.3f}')
        plt.tight_layout()
        if '/' in pheno_name:
            pheno_name = pheno_name.replace('/', '_')
        output_path = os.path.join(
            OUTPUT_FOLDER, f"Supplementary_Figure_{counter}_{ancestry}_{pheno}.png")
        plt.savefig(output_path, dpi=300)
        plt.close()

    # No file output for lambda; printed above per ancestry/pheno
