import os
import pandas as pd

HIT_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'
MODELS_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/models'
PROBS_FOLDER = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/probs'

ANCESTRY_LIST = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS', 'NAT']
SCENARIOS = {
    'ancestry': 'Local Ancestry',
    'add': 'Genotype (ADD)',
    'environment': 'Environmental Covariates',
}

PHENO_TABLE = 'ukbb_v1.xlsx'
PHENO_SHEET = 'first_batch'

TABLE_COLUMNS = [
    'Phenotype',
    'ancestry tested',
    'ancestry of the population',
    'Start Position of the reduced significant region',
    'End Position of the reduced significant region',
    'Number of significant SNPs',
    'ID',
    'allele',
    'OR (CI = 95%)',
    'p value',
    'chr',
    'Delta_P_mean_ancestry',
    'Delta_P_median_ancestry',
    'Delta_P_std_ancestry',
    'Delta_P_mean_add',
    'Delta_P_median_add',
    'Delta_P_std_add',
    'Delta_P_mean_environment',
    'Delta_P_median_environment',
    'Delta_P_std_environment',
    'Delta_P_samples',
]


def load_pheno_table(path, sheet_name):
    if not os.path.exists(path):
        raise FileNotFoundError(f"Phenotype table not found: {path}")
    return pd.read_excel(path, sheet_name=sheet_name)


def build_table():
    df_first_batch = load_pheno_table(PHENO_TABLE, PHENO_SHEET)
    rows = []

    if not os.path.isdir(HIT_FOLDER):
        raise FileNotFoundError(f"Hit folder not found: {HIT_FOLDER}")

    hits_list = os.listdir(HIT_FOLDER)

    for hit in hits_list:
        print(hit)
        hit_parts = hit.split('_')
        if len(hit_parts) < 3:
            continue

        imp_ancestry = hit_parts[0]
        pheno = hit_parts[1]
        chrom = hit_parts[-1]

        pheno_row = df_first_batch[df_first_batch['ID'] == pheno]['ID2']
        if not pheno_row.empty:
            pheno_name = pheno_row.iloc[0]
        else:
            pheno_name = 'Unknown'

        glm_base_path = os.path.join(HIT_FOLDER, hit)
        model_base_path = os.path.join(MODELS_FOLDER, hit)

        for ancestry in ANCESTRY_LIST:
            glm_path = os.path.join(glm_base_path, f"{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid")
            if not os.path.exists(glm_path):
                continue

            glm_df = pd.read_csv(glm_path, sep='\t')

            model_folder = os.path.join(model_base_path, ancestry)
            if not os.path.isdir(model_folder):
                continue

            snps_list = [
                filename.replace('.csv', '')
                for filename in os.listdir(model_folder)
                if filename.endswith('.csv')
            ]
            if not snps_list:
                continue

            filtered_df = glm_df[(glm_df['TEST'] == imp_ancestry) & (glm_df['ID'].isin(snps_list))]
            filtered_df = filtered_df.dropna(subset=['P'])
            if filtered_df.empty:
                continue

            start_pos = filtered_df['POS'].min()
            end_pos = filtered_df['POS'].max()
            most_significant_snp = filtered_df.loc[filtered_df['P'].idxmin()]['ID']

            add_row = glm_df[(glm_df['ID'] == most_significant_snp) & (glm_df['TEST'] == 'ADD')]
            if add_row.empty:
                continue

            p_value = add_row['P'].values[0]
            odds_ratio = add_row['OR'].values[0]
            allele = add_row['A1'].values[0]
            l95 = add_row['L95'].values[0]
            u95 = add_row['U95'].values[0]

            scenario_summaries = {}
            ancestry_path = os.path.join(
                PROBS_FOLDER,
                'ancestry',
                hit,
                ancestry,
                f"{most_significant_snp}.tsv",
            )
            if not os.path.exists(ancestry_path):
                continue

            for scenario_key in SCENARIOS.keys():
                prob_path = os.path.join(
                    PROBS_FOLDER,
                    scenario_key,
                    hit,
                    ancestry,
                    f"{most_significant_snp}.tsv",
                )
                if not os.path.exists(prob_path):
                    scenario_summaries[scenario_key] = None
                    continue

                delta_df = pd.read_csv(prob_path, sep='\t')
                if delta_df.empty or 'delta_P' not in delta_df.columns:
                    scenario_summaries[scenario_key] = None
                    continue

                scenario_summaries[scenario_key] = {
                    'mean': delta_df['delta_P'].mean(),
                    'median': delta_df['delta_P'].median(),
                    'std': delta_df['delta_P'].std(),
                    'samples': delta_df['delta_P'].count(),
                }

            ancestry_summary = scenario_summaries.get('ancestry')
            if ancestry_summary is None:
                continue

            rows.append({
                'Phenotype': pheno_name,
                'ancestry tested': imp_ancestry,
                'ancestry of the population': ancestry,
                'Start Position of the reduced significant region': start_pos,
                'End Position of the reduced significant region': end_pos,
                'Number of significant SNPs': len(snps_list),
                'ID': most_significant_snp,
                'allele': allele,
                'OR (CI = 95%)': f'{odds_ratio} ({l95}, {u95})',
                'p value': p_value,
                'chr': chrom,
                'Delta_P_mean_ancestry': ancestry_summary['mean'],
                'Delta_P_median_ancestry': ancestry_summary['median'],
                'Delta_P_std_ancestry': ancestry_summary['std'],
                'Delta_P_mean_add': scenario_summaries.get('add', {}).get('mean') if scenario_summaries.get('add') else None,
                'Delta_P_median_add': scenario_summaries.get('add', {}).get('median') if scenario_summaries.get('add') else None,
                'Delta_P_std_add': scenario_summaries.get('add', {}).get('std') if scenario_summaries.get('add') else None,
                'Delta_P_mean_environment': scenario_summaries.get('environment', {}).get('mean') if scenario_summaries.get('environment') else None,
                'Delta_P_median_environment': scenario_summaries.get('environment', {}).get('median') if scenario_summaries.get('environment') else None,
                'Delta_P_std_environment': scenario_summaries.get('environment', {}).get('std') if scenario_summaries.get('environment') else None,
                'Delta_P_samples': ancestry_summary['samples'],
            })

    if not rows:
        return pd.DataFrame(columns=TABLE_COLUMNS)

    return pd.DataFrame(rows, columns=TABLE_COLUMNS)


def main():
    out_df = build_table()
    out_df.to_excel('table_fine_mapping.xlsx', index=False)


if __name__ == "__main__":
    main()
