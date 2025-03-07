import os
import pandas as pd

out1 = pd.DataFrame(columns=[
    'Phenotype',
    'ID',
    'ancestry tested',
    'ancestry of the population',
    'OR (CI = 95%)',
    'p value',
    'chr',
    'Delta_P_mean',
    'Delta_P_std',
    'Delta_P_median',
    'Delta_P_samples',
])

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'
dataset_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/samples'
models_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/models'
plots_folder_genearl = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/plots'
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/results'
probs_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/probs'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS', 'NAT']

df_first_batch = pd.read_excel("ukbb_v1.xlsx", sheet_name="first_batch")

hits_list = os.listdir(hit_folder_name)

for hit in hits_list:
    print(hit)
    pheno = hit.split('_')[1]
    imp_ancestry = hit.split('_')[0]
    most_significant_SNPs = []
    chrom = hit.split('/')[-1].split('_')[2]
    pheno_name = df_first_batch[df_first_batch['ID'] == pheno]['ID2']
    for ancestry in ancestry_list:

        file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'
        df = pd.read_csv(file, sep='\t')

        current_model_folder = f'{models_folder}/{hit}/{ancestry}'
        current_dataset_folder = f'{dataset_folder}/{hit}/{ancestry}'

        result_folder = os.path.join(output_folder, hit, ancestry)
        os.makedirs(result_folder, exist_ok=True)

        snps_list = os.listdir(current_model_folder)
        snps_list = [file_name.replace('.csv', '') for file_name in snps_list if file_name.endswith('.csv')]
        if len(snps_list) != 0:

            filtered_df = df[(df['TEST'] == imp_ancestry) & (df['ID'].isin(snps_list))]
            most_significant_SNP = filtered_df.loc[filtered_df['P'].idxmin()]['ID']
            most_significant_SNPs.append(most_significant_SNP)
            row_with_min_p = df[(df['ID'] == most_significant_SNP) & (df['TEST'] == 'ADD')]
            p = row_with_min_p['P'].values[0]
            oddr = row_with_min_p['OR'].values[0]
            snpid = row_with_min_p['ID'].values[0]

            delta_p_df = pd.read_csv(f'{probs_folder}/{hit}/{ancestry}/{most_significant_SNP}.tsv', sep='\t')

            out1 = pd.concat([
                out1,
                pd.DataFrame([{
                    'Phenotype': pheno_name.iloc[0] if not pheno_name.empty else 'Unknown',  # Nome del fenotipo
                    'ID': snpid,
                    'ancestry tested': imp_ancestry,  
                    'ancestry of the population': ancestry,
                    'OR (CI = 95%)': f'{oddr} ({row_with_min_p["L95"].values[0]}, {row_with_min_p["U95"].values[0]})', 
                    'p value': p,
                    'chr': chrom,
                    'Delta_P_mean': delta_p_df['delta_P'].mean(),
                    'Delta_P_std': delta_p_df['delta_P'].std(),
                    'Delta_P_median': delta_p_df['delta_P'].median(),
                    'Delta_P_samples': len(delta_p_df['delta_P'])
                }])

            ], ignore_index=True)

out1.to_excel('table_fine_mapping.xlsx', index=False)
