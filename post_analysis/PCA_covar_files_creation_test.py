import os
import pandas as pd
import sys

ancestries = ['AFR', 'EAS', 'EUR']
# === Path settings ===

results_folder = sys.argv[1]


covar_path = os.path.join(results_folder, 'wind_covar_files_processed')
eigenvec_base_path = os.path.join(results_folder, 'PCA_files', 'PCA_res', 'PCA_*.out.eigenvec')
output_folder = os.path.join(results_folder, 'PCA_files', 'PCA_covar_files')

os.makedirs(output_folder, exist_ok=True)

# === Loop over all covariate files ===
for covar_file in os.listdir(covar_path):

    covar_df = pd.read_csv(os.path.join(covar_path, covar_file), sep='\t')

    for ancestry in ancestries:
        eigenvec_path = eigenvec_base_path.replace('*', ancestry)
        if os.path.exists(eigenvec_path):
            eigenvec_df = pd.read_csv(eigenvec_path, sep='\t')
        else:
            print(f"⚠️ Eigenvec file for ancestry {ancestry} not found: {eigenvec_path}")
            continue
        eigenvec_df.rename(columns={'#IID': 'IID'}, inplace=True)

        # Rimuovi colonne PC globali dal file covar base
        covar_df_filtered = covar_df.drop(columns=[col for col in covar_df.columns if col.startswith("Global_PC")])

        covar_basic = covar_df_filtered[['IID', 'age', 'sex', 'BMI']]
        covar_rest = covar_df_filtered.drop(columns=['age', 'sex', 'BMI'])

        # Merge triplo (basic + eigenvec + rest)
        merged_df = pd.merge(pd.merge(covar_basic, eigenvec_df, on='IID', how='inner'), covar_rest, on='IID')

        output_path = os.path.join(output_folder, f"{covar_file.split('.')[0]}_{ancestry}.tsv")
        merged_df.to_csv(output_path, sep='\t', index=False)

print("✅ PCA covariate files created for all ancestries.")
