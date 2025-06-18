import os
import pandas as pd
import numpy as np
import sys

# Cartella input/output secondo pipeline
result_folder = sys.argv[1]

input_root = os.path.join(result_folder, 'fine_mapping_verbose')
output_root = os.path.join(result_folder, 'fine_mapping_models')

ancestry_list = ['AFR', 'EAS', 'EUR']


for pheno in os.listdir(input_root):
    pheno_path = os.path.join(input_root, pheno)
    if not os.path.isdir(pheno_path):
        continue

    for chr_folder in os.listdir(pheno_path):
        chr_path = os.path.join(pheno_path, chr_folder)

        for ancestry in ancestry_list:
            glm_file = os.path.join(chr_path, ancestry, 'output.glm.logistic.hybrid')

            if not os.path.exists(glm_file):
                continue

            print(f"Processing: {pheno}, {chr_folder}, {ancestry}")
            df = pd.read_csv(glm_file, sep='\t')

            # Determina i test presenti
            tests = df['TEST'].unique().tolist()
            ids = df['ID'].unique().tolist()

            # Identifica SNP significativi per ADD + ancestry
            significant_ids = []
            for id in ids:
                add_p = df[(df['ID'] == id) & (df['TEST'] == 'ADD')]['P']
                anc_p = df[(df['ID'] == id) & (df['TEST'] == ancestry)]['P']

                if not add_p.empty and not anc_p.empty:
                    if add_p.values[0] < 0.05 and anc_p.values[0] < 0.05:
                        significant_ids.append(id)

            # Salva modelli
            for id in significant_ids:
                df_id = df[df['ID'] == id]
                df_beta = pd.DataFrame({'ID': [id]})

                chrom_values = df_id[df_id['TEST'] == 'ADD']['#CHROM'].values
                if len(chrom_values) > 0:
                    df_beta['CHROM'] = chrom_values[0]

                for test in tests:
                    test_row = df_id[df_id['TEST'] == test]
                    if not test_row.empty:
                        p = test_row['P'].values[0]
                        or_ = test_row['OR'].values[0]
                        if p < 0.05:
                            df_beta[test] = np.log(or_)

                # Salva file modello
                model_output_folder = os.path.join(output_root, pheno, chr_folder, ancestry)
                os.makedirs(model_output_folder, exist_ok=True)
                df_beta.to_csv(os.path.join(model_output_folder, f'{id}.csv'), index=False)

