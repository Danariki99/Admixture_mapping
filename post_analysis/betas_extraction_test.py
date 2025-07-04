import pandas as pd
import os
import numpy as np
import sys

result_folder = sys.argv[1]

hit_folder_name = os.path.join(result_folder, 'fine_mapping_ancestries_PCA_verbose')

ancestry_list = ['AFR', 'EAS', 'EUR']

hits_list = os.listdir(hit_folder_name)

for hit in hits_list:
    print(f'Processing hit: {hit}')
    pheno = hit.split('_')[1]
    imp_ancestry = hit.split('_')[0]
    for ancestry in ancestry_list:
        output_folder = os.path.join(result_folder, f'probabities_pipeline/models/{hit}/{ancestry}')
        os.makedirs(output_folder, exist_ok=True)


        file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'
        if not os.path.exists(file):
            print(f'File {file} does not exist, skipping.')
            continue

        df = pd.read_csv(file, sep='\t')

        tests = df['TEST'].unique().tolist()

        ids = df['ID'].unique().tolist()
        significant_ids = []
        for id in ids:
            add_p_value = df[(df['ID'] == id) & (df['TEST'] == 'ADD')]['P']
            # Filtra per imp_ancestry
            imp_p_value = df[(df['ID'] == id) & (df['TEST'] == imp_ancestry)]['P']

            # Controlla se entrambi i test hanno P < 0.05
            if not add_p_value.empty and not imp_p_value.empty:
                print(add_p_value.values[0], imp_p_value.values[0])
                if add_p_value.values[0] < 0.7 and imp_p_value.values[0] < 0.1:
                    print(f'Adding significant ID: {id} with ADD P-value: {add_p_value.values[0]} and {imp_ancestry} P-value: {imp_p_value.values[0]}')
                    significant_ids.append(id)

        #now that we have all the informations needed, we can create and save the models
        # we can create a list of dataframes, one for each significant hit and windows to have the correct betas.

        for id in significant_ids:
            # Filtra per ID specifico
            df_id = df[df['ID'] == id]

            # Inizializza un DataFrame vuoto per i betas
            df_beta = pd.DataFrame({'ID': [id]})  # Usa una lista per evitare problemi di broadcasting

            # Aggiungi la colonna CHROM
            chrom_values = df_id[df_id['TEST'] == 'ADD']['#CHROM'].values
            if len(chrom_values) > 0:
                df_beta['CHROM'] = chrom_values[0]  # Prendi il primo valore, se esiste

            for test in tests:
                test_row = df_id[df_id['TEST'] == test]
                
                if not test_row.empty:
                    p_value = test_row['P'].values[0]  # Estrai il valore P
                    or_value = test_row['OR'].values[0]  # Estrai l'OR

                    if p_value < 0.05:
                        try:
                            df_beta[test] = np.log(or_value)  # Prova a calcolare il logaritmo
                        except (ValueError, TypeError) as e:
                            print(f"Errore nel calcolare il logaritmo per {test}: {e}. Il file sarà saltato.")
                            continue  # Salta questa iterazione e non salva il file
                    elif test == 'ADD' or test == imp_ancestry:
                        try:
                            df_beta[test] = np.log(or_value)  # Prova a calcolare il logaritmo
                        except (ValueError, TypeError) as e:
                            print(f"Errore nel calcolare il logaritmo per {test}: {e}. Il file sarà saltato.")
                            continue  # Salta questa iterazione e non salva il file

                # Se non ci sono errori, salva il file
                if 'ADD' in df_beta.columns and imp_ancestry in df_beta.columns:
                    if df_beta.select_dtypes(include=[np.number]).applymap(np.isfinite).all().all():
                        df_beta.to_csv(f'{output_folder}/{id}.csv', index=False)
                    else:
                        print(f"⚠️  Skip saving {id} due to non-finite values in numeric columns.")
                else:
                    print(f"⚠️  Skip saving {id} because required columns are missing: ADD and/or {imp_ancestry}")

