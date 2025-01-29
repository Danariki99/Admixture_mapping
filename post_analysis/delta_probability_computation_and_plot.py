import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit

hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'
dataset_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/samples'
models_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/models'
plots_folder_genearl = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/plots'
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/results'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS', 'NAT']

hits_list = os.listdir(hit_folder_name)

for hit in hits_list:
    print(hit)
    pheno = hit.split('_')[1]
    imp_ancestry = hit.split('_')[0]
    deltas = []
    abs_deltas = []
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
            for snp in snps_list:
                dataset = pd.read_csv(f'{current_dataset_folder}/{snp}.tsv', sep='\t')
                dataset = dataset[dataset[imp_ancestry] == 1]

                model = pd.read_csv(f'{current_model_folder}/{snp}.csv')

                no_standardize_covars = ['IID', 'sex', imp_ancestry, 'ADD']

                columns_to_standardize = [col for col in dataset.columns if col not in no_standardize_covars]

                dataset[columns_to_standardize] = (dataset[columns_to_standardize] - dataset[columns_to_standardize].mean()) / dataset[columns_to_standardize].std()

                interesting_columns = model.columns.tolist()[3:]

                intercept = model['INTERCEPT']

                # Applica i coefficienti e l'intercetta a ogni riga del dataset
                dataset['z_with'] = intercept.values[0]  # Inizia con l'intercetta
                dataset['z_without'] = intercept.values[0]  # Inizia con l'intercetta

                for column in interesting_columns:
                    if column in dataset.columns:
                        if column == imp_ancestry:
                            dataset['z_without'] += dataset[column] * model[column].values[0]  # Escludi imp_ancestry
                        dataset['z_with'] += dataset[column] * model[column].values[0]  # Aggiungi il contributo

                # Calcola le probabilità e la differenza
                dataset['P_with'] = expit(dataset['z_with'])
                dataset['P_without'] = expit(dataset['z_without'])
                dataset['delta_P'] = dataset['P_with'] - dataset['P_without']
                dataset['delta_P_abs'] = abs(dataset['P_with'] - dataset['P_without'])

                df.to_csv(os.path.join(result_folder, f'{snp}.tsv'), sep='\t', index=False)

                if snp == most_significant_SNP:
                    deltas.append(dataset['delta_P'].tolist())
                    abs_deltas.append(dataset['delta_P_abs'].tolist())

        else:
                deltas.append([])
                abs_deltas.append([])
        
    # let's make the plot at the end of the iteration for all the ancestries
    plots_folder = os.path.join(plots_folder_genearl, hit)
    os.makedirs(plots_folder, exist_ok=True)

   # graphs
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.boxplot(deltas, tick_labels=ancestry_list, patch_artist=True)

    ax.set_xlabel("Ancestry")
    ax.set_ylabel("Delta Probabilities")
    ax.set_ylim(-1, 1)
    ax.set_title("Boxplot of Delta Probabilities by Ancestry")

    plt.savefig(os.path.join(plots_folder, f"delta_probabilities_boxplot.png"))

    # abs graph
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.boxplot(abs_deltas, tick_labels=ancestry_list, patch_artist=True)

    # Etichette e limiti
    ax.set_xlabel("Ancestry")
    ax.set_ylabel("Delta Probabilities")
    ax.set_ylim(0, 1)
    ax.set_title("Boxplot of |Delta Probabilities| by Ancestry")

    # Salva il grafico
    plt.savefig(os.path.join(plots_folder, f"abs_delta_probabilities_boxplot.png"))







            












