import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit
from matplotlib.patches import Patch
import numpy as np

# Colori e ancestry
ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS', 'NAT']
colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown']
color_map = {ancestry: colors[i % len(colors)] for i, ancestry in enumerate(ancestry_list)}

# Cartelle aggiornate
hit_folder_name = './results/fine_mapping_ancestries_PCA_verbose'
dataset_folder = './results/fine_mapping_samples'
models_folder = './results/fine_mapping_models'
plots_folder_general = './results/plots'
output_folder = './results/probabilities'
probs_folder = './results/probabilities_probs'
pheno_file = './data/phenotypes/ukbb_v1.xlsx'

os.makedirs(output_folder, exist_ok=True)
os.makedirs(plots_folder_general, exist_ok=True)
os.makedirs(probs_folder, exist_ok=True)

# Carica fenotipi
df_first_batch = pd.read_excel(pheno_file, sheet_name="first_batch")

# Contenitore boxplot
boxplot_data = []

def data_processing():
    for hit in os.listdir(hit_folder_name):
        print(f"Processing {hit}")
        hit_parts = hit.split('_')
        imp_ancestry = hit_parts[0]
        pheno = hit_parts[1]
        chrom = hit_parts[-1]

        # Nome leggibile
        pheno_name = df_first_batch[df_first_batch['ID'] == pheno]['ID2']
        pheno_name = '_'.join(pheno_name.iloc[0].split('_')[1:]) if not pheno_name.empty else 'Unknown'
        if pheno_name == 'TTE_acute_upper_respiratory_infections_of_multiple_and_unspecified_sites':
            pheno_name = 'TTE_acute_upper_respiratory_infections'
        pheno_name = pheno_name.replace('_', ' ')
        label = f"{pheno_name} ({imp_ancestry}, {chrom})"

        for ancestry in ancestry_list:
            glm_file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'
            if not os.path.isfile(glm_file):
                continue

            df = pd.read_csv(glm_file, sep='\t')
            current_model_folder = f'{models_folder}/{hit}/{ancestry}'
            current_dataset_folder = f'{dataset_folder}/{hit}/{ancestry}'

            result_folder = os.path.join(output_folder, hit, ancestry)
            os.makedirs(result_folder, exist_ok=True)

            output_probs_folder = os.path.join(probs_folder, hit, ancestry)
            os.makedirs(output_probs_folder, exist_ok=True)

            snps_list = [f.replace('.csv', '') for f in os.listdir(current_model_folder) if f.endswith('.csv')]
            if not snps_list:
                continue

            filtered_df = df[(df['TEST'] == imp_ancestry) & (df['ID'].isin(snps_list))]
            if filtered_df.empty:
                continue

            most_significant_SNP = filtered_df.loc[filtered_df['P'].idxmin()]['ID']
            for snp in snps_list:
                dataset_file = f'{current_dataset_folder}/{snp}.tsv'
                model_file = f'{current_model_folder}/{snp}.csv'
                if not os.path.isfile(dataset_file) or not os.path.isfile(model_file):
                    continue

                dataset = pd.read_csv(dataset_file, sep='\t')
                dataset = dataset[dataset[imp_ancestry] == 1]
                model = pd.read_csv(model_file)

                if 'INTERCEPT' not in model.columns:
                    continue

                no_standardize_covars = ['IID', 'sex', imp_ancestry, 'ADD']
                columns_to_standardize = [col for col in dataset.columns if col not in no_standardize_covars]
                dataset[columns_to_standardize] = (dataset[columns_to_standardize] - dataset[columns_to_standardize].mean()) / dataset[columns_to_standardize].std()

                intercept = model['INTERCEPT'].values[0]
                dataset['z_with'] = intercept
                dataset['z_without'] = intercept
                for column in model.columns[3:]:
                    if column in dataset.columns:
                        if column == imp_ancestry:
                            dataset['z_without'] += dataset[column] * model[column].values[0]
                        dataset['z_with'] += dataset[column] * model[column].values[0]

                dataset['P_with'] = expit(dataset['z_with'])
                dataset['P_without'] = expit(dataset['z_without'])
                dataset['delta_P'] = dataset['P_with'] - dataset['P_without']
                dataset['delta_P_abs'] = abs(dataset['delta_P'])

                dataset.to_csv(os.path.join(result_folder, f'{snp}.tsv'), sep='\t', index=False)

                if snp == most_significant_SNP:
                    boxplot_data.append({
                        'label': label,
                        'ancestry': 'AMR' if ancestry == 'NAT' else ancestry,
                        'snp': snp,
                        'delta': dataset['delta_P'].tolist(),
                        'delta_abs': dataset['delta_P_abs'].tolist(),
                        'color': color_map[ancestry]
                    })
                    dataset.to_csv(os.path.join(output_probs_folder, f'{snp}.tsv'), sep='\t', index=False)

def plot_filtered_boxplot(data_key, ylabel, title, filename):
    filtered = [b for b in boxplot_data if len(b[data_key]) > 0]
    if not filtered:
        print("Nessun dato da plottare.")
        return

    filtered.sort(key=lambda x: np.median(x[data_key]))
    spacing = 2
    box_positions = [i * spacing for i in range(1, len(filtered) + 1)]

    fig, ax = plt.subplots(figsize=(max(14, len(filtered) * 0.75), 10))
    box_data = [b[data_key] for b in filtered]
    bp = ax.boxplot(box_data, positions=box_positions, patch_artist=True)

    for patch, b in zip(bp['boxes'], filtered):
        ancestry = b['ancestry']
        patch.set_facecolor(color_map.get(ancestry, 'gray'))

    x_labels = [b['label'] for b in filtered]
    ax.set_xticks(box_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.set_ylim(-1 if ylabel == "Delta Probabilities" else 0, 1.1)

    max_values = [max(b[data_key]) for b in filtered]
    for i, (pos, b, max_val) in enumerate(zip(box_positions, filtered, max_values)):
        y_pos = min(max_val + (0.25 if i % 2 == 0 else 0.05), 1.08)
        ax.text(pos, y_pos, b['snp'], ha='center', va='bottom', fontsize=11)

    legend_elements = [Patch(facecolor=color_map.get(a, 'gray'), label=a) for a in sorted(set(b['ancestry'] for b in filtered))]
    ax.legend(handles=legend_elements, title="Global Ancestry", loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=11, title_fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(plots_folder_general, filename), bbox_inches='tight')
    plt.close()

# Esegui
if __name__ == "__main__":
    data_processing()
    plot_filtered_boxplot('delta', "Delta Probabilities", "Boxplot of Delta Probabilities", "delta_probabilities.pdf")
    plot_filtered_boxplot('delta_abs', "|Delta Probabilities|", "Boxplot of |Delta Probabilities|", "abs_delta_probabilities.pdf")
