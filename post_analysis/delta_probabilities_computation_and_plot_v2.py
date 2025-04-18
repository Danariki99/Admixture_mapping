import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit
from matplotlib.patches import Patch
import random

# Carica la tabella dei fenotipi
df_first_batch = pd.read_excel("../tables_plots/ukbb_v1.xlsx", sheet_name="first_batch")

# Cartelle
hit_folder_name = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/fine_mapping_ancestries_PCA_verbose'
dataset_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/samples'
models_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/models'
plots_folder_general = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/plots'
output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/results'
probs_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/probabilities_pipeline/probs'

ancestry_list = ['AFR', 'EAS', 'EUR', 'SAS', 'WAS', 'NAT']
colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown']
color_map = {ancestry: colors[i % len(colors)] for i, ancestry in enumerate(ancestry_list)}

# Contenitore unificato
boxplot_data = []

hits_list = os.listdir(hit_folder_name)

for hit in hits_list:
    print(hit)
    hit_parts = hit.split('_')
    imp_ancestry = hit_parts[0]
    pheno = hit_parts[1]
    chrom = hit_parts[-1]

    # Ottieni nome fenotipo leggibile
    pheno_name = df_first_batch[df_first_batch['ID'] == pheno]['ID2']
    pheno_name = '_'.join(pheno_name.iloc[0].split('_')[1:]) if not pheno_name.empty else 'Unknown'
    if pheno_name == 'TTE_acute_upper_respiratory_infections_of_multiple_and_unspecified_sites':
        pheno_name = 'TTE_acute_upper_respiratory_infections'
    pheno_name = pheno_name.replace('_', ' ')

    # Etichetta fenotipo arricchita con ancestry e chrom
    label = f"{pheno_name} ({imp_ancestry}, {chrom})"

    for ancestry in ancestry_list:
        file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'
        df = pd.read_csv(file, sep='\t')

        current_model_folder = f'{models_folder}/{hit}/{ancestry}'
        current_dataset_folder = f'{dataset_folder}/{hit}/{ancestry}'

        result_folder = os.path.join(output_folder, hit, ancestry)
        os.makedirs(result_folder, exist_ok=True)

        output_probs_folder = os.path.join(probs_folder, hit, ancestry)
        os.makedirs(output_probs_folder, exist_ok=True)

        snps_list = [f.replace('.csv', '') for f in os.listdir(current_model_folder) if f.endswith('.csv')]

        if len(snps_list) == 0:
            continue

        filtered_df = df[(df['TEST'] == imp_ancestry) & (df['ID'].isin(snps_list))]
        if filtered_df.empty:
            continue

        most_significant_SNP = filtered_df.loc[filtered_df['P'].idxmin()]['ID']
        for snp in snps_list:
            dataset = pd.read_csv(f'{current_dataset_folder}/{snp}.tsv', sep='\t')
            dataset = dataset[dataset[imp_ancestry] == 1]
            model = pd.read_csv(f'{current_model_folder}/{snp}.csv')

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

# Funzione per plottare i risultati filtrati
def plot_filtered_boxplot(data_key, ylabel, title, filename):
    filtered = [b for b in boxplot_data if len(b[data_key]) > 0]
    if not filtered:
        print("Nessun dato da plottare.")
        return

    fig, ax = plt.subplots(figsize=(max(14, len(filtered)*0.6), 10))
    box_data = [b[data_key] for b in filtered]
    box_positions = list(range(1, len(filtered) + 1))
    bp = ax.boxplot(box_data, positions=box_positions, patch_artist=True)

    for patch, b in zip(bp['boxes'], filtered):
        patch.set_facecolor(b['color'])

    x_labels = [b['label'] for b in filtered]
    ax.set_xticks(box_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=12)

    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.set_ylim(-1 if ylabel == "Delta Probabilities" else 0, 1.1)

    # Mostra SNP sopra il boxplot
    max_values = [max(b[data_key]) for b in filtered]
    for pos, b, max_val in zip(box_positions, filtered, max_values):
        y_pos = min(max_val + random.uniform(0.05, 0.15), 1.08)
        ax.text(pos, y_pos, b['snp'], ha='center', va='bottom', fontsize=11, rotation=0)

    # Legenda ancestry
    legend_labels = {b['ancestry']: b['color'] for b in filtered}
    legend_elements = [Patch(facecolor=color, label=label) for label, color in legend_labels.items()]
    ax.legend(handles=legend_elements, title="Ancestry", loc='upper right', fontsize=11, title_fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(plots_folder_general, filename))
    plt.close()

# Plotta risultati
plot_filtered_boxplot('delta', "Delta Probabilities", "Boxplot of Delta Probabilities by Ancestry for UKBB", "delta_probabilities_UKBB.png")
plot_filtered_boxplot('delta_abs', "|Delta Probabilities|", "Boxplot of |Delta Probabilities| by Ancestry for UKBB", "abs_delta_probabilities_UKBB.png")
