import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit
from matplotlib.patches import Patch
import random
import numpy as np
import sys

def data_processing(result_folder):

    # Cartelle
    hit_folder_name = os.path.join(result_folder, 'fine_mapping_ancestries_PCA_verbose')
    dataset_folder = os.path.join(result_folder, 'probabities_pipeline/samples')
    models_folder = os.path.join(result_folder, 'probabities_pipeline/models')
    plots_folder_general = os.path.join(result_folder, 'plots')
    output_folder = os.path.join(result_folder, 'probabities_pipeline/results')
    probs_folder = os.path.join(result_folder, 'probabities_pipeline/probs')

    ancestry_list = ['AFR', 'EAS', 'EUR']

    colors = ['blue', 'green', 'red']
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

        # Etichetta fenotipo arricchita con ancestry e chrom
        label = f"{pheno} ({imp_ancestry}, {chrom})"

        for ancestry in ancestry_list:
            file = f'{hit_folder_name}/{hit}/{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid'
            if not os.path.exists(file):
                print(f"File not found: {file}")
                continue
            df = pd.read_csv(file, sep='\t')

            current_model_folder = f'{models_folder}/{hit}/{ancestry}'
            current_dataset_folder = f'{dataset_folder}/{hit}/{ancestry}'

            results_folder = os.path.join(output_folder, hit, ancestry)
            os.makedirs(results_folder, exist_ok=True)

            output_probs_folder = os.path.join(probs_folder, hit, ancestry)
            os.makedirs(output_probs_folder, exist_ok=True)

            snps_list = [f.replace('.csv', '') for f in os.listdir(current_model_folder) if f.endswith('.csv')]

            if len(snps_list) == 0:
                continue

            filtered_df = df[(df['TEST'] == imp_ancestry) & (df['ID'].isin(snps_list))]
            if filtered_df.empty:
                continue
            
            ordered_snps = filtered_df.sort_values(by='P')['ID'].tolist()

            # Se la lista è vuota, passa oltre
            if not ordered_snps:
                print(f"No SNPs found for {hit} with ancestry {ancestry}.")
                continue

            valid_snp_found = False
            i = 0

            while not valid_snp_found and i < len(ordered_snps):
                snp = ordered_snps[i]
                i += 1

                dataset_path = os.path.join(current_dataset_folder, f'{snp}.tsv')
                model_path = os.path.join(current_model_folder, f'{snp}.csv')

                if not os.path.exists(dataset_path) or not os.path.exists(model_path):
                    continue

                dataset = pd.read_csv(dataset_path, sep='\t')
                dataset = dataset[dataset[imp_ancestry] == 1]
                model = pd.read_csv(model_path)

                if 'INTERCEPT' not in model.columns:
                    continue

                no_standardize_covars = ['IID', 'sex', imp_ancestry, 'ADD']
                columns_to_standardize = [col for col in dataset.columns if col not in no_standardize_covars]
                dataset[columns_to_standardize] = (
                    (dataset[columns_to_standardize] - dataset[columns_to_standardize].mean())
                    / dataset[columns_to_standardize].std()
                )

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

                if not dataset.isnull().values.any() and np.isfinite(dataset.select_dtypes(include=[np.number])).all().all():
                    dataset.to_csv(os.path.join(results_folder, f'{snp}.tsv'), sep='\t', index=False)
                    boxplot_data.append({
                        'label': label,
                        'ancestry': 'AMR' if ancestry == 'NAT' else ancestry,
                        'snp': snp,
                        'delta': dataset['delta_P'].tolist(),
                        'delta_abs': dataset['delta_P_abs'].tolist(),
                        'color': color_map[ancestry]
                    })
                    dataset.to_csv(os.path.join(output_probs_folder, f'{snp}.tsv'), sep='\t', index=False)
                    print(f"✅ Added SNP: {snp} for boxplot.")
                    valid_snp_found = True
                else:
                    print(f"⚠️  Skip {snp} due to NaN/Inf values.")
                
    return boxplot_data


                

# Funzione per plottare i risultati filtrati
def plot_filtered_boxplot(boxplot_data, data_key, ylabel, title, filename):
    filtered = [b for b in boxplot_data if len(b[data_key]) > 0]
    if not filtered:
        print("Nessun dato da plottare.")
        return

    # Ordina in base alla mediana dei dati
    filtered.sort(key=lambda x: np.median(x[data_key]))

    spacing = 2  # maggiore separazione
    box_positions = [i * spacing for i in range(1, len(filtered) + 1)]

    fig, ax = plt.subplots(figsize=(max(14, len(filtered) * 0.75), 10))

    # Lista di liste dei dati per il boxplot
    box_data = [b[data_key] for b in filtered]

    # ✅ Calcola e stampa la media per ciascun SNP
    for b in filtered:
        snp_label = b['snp']
        mean_value = np.mean(b[data_key])
        print(f"Mean for SNP {snp_label}: {mean_value:.4f}")

    # ⏹️ Costruisci il boxplot come prima
    bp = ax.boxplot(box_data, positions=box_positions, patch_artist=True)


    for patch, b in zip(bp['boxes'], filtered):
        ancestry = b['ancestry']
        patch.set_facecolor(color_map.get('NAT') if ancestry == 'AMR' else color_map.get(ancestry, 'gray'))

    x_labels = [b['label'] for b in filtered]
    ax.set_xticks(box_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=12)

    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.set_ylim(-1 if ylabel == "Delta Probabilities" else 0, 1.1)

    # Mostra SNP sopra il boxplot alternando l'altezza
    max_values = [max(b[data_key]) for b in filtered]
    for i, (pos, b, max_val) in enumerate(zip(box_positions, filtered, max_values)):
        if i % 2 == 0:
            y_pos = min(max_val + 0.25, 1.08)
        else:
            y_pos = min(max_val + 0.05, 1.08)
        ax.text(pos, y_pos, b['snp'], ha='center', va='bottom', fontsize=11)

    # 🔥 Legenda corretta: costruisci ancestry -> colore
    ancestry_set = sorted(set(b['ancestry'] for b in filtered))
    legend_elements = []

    for ancestry in ancestry_set:
        color = color_map.get('NAT') if ancestry == 'AMR' else color_map.get(ancestry, 'gray')
        legend_elements.append(Patch(facecolor=color, label=ancestry))

    ax.legend(handles=legend_elements, title="Global Ancestry", loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=11, title_fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(plots_folder_general, filename), bbox_inches='tight')
    plt.close()




# Plotta risultati
if __name__ == "__main__":
    result_folder = sys.argv[1]

    boxplot_data = data_processing(result_folder)
    # Cartelle
    probs_folder = os.path.join(result_folder, 'probabities_pipeline/probs')
    plots_folder_general = os.path.join(result_folder, 'plots')

    # Colori per ancestry
    ancestry_list = ['AFR', 'EAS', 'EUR']
    colors = ['blue', 'green', 'red']
    color_map = {ancestry: colors[i % len(colors)] for i, ancestry in enumerate(ancestry_list)}
    # Plotta i grafici
    plot_filtered_boxplot(boxplot_data, 'delta', "Delta Probabilities", "Boxplot of Delta Probabilities by Ancestry for UKBB", "delta_probabilities_UKBB.pdf")
    plot_filtered_boxplot(boxplot_data, 'delta_abs', "|Delta Probabilities|", "Boxplot of |Delta Probabilities| by Ancestry for UKBB", "abs_delta_probabilities_UKBB.pdf")

