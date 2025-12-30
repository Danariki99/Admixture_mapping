import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit
from matplotlib.patches import Patch
import numpy as np
import sys

# Filtri dichiarati nel paper
MIN_SAMPLE_THRESHOLD = 0
MIN_DELTA_THRESHOLD = 0.0

ANCESTRY_LIST = ['AFR', 'EAS', 'EUR']
COLORS = ['blue', 'green', 'red']
COLOR_MAP = {ancestry: COLORS[i % len(COLORS)] for i, ancestry in enumerate(ANCESTRY_LIST)}


def exclude_ancestry(imp_ancestry, env_columns):
    return {imp_ancestry}


def exclude_add(imp_ancestry, env_columns):
    return {'ADD'}


def exclude_environment(imp_ancestry, env_columns):
    return set(env_columns)


SCENARIOS = {
    'ancestry': {
        'exclude_fn': exclude_ancestry,
        'plot_label': 'Local Ancestry',
        'file_suffix': 'ancestry',
    },
    'add': {
        'exclude_fn': exclude_add,
        'plot_label': 'Genotype (ADD)',
        'file_suffix': 'add',
    },
    'environment': {
        'exclude_fn': exclude_environment,
        'plot_label': 'Environmental Covariates',
        'file_suffix': 'environment',
    },
}


def standardize_covariates(dataset, imp_ancestry):
    protected_cols = {'IID', imp_ancestry, 'ADD'}
    if 'sex' in dataset.columns:
        protected_cols.add('sex')

    candidate_cols = [
        col for col in dataset.columns
        if col not in protected_cols and pd.api.types.is_numeric_dtype(dataset[col])
    ]

    if not candidate_cols:
        return candidate_cols

    means = dataset[candidate_cols].mean()
    stds = dataset[candidate_cols].std(ddof=0)
    stds.replace(to_replace=0, value=1, inplace=True)
    dataset[candidate_cols] = (dataset[candidate_cols] - means) / stds
    return candidate_cols


def has_valid_numeric_values(dataframe):
    numeric_df = dataframe.select_dtypes(include=[np.number])
    if numeric_df.empty:
        return True
    return np.isfinite(numeric_df.values).all() and not numeric_df.isnull().values.any()


def build_env_columns(dataset, imp_ancestry):
    base_exclusions = {'IID', imp_ancestry, 'ADD'}
    return [col for col in dataset.columns if col not in base_exclusions]


def compute_linear_predictor(dataset, intercept, coefficients, allowed_columns):
    z = np.full(len(dataset), intercept, dtype=float)
    for column, beta in coefficients.items():
        if column in allowed_columns:
            z += dataset[column].astype(float).values * beta
    return z


def list_cached_snps(result_dirs):
    cached = {scenario: set() for scenario in result_dirs}
    for scenario, dir_path in result_dirs.items():
        if not os.path.isdir(dir_path):
            continue
        for file_name in os.listdir(dir_path):
            if file_name.endswith('.tsv'):
                cached[scenario].add(file_name.replace('.tsv', ''))
    return cached


def append_cached_boxplot_data(prob_dirs, boxplot_data_by_scenario, snp, ancestry, label):
    for scenario_name, prob_dir in prob_dirs.items():
        prob_path = os.path.join(prob_dir, f'{snp}.tsv')
        if not os.path.exists(prob_path):
            continue

        df = pd.read_csv(prob_path, sep='\t')
        if 'delta_P' not in df.columns or 'delta_P_abs' not in df.columns:
            continue

        boxplot_data_by_scenario[scenario_name].append({
            'label': label,
            'ancestry': 'AMR' if ancestry == 'NAT' else ancestry,
            'snp': snp,
            'delta': df['delta_P'].tolist(),
            'delta_abs': df['delta_P_abs'].tolist(),
            'color': COLOR_MAP.get(ancestry, 'gray'),
        })


def data_processing(result_folder):
    hit_folder_name = os.path.join(result_folder, 'fine_mapping_ancestries_PCA_verbose')
    dataset_folder = os.path.join(result_folder, 'probabities_pipeline/samples')
    models_folder = os.path.join(result_folder, 'probabities_pipeline/models')
    output_folder = os.path.join(result_folder, 'probabities_pipeline/results')
    probs_folder = os.path.join(result_folder, 'probabities_pipeline/probs')

    boxplot_data_by_scenario = {name: [] for name in SCENARIOS}

    if not os.path.isdir(hit_folder_name):
        print(f"Hit folder not found: {hit_folder_name}")
        return boxplot_data_by_scenario

    hits_list = os.listdir(hit_folder_name)

    for hit in hits_list:
        print(hit)
        hit_parts = hit.split('_')
        if len(hit_parts) < 3:
            continue

        imp_ancestry = hit_parts[0]
        pheno = hit_parts[1]
        chrom = hit_parts[-1]

        label = f"{pheno} ({imp_ancestry}, {chrom})"

        for ancestry in ANCESTRY_LIST:
            glm_path = os.path.join(
                hit_folder_name,
                hit,
                f"{hit}_output.{ancestry}.{pheno}.glm.logistic.hybrid",
            )
            if not os.path.exists(glm_path):
                continue
            df = pd.read_csv(glm_path, sep='\t')

            current_model_folder = os.path.join(models_folder, hit, ancestry)
            current_dataset_folder = os.path.join(dataset_folder, hit, ancestry)

            if not os.path.isdir(current_model_folder) or not os.path.isdir(current_dataset_folder):
                continue

            scenario_result_dirs = {}
            scenario_prob_dirs = {}
            for scenario_name in SCENARIOS:
                result_dir = os.path.join(output_folder, scenario_name, hit, ancestry)
                prob_dir = os.path.join(probs_folder, scenario_name, hit, ancestry)
                os.makedirs(result_dir, exist_ok=True)
                os.makedirs(prob_dir, exist_ok=True)
                scenario_result_dirs[scenario_name] = result_dir
                scenario_prob_dirs[scenario_name] = prob_dir

            cached_snps = list_cached_snps(scenario_result_dirs)

            snps_list = [
                f.replace('.csv', '')
                for f in os.listdir(current_model_folder)
                if f.endswith('.csv')
            ]
            if not snps_list:
                continue

            filtered_df = df[(df['TEST'] == imp_ancestry) & (df['ID'].isin(snps_list))]
            filtered_df = filtered_df.dropna(subset=['P'])
            if filtered_df.empty:
                continue

            most_significant_SNP = filtered_df.loc[filtered_df['P'].idxmin()]['ID']

            for snp in snps_list:
                dataset_path = os.path.join(current_dataset_folder, f'{snp}.tsv')
                model_path = os.path.join(current_model_folder, f'{snp}.csv')

                if not os.path.exists(dataset_path) or not os.path.exists(model_path):
                    continue

                already_done = all(
                    snp in cached_snps.get(scenario_name, set())
                    for scenario_name in SCENARIOS
                )
                if already_done:
                    if snp == most_significant_SNP:
                        append_cached_boxplot_data(
                            scenario_prob_dirs,
                            boxplot_data_by_scenario,
                            snp,
                            ancestry,
                            label,
                        )
                    continue

                dataset = pd.read_csv(dataset_path, sep='\t')
                dataset = dataset[dataset[imp_ancestry] == 1]
                if len(dataset) < MIN_SAMPLE_THRESHOLD:
                    continue

                model = pd.read_csv(model_path)

                if 'INTERCEPT' not in model.columns:
                    continue

                standardize_covariates(dataset, imp_ancestry)

                env_columns = build_env_columns(dataset, imp_ancestry)

                intercept = float(model['INTERCEPT'].values[0])
                coefficient_columns = [
                    col for col in model.columns
                    if col not in {'ID', 'CHROM', 'INTERCEPT'}
                    and col in dataset.columns
                    and not pd.isnull(model[col].values[0])
                ]
                coefficients = {
                    col: float(model[col].values[0])
                    for col in coefficient_columns
                }

                if not coefficients:
                    continue

                z_full = compute_linear_predictor(dataset, intercept, coefficients, coefficients.keys())
                p_full = expit(z_full)

                if not np.isfinite(p_full).all():
                    continue

                scenario_outputs = {}
                ancestry_passed = False

                for scenario_name, scenario_conf in SCENARIOS.items():
                    exclude_candidates = scenario_conf['exclude_fn'](imp_ancestry, env_columns)
                    exclude_columns = {col for col in exclude_candidates if col in coefficients}
                    allowed_columns = set(coefficients.keys()) - exclude_columns

                    z_without = compute_linear_predictor(dataset, intercept, coefficients, allowed_columns)
                    p_without = expit(z_without)
                    delta = p_full - p_without
                    delta_abs = np.abs(delta)

                    max_delta_abs = np.max(delta_abs) if len(delta_abs) else 0.0

                    if scenario_name == 'ancestry':
                        if max_delta_abs < MIN_DELTA_THRESHOLD:
                            ancestry_passed = False
                            scenario_outputs = {}
                            break
                        ancestry_passed = True

                    scenario_df = dataset.copy()
                    scenario_df['P_with'] = p_full
                    scenario_df['P_without'] = p_without
                    scenario_df['delta_P'] = delta
                    scenario_df['delta_P_abs'] = delta_abs
                    scenario_df['scenario'] = scenario_name

                    if not has_valid_numeric_values(scenario_df):
                        continue

                    scenario_outputs[scenario_name] = {
                        'dataframe': scenario_df,
                        'delta': delta,
                        'delta_abs': delta_abs,
                    }

                if not ancestry_passed or 'ancestry' not in scenario_outputs:
                    continue

                for scenario_name, scenario_data in scenario_outputs.items():
                    scenario_df = scenario_data['dataframe']
                    delta = scenario_data['delta']
                    delta_abs = scenario_data['delta_abs']

                    out_path = os.path.join(scenario_result_dirs[scenario_name], f'{snp}.tsv')
                    scenario_df.to_csv(out_path, sep='\t', index=False)

                    if snp == most_significant_SNP:
                        boxplot_data_by_scenario[scenario_name].append({
                            'label': label,
                            'ancestry': 'AMR' if ancestry == 'NAT' else ancestry,
                            'snp': snp,
                            'delta': delta.tolist(),
                            'delta_abs': delta_abs.tolist(),
                            'color': COLOR_MAP.get(ancestry, 'gray'),
                        })
                        prob_path = os.path.join(scenario_prob_dirs[scenario_name], f'{snp}.tsv')
                        scenario_df.to_csv(prob_path, sep='\t', index=False)

    return boxplot_data_by_scenario


def plot_filtered_boxplot(boxplot_data, data_key, ylabel, title, output_path):
    if not boxplot_data:
        print(f"Nessun dato da plottare per {title}.")
        return

    filtered = [b for b in boxplot_data if len(b[data_key]) > 0]
    if not filtered:
        print(f"Nessun dato valido per {title}.")
        return

    filtered.sort(key=lambda x: np.median(x[data_key]))
    spacing = 2
    box_positions = [i * spacing for i in range(1, len(filtered) + 1)]

    fig, ax = plt.subplots(figsize=(max(14, len(filtered) * 0.75), 10))
    box_data = [b[data_key] for b in filtered]
    bp = ax.boxplot(box_data, positions=box_positions, patch_artist=True)

    for patch, b in zip(bp['boxes'], filtered):
        ancestry = b['ancestry']
        patch.set_facecolor(COLOR_MAP.get('NAT') if ancestry == 'AMR' else COLOR_MAP.get(ancestry, 'gray'))

    x_labels = [b['label'] for b in filtered]
    ax.set_xticks(box_positions)
    ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=12)

    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.set_ylim(-1 if ylabel == "Delta Probabilities" else 0, 1.1)

    max_values = [max(b[data_key]) for b in filtered]
    for i, (pos, b, max_val) in enumerate(zip(box_positions, filtered, max_values)):
        y_pos = min(max_val + (0.30 if i % 2 == 0 else 0.05), 1.08)
        ax.text(pos, y_pos, b['snp'], ha='center', va='bottom', fontsize=11)

    ancestry_set = sorted(set(b['ancestry'] for b in filtered))
    legend_elements = []
    for ancestry in ancestry_set:
        color = COLOR_MAP.get('NAT') if ancestry == 'AMR' else COLOR_MAP.get(ancestry, 'gray')
        legend_elements.append(Patch(facecolor=color, label=ancestry))

    ax.legend(
        handles=legend_elements,
        title="Global Ancestry",
        loc='center left',
        bbox_to_anchor=(1.01, 0.5),
        fontsize=11,
        title_fontsize=12,
    )

    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    result_folder = sys.argv[1]

    boxplot_data_by_scenario = data_processing(result_folder)

    plots_folder_general = os.path.join(result_folder, 'plots')

    for scenario_name, boxplot_data in boxplot_data_by_scenario.items():
        scenario_conf = SCENARIOS[scenario_name]
        scenario_plot_folder = os.path.join(plots_folder_general, scenario_name)
        os.makedirs(scenario_plot_folder, exist_ok=True)

        delta_title = "Boxplot of Delta Probabilities by Ancestry for UKBB"
        delta_path = os.path.join(
            scenario_plot_folder,
            f"delta_probabilities_UKBB_{scenario_conf['file_suffix']}.pdf",
        )
        plot_filtered_boxplot(boxplot_data, 'delta', "Delta Probabilities", delta_title, delta_path)

        abs_delta_title = "Boxplot of |Delta Probabilities| by Ancestry for UKBB"
        abs_delta_path = os.path.join(
            scenario_plot_folder,
            f"abs_delta_probabilities_UKBB_{scenario_conf['file_suffix']}.pdf",
        )
        plot_filtered_boxplot(boxplot_data, 'delta_abs', "|Delta Probabilities|", abs_delta_title, abs_delta_path)
