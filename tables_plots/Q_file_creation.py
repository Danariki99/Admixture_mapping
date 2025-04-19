import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

def reorder_admixture(Q_mat):
    n_samples, K = Q_mat.shape
    row_groups = []
    boundary_list = [0]
    for k in range(K):
        rows_k = np.where(np.argmax(Q_mat, axis=1) == k)[0]
        rows_k_sorted = rows_k[np.argsort(Q_mat[rows_k, k])[::-1]]
        row_groups.append(rows_k_sorted)
        boundary_list.append(boundary_list[-1] + len(rows_k_sorted))

    row_order = np.concatenate(row_groups)
    Q_mat_sorted = Q_mat[row_order, :]
    return Q_mat_sorted, row_order, boundary_list, None  # no col_order

def plot_admixture(ax, Q_mat_sorted, boundary_list, col_order=None, colors=None):
    n_samples, K = Q_mat_sorted.shape
    Q_cum = np.cumsum(Q_mat_sorted, axis=1)
    x_vals = np.arange(n_samples)
    step_mode = "post"
    ax.step(x_vals, Q_cum, linewidth=0.0, where=step_mode)
    for j in range(K):
        c = colors[j] if colors is not None else None
        if j == 0:
            ax.fill_between(x_vals, 0, Q_cum[:, j], step=step_mode, color=c)
        else:
            ax.fill_between(x_vals, Q_cum[:, j - 1], Q_cum[:, j], step=step_mode, color=c)
    for boundary in boundary_list:
        ax.axvline(boundary, color='black', ls='--', lw=1.0)
    ax.set_xlim(0, n_samples - 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Samples")
    ax.set_ylabel("Ancestry Proportion")

def create_q_file(input_folder, output_folder, keep_file_path):
    keep_file = pd.read_csv(keep_file_path, sep="\t")
    keep_ids = keep_file["#IID"]
    files_list = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith(".csv")]
    q_file = pd.read_csv(files_list[0])
    for file in files_list[1:]:
        temp_file = pd.read_csv(file)
        q_file.iloc[:, 1:] += temp_file.iloc[:, 1:]
    q_file = q_file.rename(columns={"NAT": "AMR"})
    q_file = q_file[q_file["#IID"].isin(keep_ids)]
    id_col = q_file["#IID"]
    ancestry_vals = q_file.drop(columns=["#IID"])
    ancestry_proportions = ancestry_vals.div(ancestry_vals.sum(axis=1), axis=0)
    q_file = pd.concat([id_col, ancestry_proportions], axis=1)
    os.makedirs(output_folder, exist_ok=True)
    q_file_path = os.path.join(output_folder, "output_normalized_with_ids.Q")
    q_file.to_csv(q_file_path, index=False)
    return q_file_path

def plot_q_file(q_file_path, output_folder, ancestry_color_map):
    q_df = pd.read_csv(q_file_path)
    ancestry_labels = q_df.columns[1:].tolist()
    Q_mat = q_df.drop(columns=["#IID"]).to_numpy()
    Q_sorted, row_order, boundary_list, _ = reorder_admixture(Q_mat)
    colors_for_plot = [ancestry_color_map[anc] for anc in ancestry_labels]

    fig, ax = plt.subplots(figsize=(14, 5))
    plot_admixture(ax, Q_sorted, boundary_list, col_order=None, colors=colors_for_plot)

    handles = [plt.Rectangle((0, 0), 1, 1, color=ancestry_color_map[anc]) for anc in ancestry_labels]
    ax.legend(handles, ancestry_labels, loc='upper center', bbox_to_anchor=(0.5, -0.12),
              ncol=len(ancestry_labels), frameon=False)

    plot_path = os.path.join(output_folder, "admixture_plot_with_legend.png")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to: {plot_path}")

if __name__ == "__main__":
    input_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/counts'
    output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/'
    keep_file_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt'

    ancestry_color_map = {
        'AFR': '#1f77b4',  # blue
        'AHG': '#ff7f0e',  # orange
        'EAS': '#2ca02c',  # green
        'EUR': '#d62728',  # red
        'AMR': '#8c564b',  # brown
        'OCE': '#9467bd',  # purple
        'SAS': '#e377c2',  # pink
        'WAS': '#7f7f7f'   # gray
    }

    q_file_path = create_q_file(input_folder, output_folder, keep_file_path)
    plot_q_file(q_file_path, output_folder, ancestry_color_map)
