import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# ============================
# Reorder stile "David"
# ============================
def reorder_admixture(Q_mat, use_secondary=True):
    """
    - Ordina le colonne (ancestry) per media decrescente.
    - Ordina le righe (individui) per ancestry dominante,
      con tie-break sulla seconda ancestry se use_secondary=True.
    """
    n_samples, K = Q_mat.shape

    # 1) Ordina colonne per media
    col_means = Q_mat.mean(axis=0)
    col_order = np.argsort(col_means)[::-1]
    Qc = Q_mat[:, col_order]

    # 2) Ordina righe per ancestry dominante
    row_groups = []
    boundary_list = [0]
    argmax_all = np.argmax(Qc, axis=1)

    for k in range(K):
        rows_k = np.where(argmax_all == k)[0]
        if rows_k.size == 0:
            boundary_list.append(boundary_list[-1])
            continue

        if use_secondary:
            top1 = Qc[rows_k, k]
            other = np.delete(Qc[rows_k], k, axis=1)
            top2 = other.max(axis=1)
            # ordina per top1 desc, poi top2 desc
            order_k = np.lexsort((-top2, -top1))
            rows_k_sorted = rows_k[order_k]
        else:
            rows_k_sorted = rows_k[np.argsort(Qc[rows_k, k])[::-1]]

        row_groups.append(rows_k_sorted)
        boundary_list.append(boundary_list[-1] + rows_k_sorted.size)

    if row_groups:
        row_order = np.concatenate(row_groups)
    else:
        row_order = np.arange(n_samples)

    Q_sorted = Qc[row_order, :]
    return Q_sorted, row_order, boundary_list, col_order


# ============================
# Crea Q file (normalizzato "alla David")
# ============================
def create_q_file(input_folder, output_folder, keep_file_path):
    # keep file (IID da tenere)
    keep_file = pd.read_csv(keep_file_path, sep="\t")
    keep_ids = keep_file["#IID"].astype(str)

    # lista di file di conteggi (.csv)
    files_list = [
        os.path.join(input_folder, f)
        for f in os.listdir(input_folder)
        if f.endswith(".csv")
    ]
    if not files_list:
        raise FileNotFoundError(f"Nessun .csv trovato in {input_folder}")

    # somma di tutte le matrici di conteggi
    q_file = pd.read_csv(files_list[0])
    for file in files_list[1:]:
        temp_file = pd.read_csv(file)
        q_file.iloc[:, 1:] += temp_file.iloc[:, 1:]

    # rinomina NAT -> AMR (come nel tuo codice)
    if "NAT" in q_file.columns:
        q_file = q_file.rename(columns={"NAT": "AMR"})

    # filtra per keep
    q_file["#IID"] = q_file["#IID"].astype(str)
    q_file = q_file[q_file["#IID"].isin(keep_ids)].copy()

    # separa colonne
    id_col = q_file["#IID"]
    ancestry_vals = q_file.drop(columns=["#IID"])

    # normalizzazione per riga (somme = 1, stile David)
    row_sums = ancestry_vals.sum(axis=1)
    nonzero = row_sums > 0
    ancestry_proportions = ancestry_vals.copy()
    ancestry_proportions.loc[nonzero] = ancestry_vals.loc[nonzero].div(
        row_sums[nonzero], axis=0
    )

    # ricompone Q file
    q_file_norm = pd.concat([id_col.reset_index(drop=True),
                             ancestry_proportions.reset_index(drop=True)], axis=1)

    os.makedirs(output_folder, exist_ok=True)
    q_file_path = os.path.join(output_folder, "output_normalized_with_ids.Q")
    q_file_norm.to_csv(q_file_path, index=False)
    print(f"Q file saved to: {q_file_path}")
    return q_file_path


# ============================
# Plot stile "David"
# ============================
def plot_q_file(q_file_path, output_folder, ancestry_color_map):
    q_df = pd.read_csv(q_file_path)

    # colonne ancestry (tutte tranne #IID)
    ancestry_labels = q_df.columns[1:].tolist()
    Q_mat = q_df[ancestry_labels].to_numpy()

    # rinormalizza per sicurezza (somme ~1)
    row_sums = Q_mat.sum(axis=1, keepdims=True)
    nonzero = row_sums.squeeze() > 0
    Q_mat[nonzero] = Q_mat[nonzero] / row_sums[nonzero]

    # reorder alla David (colonne + righe)
    Q_sorted, row_order, boundary_list, col_order = reorder_admixture(Q_mat, use_secondary=True)

    # ancestry e colori nell'ordine riordinato
    ancestry_labels_ordered = [ancestry_labels[i] for i in col_order]
    colors_for_plot = [ancestry_color_map.get(anc, "#CCCCCC") for anc in ancestry_labels_ordered]

    # cumulata per stack
    Q_cum = np.cumsum(Q_sorted, axis=1)
    n_samples, K = Q_sorted.shape

    # edges per fill_between step="post" (stile David)
    x_edges = np.arange(n_samples + 1)
    Q_pad = np.vstack([Q_cum, Q_cum[-1]])

    fig, ax = plt.subplots(figsize=(18, 4))

    for j in range(K):
        lower = Q_pad[:, j - 1] if j > 0 else np.zeros(n_samples + 1)
        upper = Q_pad[:, j]
        ax.fill_between(
            x_edges,
            lower,
            upper,
            step="post",
            color=colors_for_plot[j],
            linewidth=0,
        )

    # linee verticali tra cluster di ancestry dominante
    # (boundary_list è in termini di numero di individui)
    for b in boundary_list[1:-1]:  # salta 0 e l'ultima
        ax.axvline(b, color="black", ls="--", lw=1.0)

    # stile più "clean", tipo admixture plot classico
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_ylim(0, 1)

    ax.set_title("Admixture Plot for UKBB", fontsize=16)

    # legenda con ancestries nell'ordine riordinato
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors_for_plot[j])
        for j in range(K)
    ]
    ax.legend(
        handles,
        ancestry_labels_ordered,
        loc='upper center',
        bbox_to_anchor=(0.5, -0.10),
        ncol=len(ancestry_labels_ordered),
        frameon=False,
        prop={'size': 10},
    )

    plot_path = os.path.join(output_folder, "admixture_plot_with_legend.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Plot saved to: {plot_path}")


# ============================
# MAIN
# ============================
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
