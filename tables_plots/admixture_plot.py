import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({
    'font.size':          6,
    'axes.titlesize':     6,
    'axes.labelsize':     6,
    'xtick.labelsize':    5,
    'ytick.labelsize':    5,
    'pdf.fonttype':       42,
    'savefig.dpi':        600,
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Liberation Sans', 'Helvetica', 'DejaVu Sans'],
})

Q_FILE   = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files_old/ukbb/output_normalized_with_ids.Q'
OUT_BASE = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files_old/ukbb/admixture_files/run1/admixture_plot'

ANCESTRY_COLORS = {
    'EUR': '#d62728',
    'AFR': '#1f77b4',
    'EAS': '#2ca02c',
    'SAS': '#e377c2',
    'WAS': '#7f7f7f',
    'AMR': '#8c564b',
    'AHG': '#ff7f0e',
    'OCE': '#9467bd',
}


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

    col_means = Q_mat.mean(axis=0)
    col_order = np.argsort(col_means)[::-1]
    Qc = Q_mat[:, col_order]

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
            order_k = np.lexsort((-top2, -top1))
            rows_k_sorted = rows_k[order_k]
        else:
            rows_k_sorted = rows_k[np.argsort(Qc[rows_k, k])[::-1]]

        row_groups.append(rows_k_sorted)
        boundary_list.append(boundary_list[-1] + rows_k_sorted.size)

    row_order = np.concatenate(row_groups) if row_groups else np.arange(n_samples)
    Q_sorted = Qc[row_order, :]
    return Q_sorted, row_order, boundary_list, col_order


# Load Q file
df = pd.read_csv(Q_FILE)
df = df.rename(columns={'#IID': 'IID'})
df['IID'] = df['IID'].astype(str)
print(f'Q file loaded: {len(df)} samples')

# Crop the near-pure-EUR block: keep only samples with EUR ancestry <= EUR_MAX
EUR_MAX = 0.95
df = df[df['EUR'] <= EUR_MAX].reset_index(drop=True)
print(f'After EUR <= {EUR_MAX} crop: {len(df)} samples')

ancestry_labels = [c for c in df.columns if c != 'IID']
Q_mat = df[ancestry_labels].to_numpy()

Q_sorted, row_order, boundary_list, col_order = reorder_admixture(Q_mat, use_secondary=True)
ancestry_labels_ordered = [ancestry_labels[i] for i in col_order]
colors_for_plot = [ANCESTRY_COLORS.get(a, '#CCCCCC') for a in ancestry_labels_ordered]

Q_cum = np.cumsum(Q_sorted, axis=1)
n_samples, K = Q_sorted.shape

# edges per fill_between step="post"
x_edges = np.arange(n_samples + 1)
Q_pad = np.vstack([Q_cum, Q_cum[-1]])

# Plot
MM_TO_INCH = 1 / 25.4
fig, ax = plt.subplots(figsize=(100 * MM_TO_INCH, 45 * MM_TO_INCH))

for j in range(K):
    lower = Q_pad[:, j - 1] if j > 0 else np.zeros(n_samples + 1)
    upper = Q_pad[:, j]
    ax.fill_between(x_edges, lower, upper, step='post',
                     color=colors_for_plot[j], linewidth=0)

# linee verticali tra cluster di ancestry dominante
for b in boundary_list[1:-1]:
    ax.axvline(b, color='black', ls='--', lw=0.5)

for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(0, n_samples)
ax.set_ylim(0, 1)

ax.set_title('Admixture plot UKBB', fontsize=6)

handles = [plt.Rectangle((0, 0), 1, 1, color=colors_for_plot[j]) for j in range(K)]
ax.legend(
    handles, ancestry_labels_ordered,
    loc='upper center', bbox_to_anchor=(0.5, -0.05),
    ncol=K, frameon=False, fontsize=5,
    columnspacing=1.0, handletextpad=0.4,
)

plt.tight_layout()
os.makedirs(os.path.dirname(OUT_BASE), exist_ok=True)
plt.savefig(f'{OUT_BASE}.png', bbox_inches='tight', pad_inches=0.01)
plt.savefig(f'{OUT_BASE}.pdf', bbox_inches='tight', pad_inches=0.01, format='pdf')
print(f'Saved → {OUT_BASE}.png / .pdf')
