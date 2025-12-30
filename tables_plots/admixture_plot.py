import os
import pandas as pd
import matplotlib.pyplot as plt

# Percorso del file .Q (file di output)
file_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files/run1/train_demo.8.Q'

# Carica i dati dal file .Q
df = pd.read_csv(file_path)

# Numero di gruppi di ancestry (dovrebbe corrispondere al numero di colonne)
NUM_ANCESTRIES = 8

# Creiamo una variabile per le colonne dei gruppi di ancestry (tutte tranne la prima colonna)
ancestry_columns = df.columns[1:]

# Crea un grafico a barre impilate
fig, ax = plt.subplots(figsize=(12, 8))

# Accumulare i valori per ogni gruppo di ancestry e tracciarli come barre impilate
df[ancestry_columns].plot(kind='bar', stacked=True, ax=ax, colormap='tab20', width=1.0)

# Personalizzare il grafico
ax.set_title('Admixture Plot', fontsize=16)
ax.set_xlabel('Individuals (IID)', fontsize=12)
ax.set_ylabel('Proportion', fontsize=12)

# Personalizzazione della legenda
ax.legend(title="Ancestry Groups", bbox_to_anchor=(1.05, 1), loc='upper left')

# Aggiungere etichette ai tick dell'asse x (IID degli individui)
ax.set_xticks(range(len(df)))
ax.set_xticklabels(df[df.columns[0]], rotation=90, fontsize=10)  # Usa direttamente la colonna #IID

# Visualizzare il grafico
plt.tight_layout()
plt.savefig('/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files/run1/admixture_plot.png')
