import os
import pandas as pd

# Percorso della directory dei file CSV con i conteggi
counts_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/counts/'

# Directory di output per i file .Q
output_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files/run1'

# Assicurarsi che la directory di output esista
os.makedirs(output_path, exist_ok=True)

# Elenco dei file CSV nella directory dei conteggi
counts_files = os.listdir(counts_path)

# Numero di ancestry groups
df_total = pd.read_csv(os.path.join(counts_path, counts_files[0]))
# Processare ogni file CSV
for count_file in counts_files[1:]:
    # Leggere il file CSV in un DataFrame
    df = pd.read_csv(os.path.join(counts_path, count_file))

    df_total = pd.concat([df_total, df])
    
# Creare un nuovo DataFrame vuoto con le stesse colonne del DataFrame originale
df_processed = pd.DataFrame(columns=df_total.columns)

# Copiare la prima colonna (assunta come identificatore non numerico) direttamente nel nuovo DataFrame
df_processed[df_total.columns[0]] = df_total[df_total.columns[0]]

# Processare il resto delle colonne (da 1 in poi) convertendole in float
for col in df.columns[1:]:
    df_processed[col] = df_total[col].astype(float)

# Normalizzare i dati elaborati (divisione per la somma delle righe)
df_processed.iloc[:, 1:] = df_processed.iloc[:, 1:].div(df_processed.iloc[:, 1:].sum(axis=1), axis=0)

# Verificare il numero di colonne di ancestry (escludendo l'identificatore)
num_columns = df_processed.shape[1] - 1  # Escludiamo la colonna degli identificatori

df_processed.to_csv(os.path.join(output_path, 'train_demo.all.Q'), index=False)