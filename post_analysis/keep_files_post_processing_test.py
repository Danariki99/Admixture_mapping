import os
import pandas as pd
import sys

# Percorsi coerenti con la tua pipeline
result_folder = sys.argv[1]
data_folder = sys.argv[2]

keep_file_global = os.path.join(data_folder,'input.keep')
keep_folder = os.path.join(result_folder, 'keep_files')
output_folder = os.path.join(result_folder, 'keep_files_processed')

# Crea la cartella di output se non esiste
os.makedirs(output_folder, exist_ok=True)

# Carica il file globale
keep_global_df = pd.read_csv(keep_file_global, sep='\t', header=None, names=['FID', 'IID'])

# Cicla su tutti i keep ancestry-specifici
for keep_file in os.listdir(keep_folder):
    keep_path = os.path.join(keep_folder, keep_file)

    # Carica file ancestry
    ancestry_df = pd.read_csv(keep_path, sep='\t', header=None, names=['FID', 'IID'])

    # Filtro: solo individui presenti nel keep globale
    filtered_df = ancestry_df[ancestry_df['IID'].isin(keep_global_df['IID'])]

    # Salva il file filtrato
    output_path = os.path.join(output_folder, keep_file)
    filtered_df.to_csv(output_path, sep='\t', index=False, header=False)

    print(f"Saved: {output_path} ({len(filtered_df)-1} individuals)")

print("✅ All ancestry-specific keep files filtered.")
