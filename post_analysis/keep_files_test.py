import os
import pandas as pd
import numpy as np

# Percorsi
counts_folder = './results/ancestry_counts'
output_folder = './results/keep_files'
os.makedirs(output_folder, exist_ok=True)

# Carica tutti i file counts
counts_files = [f for f in os.listdir(counts_folder) if f.startswith("final_counts_") and f.endswith(".csv")]

df_counts = None

for count_file in counts_files:
    print('Processing', count_file)
    df = pd.read_csv(os.path.join(counts_folder, count_file))
    df = df.set_index('#IID')
    df = df.dropna()
    df = df.astype(int)
    if df_counts is None:
        df_counts = df
    else:
        df_counts = df_counts.add(df, fill_value=0)

df_counts = df_counts.reset_index()

# Genera i keep file per ogni ancestry
for ancestry in df_counts.columns[1:]:
    max_ancestry_samples = df_counts.loc[
        df_counts[ancestry] == df_counts.iloc[:, 1:].max(axis=1), '#IID'
    ]
    keep_file = os.path.join(output_folder, f"{ancestry}_keep.txt")
    keep_df = pd.DataFrame(max_ancestry_samples, columns=['FID', 'IID'])
    keep_df['FID'] = keep_df['IID']  # FID = IID
    keep_df.to_csv(keep_file, index=False, header=False, sep='\t')
    print(f"Saved: {keep_file} ({len(keep_df)} individuals)")

print("Processing complete. Keep files created in:", output_folder)
