import pandas as pd

# Carica il file originale
file_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files/run1/train_demo.all.Q'  # Sostituisci con il percorso del tuo file
df = pd.read_csv(file_path)

# Rimuovi la colonna IID (non necessaria per la matrice Q)
df.drop(columns=['#IID'], inplace=True)

# Salva la matrice Q senza intestazione e con valori numerici
output_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files/run1/train_demo.8.Q'  # Sostituisci con il percorso di output desiderato
df.to_csv(output_path, index=False, header=False, sep=' ')
