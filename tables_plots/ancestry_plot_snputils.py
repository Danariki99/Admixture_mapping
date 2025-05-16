import os
from snputils.visualization.admixture_viz import pong_viz

# Configurazione
folder_runs = "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files"  # Cartella principale
output_dir = "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_visualizations"  # Cartella per i risultati
k_value = 8  # Numero di ancestry
run_prefix = "train_demo"  # Prefisso richiesto dallo script
runs = [1]  # Un singolo run
verbose = True

# Creazione struttura di esempio
run_folder = os.path.join(folder_runs, "run1")
os.makedirs(run_folder, exist_ok=True)

# File di input
input_file = "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/admixture_files/run1/train_demo.8.Q"  # Il tuo file originale
output_file = os.path.join(run_folder, f"{run_prefix}.{k_value}.Q")

# Se il file non esiste nella cartella `run1`, rinominalo e spostalo
if not os.path.exists(output_file):
    os.rename(input_file, output_file)

# Richiamo della funzione pong_viz
try:
    pong_viz(
        folder_runs=folder_runs,
        output_dir=output_dir,
        k=k_value,
        runs=runs,
        run_prefix=run_prefix,
        verbose=verbose
    )
    print(f"Elaborazione completata. Risultati salvati in: {output_dir}")
except Exception as e:
    print(f"Errore durante l'esecuzione di pong_viz: {e}")