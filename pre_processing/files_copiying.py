import shutil
import glob
import os

# Definisci la cartella di origine e la cartella di destinazione
source_folder = '/private/groups/ioannidislab/ukbb_analysis/rfmix/output'
destination_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb'

# Crea la cartella di destinazione se non esiste
os.makedirs(destination_folder, exist_ok=True)

# Trova tutti i file che corrispondono al modello
for file_name in glob.glob(os.path.join(source_folder, 'ukb_hap_chr*_v2_rfmix.msp.tsv')):
    # Sposta ogni file nella cartella di destinazione
    shutil.copy(file_name, destination_folder)