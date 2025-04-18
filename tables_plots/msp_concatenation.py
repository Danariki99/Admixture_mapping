from snputils.ancestry.io.local.read.msp import MSPReader
from snputils.ancestry.genobj.local import LocalAncestryObject

import glob
import numpy as np

# Percorsi a tutti i file .msp.tsv (uno per cromosoma)
paths = sorted(glob.glob("/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb/ukb_hap_chr*_v2_rfmix.msp.tsv"))

all_lai = []
all_physical_pos = []
all_centimorgan_pos = []
all_chromosomes = []
all_window_sizes = []

for path in paths:
    print(f"Processing {path}")
    # Usa correttamente il reader
    lai_obj = MSPReader(path).read()

    all_lai.append(lai_obj.lai)
    
    if lai_obj.physical_pos is not None:
        all_physical_pos.append(lai_obj.physical_pos)
    if lai_obj.centimorgan_pos is not None:
        all_centimorgan_pos.append(lai_obj.centimorgan_pos)
    if lai_obj.chromosomes is not None:
        all_chromosomes.append(lai_obj.chromosomes)
    if lai_obj.window_sizes is not None:
        all_window_sizes.append(lai_obj.window_sizes)

# Unisci la matrice LAI lungo le righe (n_windows total x n_haplotypes)
merged_lai = np.concatenate(all_lai, axis=0)

# Crea il nuovo oggetto LAI unificato
merged_lai_obj = LocalAncestryObject(
    haplotypes=lai_obj.haplotypes,  # Prendi da ultimo oggetto (identico in tutti)
    lai=merged_lai,
    samples=lai_obj.samples,
    ancestry_map=lai_obj.ancestry_map,
    window_sizes=np.concatenate(all_window_sizes) if all_window_sizes else None,
    centimorgan_pos=np.concatenate(all_centimorgan_pos) if all_centimorgan_pos else None,
    chromosomes=np.concatenate(all_chromosomes) if all_chromosomes else None,
    physical_pos=np.concatenate(all_physical_pos) if all_physical_pos else None
)

# Salva il file finale unificato
merged_lai_obj.save("/private/groups/ioannidislab/smeriglio/out_cleaned_codes/total_msp/merged_genome.msp")
