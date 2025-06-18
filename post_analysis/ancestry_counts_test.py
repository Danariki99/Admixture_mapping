import os
import pandas as pd
import numpy as np
import snputils as su
import time
import sys
import glob

# Usage: python generate_ancestry_counts_from_msp_folder.py <path_to_msp_folder>

if len(sys.argv) != 3:
    print("Usage: python generate_ancestry_counts_from_msp_folder.py <path_to_msp_folder>")
    sys.exit(1)

msp_folder = sys.argv[1]
result_folder = sys.argv[2]
output_folder = os.path.join(result_folder, 'ancestry_counts')
os.makedirs(output_folder, exist_ok=True)

# Define ancestry labels
ancestry_labels = ['AFR', 'EAS', 'EUR']
df_counts = None  # Will be initialized after reading first file

print(f"Processing all .msp files in folder: {msp_folder}")

# Find all chr_*.msp files
msp_files = sorted(glob.glob(os.path.join(msp_folder, 'chr*.msp')))

if len(msp_files) == 0:
    print("❌ No .msp files found in folder!")
    sys.exit(1)

start_time = time.time()

for msp_file_path in msp_files:
    print(f"\n→ Processing {msp_file_path}")
    laiobj = su.MSPReader(msp_file_path).read()

    # Initialize DataFrame on first file
    if df_counts is None:
        df_counts = pd.DataFrame(columns=['#IID'] + ancestry_labels)
        df_counts['#IID'] = laiobj.samples
        df_counts.iloc[:, 1:] = 0

    # Prepare count matrix
    n_wind, n_haplotypes = laiobj.lai.shape
    n_samp = len(laiobj.samples)
    hap_to_sample = np.arange(n_haplotypes) // 2
    ancestry_codes = laiobj.lai
    ancestry_to_col = {v: i for i, v in enumerate(df_counts.columns[1:])}
    temp_counts = np.zeros((n_samp, len(ancestry_labels)), dtype=int)

    for wind in range(n_wind):
        ancestry_labels_row = [laiobj.ancestry_map.get(str(code), None) for code in ancestry_codes[wind, :]]
        ancestry_indices = [ancestry_to_col.get(label, None) for label in ancestry_labels_row]
        valid_indices = [(hap_idx, col) for hap_idx, col in enumerate(ancestry_indices) if col is not None]
        for hap_idx, col in valid_indices:
            temp_counts[hap_to_sample[hap_idx], col] += 1

    df_counts.iloc[:, 1:] += temp_counts

total_time = round(time.time() - start_time, 2)
print(f"\n✅ All MSP files processed in {total_time} seconds.")

# Save final counts
output_file = os.path.join(output_folder, f'final_counts_aggregated.csv')
df_counts.to_csv(output_file, index=False)
print(f"📁 Saved aggregated counts to: {output_file}")
