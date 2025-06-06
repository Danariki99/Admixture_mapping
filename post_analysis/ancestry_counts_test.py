import os
import pandas as pd
import numpy as np
import snputils as su
import time
import sys

# Usage: python generate_ancestry_counts_from_msp.py <path_to_msp_file>

if len(sys.argv) != 2:
    print("Usage: python generate_ancestry_counts_from_msp.py <path_to_msp_file>")
    sys.exit(1)

msp_file_path = sys.argv[1]
msp_file_name = os.path.basename(msp_file_path)
msp_id = os.path.splitext(msp_file_name)[0].replace('.msp', '')  # identifier for naming

output_folder = './results/ancestry_counts'
os.makedirs(output_folder, exist_ok=True)

# Define ancestry labels
ancestry_labels = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']
df_counts = pd.DataFrame(columns=['#IID'] + ancestry_labels)

print('Processing', msp_file_path)

# Load the msp file
start_time = time.time()
laiobj = su.MSPReader(msp_file_path).read()

# Initialize the DataFrame
df_counts['#IID'] = laiobj.samples
df_counts.iloc[:, 1:] = 0

# Count
n_wind, n_haplotypes = laiobj.lai.shape
n_samp = len(laiobj.samples)

print('Start counting...')
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

print('Time to process the MSP file:', round(time.time() - start_time, 2), 'seconds')

# Save final counts
output_file = os.path.join(output_folder, f'final_counts_{msp_id}.csv')
df_counts.to_csv(output_file, index=False)
print(f"Saved: {output_file}")
