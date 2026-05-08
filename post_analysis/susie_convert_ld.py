"""Converts _ld.npy files to gzipped text format readable by R."""
import numpy as np
import os

input_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/SuSiE_inputs'

for hit_dir in sorted(os.listdir(input_folder)):
    npy_file = os.path.join(input_folder, hit_dir, f'{hit_dir}_ld.npy')
    gz_file  = os.path.join(input_folder, hit_dir, f'{hit_dir}_ld.gz')

    if not os.path.exists(npy_file):
        continue
    if os.path.exists(gz_file):
        print(f"Already done: {hit_dir}")
        continue

    R = np.load(npy_file)
    np.savetxt(gz_file, R, fmt='%.6f')
    print(f"Converted: {hit_dir} ({R.shape[0]}x{R.shape[1]})")

print("Done.")
