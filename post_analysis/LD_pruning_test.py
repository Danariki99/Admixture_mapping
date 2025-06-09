import subprocess
import os
import sys

if len(sys.argv) != 2:
    print("Usage: python generate_prune_file.py <vcf_file>")
    sys.exit(1)

vcf_file = sys.argv[1]

# Assicurati che la cartella di output esista
output_folder = './results/pca/PCA_pruned_files'
os.makedirs(output_folder, exist_ok=True)

plink_command = [
    "plink2",
    "--vcf", vcf_file,
    "--indep-pairwise", "50", "5", "0.5",
    "--out", os.path.join(output_folder, "output")
]

subprocess.run(plink_command)

