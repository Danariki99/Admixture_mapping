import os
import subprocess
import sys

# === Controlla argomento VCF ===
if len(sys.argv) != 2:
    print("Usage: python run_pca_for_ancestries.py <vcf_file.vcf.gz>")
    sys.exit(1)

vcf_file = sys.argv[1]

# === File paths ===
keep_path = './results/ancestry_keep/keep_files_processed'
output_folder = './results/pca/PCA_res'
exclude_snps_file = './results/pca/PCA_pruned_files/output.prune.out'

# Assicurati che la cartella di output esista
os.makedirs(output_folder, exist_ok=True)

# Lista dei file keep
keep_files = os.listdir(keep_path)

for keep_filename in keep_files:
    keep_file = os.path.join(keep_path, keep_filename)
    ancestry = keep_filename.split('.')[0].split('_')[0]

    output_prefix = os.path.join(output_folder, f"PCA_{ancestry}.out")

    plink_command = [
        "plink2",
        "--vcf", vcf_file,
        "--keep", keep_file,
        "--exclude", exclude_snps_file,
        "--pca", "approx", "10",
        "--out", output_prefix
    ]

    print(f"▶ Running PCA for ancestry: {ancestry}")
    subprocess.run(plink_command, check=True)

print("✅ PCA complete for all ancestries.")
