import os
import subprocess
import sys

# === Controlla argomento VCF ===
if len(sys.argv) != 3:
    print("Usage: python run_pca_for_ancestries.py <results_folder> <data_folder>")
    sys.exit(1)

result_folder = sys.argv[1]
data_folder = sys.argv[2]

vcf_file = os.path.join(data_folder, "input.vcf.gz")


# === File paths ===
keep_path = os.path.join(result_folder, "keep_files_processed")
output_folder = os.path.join(result_folder, "PCA_files", "PCA_res")
exclude_snps_file = os.path.join(result_folder, "PCA_files", "PCA_pruned_files", "output", "output.prune.out")

# Assicurati che la cartella di output esista
os.makedirs(output_folder, exist_ok=True)

# Lista dei file keep
keep_files = os.listdir(keep_path)

for keep_filename in keep_files:
    keep_file = os.path.join(keep_path, keep_filename)
    ancestry = keep_filename.split('.')[0].split('_')[0]

    # Count number of individuals
    with open(keep_file, 'r') as f:
        n_lines = sum(1 for _ in f)

    if n_lines < 50:
        print(f"⚠️ Skipping ancestry '{ancestry}': only {n_lines} individuals")
        continue

    output_prefix = os.path.join(output_folder, f"PCA_{ancestry}.out")

    plink_command = [
        "/private/home/rsmerigl/plink2",
        "--vcf", vcf_file,
        "--keep", keep_file,
        "--exclude", exclude_snps_file,
        "--pca", "approx", "10",
        "--out", output_prefix
    ]

    print(f"▶ Running PCA for ancestry: {ancestry}")
    subprocess.run(plink_command, check=True)

print("✅ PCA complete for all ancestries.")
