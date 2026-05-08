import os
import subprocess
import argparse

# Config
RFMIX_REF = "../data/panel.vcf.gz"
RFMIX_MAP = "../data/sample_map.txt"
RFMIX_GENMAP = "../data/genetic_map.txt"
OUTPUT_DIR = "../results/msp_folder"

def run_rfmix(query_path, chrom):
    output_prefix = os.path.join(OUTPUT_DIR, f"chr_{chrom}")
    command = [
        "rfmix",
        "-f", query_path,
        "-r", RFMIX_REF,
        "-m", RFMIX_MAP,
        "-g", RFMIX_GENMAP,
        "-o", output_prefix,
        f"--n-threads=8",
        f"--chromosome={chrom}"
    ]
    subprocess.run(command, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch run RFMix on chr_*.vcf.gz files.")
    parser.add_argument("--vcf_folder", required=True, help="Folder containing per-chromosome VCFs")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for filename in sorted(os.listdir(args.vcf_folder)):
        if filename.startswith("chr_") and filename.endswith(".vcf.gz"):
            chrom = filename.split("_")[1].split(".")[0]
            query_path = os.path.join(args.vcf_folder, filename)
            print(f"Running RFMix on chr {chrom}")
            try:
                run_rfmix(query_path, chrom)
                print(f"✅ Done chr {chrom}")
            except subprocess.CalledProcessError:
                print(f"❌ Failed chr {chrom}")
