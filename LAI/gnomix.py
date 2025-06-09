import os
import subprocess
import argparse

# Config
GNOMIX_REF = "../data/panel.vcf.gz"
OUTPUT_DIR = "../results/msp_folder"

def run_gnomix(query_path, chrom):
    output_path = os.path.join(OUTPUT_DIR, f"chr_{chrom}.msp")
    command = [
        "gnomix",
        "--train1", "0.8",
        "--train2", "0.15",
        "--val", "0.05",
        "--gens", "0,2,4,6,8,12,16,24",
        "--r_admixed", "1.5",
        "--window_size", "0.2",
        "--smoother_size", "115",
        "--context_ratio", "2.0",
        "--reference_panel", GNOMIX_REF,
        "--input", query_path,
        "--output", output_path,
        "--n-threads=32"
    ]
    subprocess.run(command, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch run Gnomix on chr_*.vcf.gz files.")
    parser.add_argument("--vcf_folder", required=True, help="Folder containing per-chromosome VCFs")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for filename in sorted(os.listdir(args.vcf_folder)):
        if filename.startswith("chr_") and filename.endswith(".vcf.gz"):
            chrom = filename.split("_")[1].split(".")[0]
            query_path = os.path.join(args.vcf_folder, filename)
            print(f"Running Gnomix on chr {chrom}")
            try:
                run_gnomix(query_path, chrom)
                print(f"✅ Done chr {chrom}")
            except subprocess.CalledProcessError:
                print(f"❌ Failed chr {chrom}")
