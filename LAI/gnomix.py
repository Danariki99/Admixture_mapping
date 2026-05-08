import os
import subprocess
import argparse

# Config
OUTPUT_DIR = "./results/msp_folder"

def run_gnomix(query_path, chrom, GNOMIX_REF):
    output_path = os.path.join(OUTPUT_DIR, f"chr_{chrom}.msp")
    command = [
        "python", "../gnomix/gnomix.py",
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
    parser.add_argument("--gnomix_ref", required=True, help="Path to the Gnomix reference panel")
    args = parser.parse_args()

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for filename in sorted(os.listdir(args.vcf_folder)):
        if filename.startswith("chr") and filename.endswith(".vcf.gz"):
            chrom = filename.replace("chr", "").replace(".vcf.gz", "")
            query_path = os.path.join(args.vcf_folder, filename)
            print(f"Running Gnomix on chr {chrom}")
            try:
                run_gnomix(query_path, chrom, args.gnomix_ref)
                print(f"✅ Done chr {chrom}")
                print(OUTPUT_DIR)
            except subprocess.CalledProcessError:
                print(f"❌ Failed chr {chrom}")
