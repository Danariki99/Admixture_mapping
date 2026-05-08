import os
import argparse
import subprocess

def extract_chromosomes(vcf_file):
    cmd = ["bcftools", "query", "-f", "%CHROM\n", vcf_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    chromosomes = sorted(set(result.stdout.strip().split("\n")))
    return chromosomes

def split_by_chromosome(vcf_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    chromosomes = extract_chromosomes(vcf_file)

    print(f"Found chromosomes: {chromosomes}")

    for chrom in chromosomes:
        output_file = os.path.join(output_dir, f"{chrom}.vcf.gz")
        cmd = [
            "bcftools", "view", "-r", chrom,
            "-Oz", "-o", output_file,
            vcf_file
        ]
        print(f"Extracting chr {chrom} → {output_file}")
        subprocess.run(cmd, check=True)

    print("✅ Splitting complete.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split multi-chromosome VCF into one VCF per chromosome.")
    parser.add_argument("--results", required=True, help="output results folder")
    parser.add_argument("--data", required=True, help="input data folder")
    
    args = parser.parse_args()  
    vcf_file = os.path.join(args.data, "input.vcf.gz")
    output_dir = os.path.join(args.results, "vcf_folder")
    split_by_chromosome(vcf_file, output_dir)
    print(output_dir)
