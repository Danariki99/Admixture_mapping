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
    parser.add_argument("--vcf", required=True, help="Input VCF file (can be .vcf or .vcf.gz)")

    output_dir = 'results/vcf_folder'
    args = parser.parse_args()
    split_by_chromosome(args.vcf, output_dir)
    print(output_dir)
