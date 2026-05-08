import os
import subprocess
import argparse

# === CONFIGURAZIONI FISSE ===
PHASE = "1"  
GENETIC_MAP_DIR = "data/maps"  # deve contenere chrN.gmap
SAMPLE_MAP_FILE = "data/sample_map_filtered.txt"  # file con assegnazione ancestry dei samples
OUTPUT_DIR = "/private/groups/ioannidislab/smeriglio/tests_files/gnomix_models"

def train_gnomix_model(test_vcf, chrom, gmap_file, panel_vcf, sample_map):
    output_basename = os.path.join(OUTPUT_DIR, f"{chrom}_model")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    current_dir = os.getcwd()
    gmap_file = os.path.abspath(gmap_file)
    panel_vcf = os.path.abspath(panel_vcf)
    sample_map = os.path.abspath(sample_map)
    print(f"Using gmap: {gmap_file}")
    print(f"Using panel: {panel_vcf}")
    print(f"Using sample map: {sample_map}")
    os.chdir('../gnomix')  # Assicurati di essere nella cartella corretta
    print(os.getcwd())  # Debug: Stampa la cartella corrente

    command = [
        "python3", "-u", "gnomix.py",
        test_vcf,
        output_basename,
        chrom,
        PHASE,
        gmap_file,
        panel_vcf,
        sample_map
    ]

    print(f"➡️  Training model for {chrom}...")
    subprocess.run(command, check=True)
    print(f"✅ Model {chrom} saved in {output_basename}")
    os.chdir(current_dir)  # Torna alla cartella originale

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train one Gnomix model per chromosome")
    parser.add_argument("--vcf_folder", required=True, help="Folder with per-chromosome test VCFs (chrN.vcf.gz)")
    parser.add_argument("--panel_folder", required=True, help="Folder with per-chromosome reference VCFs")
    args = parser.parse_args()
    print(f"Using VCF folder: {args.vcf_folder}")

    for fname in sorted(os.listdir(args.vcf_folder)):
        if fname.startswith("chr") and fname.endswith(".vcf.gz"):
            chrom = fname.replace(".vcf.gz", "")
            test_vcf = os.path.join(args.vcf_folder, fname)
            panel_vcf = os.path.join(args.panel_folder, fname)
            gmap_file = os.path.join(GENETIC_MAP_DIR, f"{chrom}.gmap")

            if not os.path.exists(panel_vcf):
                print(f"❌ Missing panel for {chrom}, skipping.")
                continue
            if not os.path.exists(gmap_file):
                print(f"❌ Missing gmap for {chrom}, skipping.")
                continue
            if not os.path.exists(SAMPLE_MAP_FILE):
                print(f"❌ Missing sample_map file, aborting.")
                break

            try:
                train_gnomix_model(test_vcf, chrom, gmap_file, panel_vcf, SAMPLE_MAP_FILE)
            except subprocess.CalledProcessError:
                print(f"❌ Failed training model for {chrom}")
