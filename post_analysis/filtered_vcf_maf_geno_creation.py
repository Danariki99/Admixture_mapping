import subprocess

# Definizione del comando PLINK
plink_command = [
    "/private/home/rsmerigl/plink2",
    "--vcf", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb.vcf.gz",
    "--maf", "0.05",
    "--geno", "0.1",
    "--recode", "vcf",
    "--keep", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt",
    "--out", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb_filtered"
]

subprocess.run(plink_command)

