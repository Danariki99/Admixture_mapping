import subprocess

# Definizione del comando PLINK
plink_command = [
    "/private/home/rsmerigl/plink2",
    "--vcf", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb_filtered.vcf.gz",
    "--indep-pairwise", "50", "5", "0.5",
    "--out", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_pruned_files/output"
]

subprocess.run(plink_command)

