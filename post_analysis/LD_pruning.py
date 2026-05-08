import os
import subprocess

os.makedirs('/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_pruned_files', exist_ok=True)

plink_command = [
    "/private/home/rsmerigl/plink2",
    "--vcf", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb_filtered.vcf.gz",
    "--keep", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/kcutoff_177/kcutoff_177.king.cutoff.in.id",
    "--indep-pairwise", "50", "5", "0.5",
    "--out", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_pruned_files/output"
]

subprocess.run(plink_command)

