import os
import subprocess

# File paths
keep_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/ancestry_keep_files/ukbb/keep_files_processed'

output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_res'

keep_files = os.listdir(keep_path)

for keep_filname in keep_files:

    keep_file = os.path.join(keep_path, keep_filname)

    ancestry = keep_filname.split('.')[0].split('_')[0]
    
    plink_command = [
        "/private/home/rsmerigl/plink2", 
        "--vcf", "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/vcf_file/ukbb_filtered.vcf.gz",  
        "--keep", keep_file,
        "--exclude",  "/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_pruned_files/output.prune.out",
        "--pca", "approx", "10", 
        "--out", os.path.join(output_folder, f"PCA_{ancestry}.out")
    ]

    subprocess.run(plink_command)