import subprocess
import os
import sys

if len(sys.argv) != 3:
    print("Usage: python generate_prune_file.py <vcf_file>")
    sys.exit(1)

result_folder = sys.argv[1]
data_folder = sys.argv[2]
vcf_file = os.path.join(data_folder, "input.vcf.gz")
output_folder = os.path.join(result_folder, "PCA_files/PCA_pruned_files/output")

os.makedirs(output_folder, exist_ok=True)

plink_command = [
    "/private/home/rsmerigl/plink2",
    "--vcf", vcf_file,
    "--set-all-var-ids", "@:#:$r:$a",
    "--new-id-max-allele-len", "500", "missing",
    "--rm-dup", "force-first",
    "--indep-pairwise", "50", "5", "0.2",
    "--out", os.path.join(output_folder, "output")
]

subprocess.run(plink_command)

