import os
import subprocess

ancestries = ['AFR', 'AHG', 'EAS', 'EUR', 'NAT', 'OCE', 'SAS', 'WAS']
generic_vcf_file = '/private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_*.vcf'
generic_output_file = '/private/groups/ioannidislab/smeriglio/LD/*_LD'

for ancestry in ancestries:

    print(f'Processing {ancestry}...')
    input_file = generic_vcf_file.replace('*', ancestry)
    output_file = generic_output_file.replace('*', ancestry)

    command = [
        "/private/home/rsmerigl/plink2",
        "--vcf",
        input_file,
        "--r2-unphased",
        "--ld-window-kb",
        "1000",
        "--ld-window-r2",
        "0.1",
        "--out",
        output_file
    ]

    subprocess.run(command, check=True)

