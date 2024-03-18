import os
import subprocess

num_chrom = 22

ancestry_map = {
 '0': 'AFR',
 '1': 'AHG',
 '2': 'EAS',
 '3': 'EUR',
 '4': 'NAT',
 '5': 'OCE',
 '6': 'SAS',
 '7': 'WAS'
}

vcf_file_dir = '/private/groups/ioannidislab/smeriglio/vcf_files/chr*'
vcf_generic_filename = 'ancestry_*.vcf'

output_dir = '/private/groups/ioannidislab/smeriglio/merged_vcfs'
generic_output_file = os.path.join(output_dir, vcf_generic_filename)

for ancestry in range(len(ancestry_map)):
    
    vcf_filenames = []

    output_file = generic_output_file.replace('*', ancestry_map[str(ancestry)])

    for chrom_n in range(num_chrom):
        chrom_n+=1
        current_chr_folder = vcf_file_dir.replace('*', str(chrom_n))

        vcf_filename = os.path.join(current_chr_folder, vcf_generic_filename.replace('*', ancestry_map[str(ancestry)]))

        vcf_filenames.append(vcf_filename)

    concat_command = ['bcftools', 'concat'] + vcf_filenames + ['-o', output_file]


    subprocess.run(concat_command)
    
    

