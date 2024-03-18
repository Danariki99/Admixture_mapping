import os
import subprocess

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

generic_output_folder = '/private/groups/ioannidislab/smeriglio/output/no_BMI/output_ancestry_*'

generic_vcf_file = '/private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_*.vcf'

keep_file = '/private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt'

covar_file = '/private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar.phe'

pheno_folder = '/private/groups/ioannidislab/smeriglio/phe_files'

pheno_files_list = os.listdir(pheno_folder)

for i in range(len(ancestry_map)):
    ancestry = ancestry_map[str(i)]
    ancestry_output_folder = generic_output_folder.replace('*', ancestry)
    if not os.path.exists(ancestry_output_folder):
        os.makedirs(ancestry_output_folder)

    vcf_file = generic_vcf_file.replace('*', ancestry)

    for pheno_file in pheno_files_list:

        pheno = pheno_file.split('.')[0]

        current_pheno_file = os.path.join(pheno_folder, pheno_file)

        current_output_folder = os.path.join(ancestry_output_folder, pheno)

        if not os.path.exists(current_output_folder):
            os.makedirs(current_output_folder)
        
        output_file = os.path.join(current_output_folder, 'output')

        association_command = [
            '/private/home/rsmerigl/plink2',
            '--vcf', f'{vcf_file}',
            '--pheno', f'{current_pheno_file}',
            '--glm', 'firth-fallback', 'hide-covar',
            '--ci', '0.95',
            '--adjust',
            '--covar', f'{covar_file}',
            '--covar-variance-standardize',
            '--keep', f'{keep_file}',
            '--out', f'{output_file}',
            '--covar-col-nums', '2,3,8-17',

        ]

        subprocess.run(association_command, check=True)


        

