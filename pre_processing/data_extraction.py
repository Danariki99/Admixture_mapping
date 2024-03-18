import subprocess
import pandas as pd
import os

def samples_number_phe(filename):
    #extract the number of samples
    N_command = f"grep -v '^#' {filename} | wc -l"
    N_output = subprocess.run(N_command, capture_output=True, text=True, shell=True)
    N = N_output.stdout.strip()
    return N

def samples_number_vcf(filename):
    # Extracting the number of samples from the VCF
    samples_cmd = f"bcftools query -l {filename} | wc -l"
    num_samples = int(subprocess.check_output(samples_cmd, shell=True).strip())
    return num_samples


def QC_unrelated_samples_phe(filename):


    usefull_samples_command = ['awk', 'NR > 1 && $3 != "DO_NOT_PASS_SQC" &&  $3 != "related" {print $2}', filename]


    usefull_samples_result = subprocess.run(usefull_samples_command, capture_output=True, text=True)
    usefull_samples_list = usefull_samples_result.stdout.strip().split('\n')
    return usefull_samples_list

def extract_sample_IDs_vcf(vcf_file, tmp_folder):

    # Create a temporary file path for storing the sample names
    tmp_sample_file = os.path.join(tmp_folder, 'sample_names.txt')

    # Run bcftools to extract sample names and save them to the temporary file
    command = f'bcftools query -l {vcf_file} > {tmp_sample_file}'
    subprocess.run(command, shell=True, check=True)

    # Read the sample IDs from the temporary file
    with open(tmp_sample_file, 'r') as f:
        sample_IDs = [line.strip() for line in f]

    # Remove the temporary file
    os.remove(tmp_sample_file)

    return sample_IDs

def extract_sample_IDs_phe(phe_file):
    command = ['awk', '{print $2}', phe_file]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode == 0:
        values_list = result.stdout.strip().split('\n')
        return values_list
    else:
        print("Error:", result.stderr)
        return []

def write_matching_FID_file(common_samples_gwas_vcf, output_filename):
    # Open the output file in write mode
    with open(output_filename, 'w') as file:
        # Write the header
        file.write('#IID\n')
        
        # Write the data from matching_FID and common_samples_gwas_vcf
        for sample in common_samples_gwas_vcf:
            file.write(f'{sample}\n')


if __name__ == "__main__":
    tmp_dir = os.getcwd()
    master_file = '/private/groups/ioannidislab/ukbb24983/phenotypedata/master_phe/master.20211020.phe'
    gwas_file = '/private/groups/ioannidislab/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
    vcf_file = '/private/groups/ioannidislab/smeriglio/vcf_files/chr21/ancestry_AFR.vcf'
    phe_file = '/private/groups/ioannidislab/ukbb24983/phenotypedata/extras/highconfidenceqc/current/phe/HC132.phe'

    master_sample_list = QC_unrelated_samples_phe(master_file)
    N_master = samples_number_phe(master_file)

    gwas_sample_list = QC_unrelated_samples_phe(gwas_file)
    N_gwas = samples_number_phe(gwas_file)

    vcf_sample_list = extract_sample_IDs_vcf(vcf_file, tmp_dir)
    N_vcf = samples_number_vcf(vcf_file)

    phe_samples_list = extract_sample_IDs_phe(phe_file)

    common_samples_master_gwas = list(set(master_sample_list) & set(gwas_sample_list))

    common_samples_gwas_vcf = list(set(master_sample_list) & set(vcf_sample_list))

    usefull_samples = list(set(common_samples_gwas_vcf) & set(phe_samples_list))

    print(f'the total number of samples in the master file {master_file} is {N_master},\nin the gwas file {gwas_file} is {N_gwas} and\nin the vcf file {vcf_file} is {N_vcf}\nin the phe file is {phe_file} is {len(phe_samples_list)}\nthe number of usefull samples in the master file is {len(master_sample_list)} and in the gwas file is {len(gwas_sample_list)}\nthe number of samples that are both in the master and in the gwas files is {len(common_samples_master_gwas)}\nthe number of samples that are both in the gwas and the vcf files are {len(common_samples_gwas_vcf)} \nthe number of usefull files compared also with the phe file is {len(usefull_samples)}')
    
    #create the keep_file:
    keep_filename = '/private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt'

    write_matching_FID_file(usefull_samples, keep_filename)