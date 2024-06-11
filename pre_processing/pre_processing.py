import os
from pre_processing_functions import vcf_creation
from pre_processing_functions import vcf_merging

if __name__ == '__main__':

    # Define the ancestry map
    ancestry_map = {
     '0': 'AFR'
    }

    #Define the folders
    msp_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/ukbb'
    output_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/ukbb'

    # Start the VCF creation
    print('Starting the VCF creation')
    #vcf_creation(ancestry_map, msp_folder, output_folder)

    # Start the VCF merging
    print('starting the VCF merging')
    vcf_merging(ancestry_map, output_folder)