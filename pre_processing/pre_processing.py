import os
import sys
from pre_processing_functions import vcf_creation
from pre_processing_functions import vcf_merging

if __name__ == '__main__':
    # Take dataset variable from terminal
    dataset = sys.argv[1]

    # Define the ancestry map
    ancestry_map = {
        '0': 'AFR',
        '1': 'AHG',
        '2': 'EAS',
        '3': 'EUR',
        '4': 'NAT',
        '6': 'SAS',
        '7': 'WAS'
    }

    # Define the folders using the dataset variable
    msp_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/msp_files/{dataset}'
    tmp_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/tmp'
    output_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files/{dataset}'

    # Start the VCF creation
    print('Starting the VCF creation')
    vcf_creation(ancestry_map, msp_folder, tmp_folder)

    # Start the VCF merging
    print('starting the VCF merging')
    vcf_merging(ancestry_map, output_folder, tmp_folder)