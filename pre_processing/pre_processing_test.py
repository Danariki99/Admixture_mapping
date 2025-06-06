import os
import sys
from pre_processing_functions import vcf_creation
from pre_processing_functions import vcf_merging

if __name__ == '__main__':
    # Take dataset variable from terminal
    msp_folder = sys.argv[1]
    output_folder = './data/vcf_files'

    os.mkdir(output_folder, exist_ok=True)

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

    # Start the VCF creation
    print('Starting the VCF creation')
    vcf_creation(ancestry_map, msp_folder, output_folder)

    vcf_merging(ancestry_map, output_folder, output_folder)