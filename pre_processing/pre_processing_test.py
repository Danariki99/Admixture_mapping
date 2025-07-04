import os
import sys
from pre_processing_functions_test import vcf_creation
from pre_processing_functions_test import vcf_merging

if __name__ == '__main__':
    # Take dataset variable from terminal
    msp_folder = sys.argv[1]
    result_folder = sys.argv[2]
    output_folder = os.path.join(result_folder, 'vcf_files')
    tmp_dir = os.path.join(result_folder, 'tmp')

    os.makedirs(output_folder, exist_ok=True)
    os.makedirs(tmp_dir, exist_ok=True)

    # Define the ancestry map
    ancestry_map = {
        '0': 'AFR',
        '1': 'EAS',
        '2': 'EUR',
    }

    # Start the VCF creation
    print('Starting the VCF creation')
    vcf_creation(ancestry_map, msp_folder, tmp_dir)

    vcf_merging(ancestry_map, output_folder, tmp_dir)