import numpy as np
import pandas as pd
import os
import glob
import subprocess
import sys
#sys.path.append('/private/home/rsmerigl/codes/cleaned_codes')

from snputils import MSPReader
from snputils import AdmixtureMappingVCFWriter


def vcf_creation(ancestry_map, msp_folder, output_folder):

    msp_file_list = os.listdir(msp_folder)
    print(msp_file_list)

    for msp_file_name in msp_file_list:

        output_file = os.path.join(output_folder, msp_file_name.replace('.msp', '.vcf'))

        # read msp file 
        msp_file = os.path.join(msp_folder, msp_file_name)

        # Create an instance of MSPReader
        msp_reader = MSPReader(msp_file)

        # Read the msp file and create a LocalAncestryObject
        print(f"Reading {msp_file}")
        lai_object = msp_reader.read()

        admix_map = AdmixtureMappingVCFWriter(lai_object, output_file, ancestry_map)

        admix_map.write()

       
    return output_folder

def vcf_merging(ancestry_map, vcf_dir, tmp_dir):

    vcf_generic_filename = '*_#.vcf'

    vcf_generic_output_filename = 'ancestry_#.vcf'

    generic_output_file = os.path.join(vcf_dir, vcf_generic_output_filename)

    for ancestry in range(len(ancestry_map)):
        
        vcf_filenames = []

        output_file = generic_output_file.replace('#', ancestry_map[str(ancestry)])

        # Find all files with the current ancestry
        vcf_filenames = glob.glob(os.path.join(tmp_dir, vcf_generic_filename.replace('#', ancestry_map[str(ancestry)])))
        print(vcf_filenames)

        # Extract chromosome numbers and sort the filenames
        vcf_filenames.sort(key=lambda filename: int(filename.split('/')[-1].split('_')[0][3:]))

        print(vcf_filenames)

        concat_command = ['bcftools', 'concat'] + vcf_filenames + ['-o', output_file]

        # Run the command and check if it was successful
        result = subprocess.run(concat_command)

        if result.returncode == 0:  # If the command was successful
            # Delete the original files
            for filename in vcf_filenames:
                os.remove(filename)
        else:
            print(f"Error merging VCF files for ancestry {ancestry_map[str(ancestry)]}")