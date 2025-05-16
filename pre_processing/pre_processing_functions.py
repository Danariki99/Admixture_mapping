import numpy as np
import pandas as pd
import os
import glob
import subprocess
import sys
sys.path.append('/private/home/rsmerigl/codes/cleaned_codes')
from snputils.snputils.genobj.snp.snpobj import SNPObject
from snputils.snputils.io.snp.write.vcf import VCFWriter
from snputils.snputils.io.ancestry.local.read.mspReader import MSPReader


def vcf_creation(ancestry_map, msp_folder, output_folder):

    msp_file_list = os.listdir(msp_folder)
    print(msp_file_list)

    for msp_file_name in msp_file_list:

        # read msp file 
        msp_file = os.path.join(msp_folder, msp_file_name)

        # Create an instance of MSPReader
        msp_reader = MSPReader(filename=msp_file)

        # Read the msp file and create a LocalAncestryObject
        print(f"Reading {msp_file}")
        lai_object = msp_reader.read()

        #processing a list of position to have both the sdtart and the end position of the windows.
        pos_list = [f"{val1}_{val2}" for val1, val2 in lai_object.physical_pos]

        
        #iterate for each ancestry
        for ancestry in range(len(ancestry_map)):

            # modify the values to simulate a SNP file 
            match = (lai_object.lai == ancestry).astype(int)
            match = match.reshape(len(lai_object.lai),int(len(lai_object.lai[0])/2), 2 )

            # create vcf object
            calldata_gt = match
            samples = np.array(lai_object.sample_IDs)
            variants_chrom = lai_object.chromosome
            variants_list = [str(i+1) for i in range(len(lai_object.window_size))]
            variants_id = np.array(variants_list)
            variants_ref = np.full(calldata_gt.shape[0], 'A', dtype='U5')
            variants_alt = np.full((calldata_gt.shape[0], 3), ['T', 'G', 'C'], dtype='U1')
            qual = np.full(variants_id.shape, 100)

            # create the SNPObject
            variant_data_obj = SNPObject(
                calldata_gt=calldata_gt,
                samples=samples,
                variants_chrom=variants_chrom,
                variants_id=variants_id,
                variants_ref = variants_ref,
                variants_alt = variants_alt,
                variants_pos = pos_list,
                variants_qual = qual
            )

            # create a filename and select the output folder
            vcf_filename = f'chr{np.unique(variants_chrom)[0]}_{ancestry_map[str(ancestry)]}.vcf'
            output_file = os.path.join(output_folder, vcf_filename)

            # write the vcf file 
            vcf_writer = VCFWriter(variant_data_obj, output_file)
            vcf_writer.write()

    # return the output folder in case of necessity
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