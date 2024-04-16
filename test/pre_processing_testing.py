import pickle
import numpy as np
import pandas as pd
import os
from snputils.genobj.snp.snpobj import SNPObject
from snputils.io.snp.write.vcf import VCFWriter

chrom_n = 1
print(chrom_n)

# read object
pickle_file = "/private/groups/ioannidislab/smeriglio/pickle_files/lai_object_chr_" + str(chrom_n) + ".pkl"

chrom = "chr" + str(chrom_n)

chrom_folder = os.path.join("/private/groups/ioannidislab/smeriglio/trial", chrom)

if not os.path.exists(chrom_folder):
    os.mkdir(chrom_folder)


with open(pickle_file, 'rb') as f:
    lai_object = pickle.load(f)

#print(lai_object.physical_pos)

pos_list = [val2 for val1, val2 in lai_object.physical_pos]
#print(pos_list)

ancestry = 3

# modify the values to simulate a SNP file 
match = (lai_object.lai == ancestry).astype(int)
match = match.reshape(len(lai_object.lai),int(len(lai_object.lai[0])/2), 2 )

# create vcf object
calldata_gt = match
samples = np.array(lai_object.sample_IDs)
variants_chrom = np.full(calldata_gt.shape[0], chrom, dtype='U5')
variants_list = [str(i+1) for i in range(len(lai_object.window_size))]
variants_id = np.array(variants_list)
variants_ref = np.full(calldata_gt.shape[0], 'A', dtype='U5')
variants_alt = np.full((calldata_gt.shape[0], 3), ['T', 'G', 'C'], dtype='U1')

variant_data_obj = SNPObject(
    calldata_gt=calldata_gt,
    samples=samples,
    variants_chrom=variants_chrom,
    variants_id=variants_id,
    variants_ref = variants_ref,
    variants_alt = variants_alt,
    variants_pos = pos_list
)

print(variant_data_obj.calldata_gt.shape)
print(variant_data_obj.calldata_gt[:,8325,1])

# save ancestry data as vcf
filename = 'ancestry_' + 'EUR' + '.vcf'
output_file = os.path.join(chrom_folder, filename) 
vcf_writer = VCFWriter(variant_data_obj, output_file)
print(output_file)

vcf_writer.write()