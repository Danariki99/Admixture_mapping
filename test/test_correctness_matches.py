import pickle
import numpy as np
import pandas as pd
import os
from snputils.genobj.snp.snpobj import SNPObject
from snputils.io.snp.write.vcf import VCFWriter

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


for chrom_n in range(num_chrom):
    chrom_n+=1

    print(f'starting chromosome {chrom_n}')
    
    # read object
    pickle_file = "/private/groups/ioannidislab/smeriglio/pickle_files/lai_object_chr_" + str(chrom_n) + ".pkl"

    chrom = "chr" + str(chrom_n)

    chrom_folder = os.path.join("/private/groups/ioannidislab/smeriglio/vcf_files", chrom)

    if not os.path.exists(chrom_folder):
        os.mkdir(chrom_folder)


    with open(pickle_file, 'rb') as f:
        lai_object = pickle.load(f)

    #print(lai_object.physical_pos)

    pos_list = [f"{val1}_{val2}" for val1, val2 in lai_object.physical_pos]
    #print(pos_list)

    
    #iterate for each ancestry
    matches = []
    for ancestry in range(len(ancestry_map)):

        # modify the values to simulate a SNP file 
        match = (lai_object.lai == ancestry).astype(int)
        match = match.reshape(len(lai_object.lai),int(len(lai_object.lai[0])/2), 2 )

        matches.append(match)
    
    print(f'finished production of matches for chrom: {chrom_n}')

    matches_array = np.array(matches)

    sum_matches = np.sum(matches_array, axis=0)

    match_test = np.ones(sum_matches.shape)
    
    if not np.array_equal(match_test, sum_matches):
        print(f"Error in chromosome {chrom_n}")
    