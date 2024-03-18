import logging
import pickle
from snputils.io.ancestry.local.read.mspReader import MSPReader

# Configure logging
logging.basicConfig(level=logging.INFO)

# Define the path to the msp file
msp_file_base = "/private/groups/ioannidislab/ukbb_analysis/rfmix/output/ukb_hap_chr*_v2_rfmix.msp.tsv"
    
# Define the path to save the pickle file
pickle_file_base = "/private/groups/ioannidislab/smeriglio/lai_object_chr_*.pkl"
num_chr = 22

for i in range(num_chr):

    chrom  = i+1

    # define current msp
    current_msp = msp_file_base.replace('*', str(chrom))

    # define current pickle file
    current_pickle_file = pickle_file_base.replace('*', str(chrom))

    # Create an instance of MSPReader
    msp_reader = MSPReader(filename=current_msp)

    # Read the msp file and create a LocalAncestryObject
    lai_object = msp_reader.read()



    # Save lai_object as pickle file
    with open(current_pickle_file, 'wb') as f:
        pickle.dump(lai_object, f)

    # Print a message indicating that the pickle file has been saved
    print(f"LocalAncestryObject saved as pickle file: {current_pickle_file}")