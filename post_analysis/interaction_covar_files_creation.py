import os
import pandas as pd

covar_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_covar_files'

interaction_folder = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/PCA_files/PCA_covar_files_interaction'

files_list = os.listdir(covar_folder)

for file in [files_list[0]]:

    df = pd.read_csv(os.path.join(covar_folder, file), sep='\t')

# thinking about this, the only way to do this thing is to take the significant SNPs according to the analysis without the interaction term and
# compute an interaction term only for those SNPS since we will have to run the association SNP by SNP
# in fact the interaction term wull be different for each of the SNPs so also the final analisys will be different.