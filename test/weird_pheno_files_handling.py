import os
import pandas as pd
import matplotlib.pyplot as plt

pheno_folder = '/private/groups/ioannidislab/smeriglio/phe_files/'
output_folder = '/private/groups/ioannidislab/smeriglio/phenotype_histograms'

keep_files = '/private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt'

keep = pd.read_csv(keep_files)

keep['#IID'] = keep['#IID'].astype(int)

pheno_file = 'FH1001.phe'

pheno = pheno_file.replace('.phe', '')

current_pheno_path = os.path.join(pheno_folder, pheno_file)

# Read the DataFrame without the first line (header)
df = pd.read_csv(current_pheno_path, sep='\t', header=None, skiprows=1)

# Manually set the column names
df.columns = ['#IID', pheno]

df['#IID'] = df['#IID'].astype(int)

print(df)