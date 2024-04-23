import pandas as pd

filename = '/private/groups/ioannidislab/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
filename_2 = '/private/groups/ioannidislab/smeriglio/trial/ukb24983_GWAS_covar.phe'

# Read the file into a pandas DataFrame
print('reading file')
df = pd.read_csv(filename, sep='\t', header=None, low_memory=False)
print('finish reading file\nstarting pre processing')

# Remove the first column and keep all the others
df = df.iloc[:, [1, 4, 9, 14,15,16,17,18,19,20,21,22,23]]

df = df.dropna()

# Write the DataFrame back to the file
df.to_csv(filename_2, sep='\t', header=False, index=False)
print(f'file saved as {filename_2}')