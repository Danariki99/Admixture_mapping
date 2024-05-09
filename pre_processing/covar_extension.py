import pandas as pd

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

#filename = '/private/groups/ioannidislab/ukbb24983/sqc/ukb24983_GWAS_covar.phe'
filename = '/private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar.phe'

generic_file = '/private/groups/ioannidislab/smeriglio/covar_file/*_freq.txt'

covar_new_file = '/private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar_new.phe'

# Read the file into a pandas DataFrame
print('reading file')
df = pd.read_csv(filename, sep='\t', header=None, low_memory=False)

print('finish reading file\nstarting pre processing')

# Remove the first column and keep all the others
#df = df.iloc[:, [1, 4, 9, 10, 14,15,16,17,18,19,20,21,22,23]]
print(df)


#df = df.dropna()
df.columns = df.iloc[0]
df = df[1:]
df['IID'] = df['IID'].astype(int)
print(df)

# Write the DataFrame back to the file
for i in range(len(ancestry_map)):

    ancestry = ancestry_map[str(i)]
    file = generic_file.replace('*', ancestry)
    print('reading', file)

    new_covar = pd.read_csv(file, sep = "\t")
    new_covar['IID'] = new_covar['IID'].astype(int)
    new_covar[ancestry] = new_covar[ancestry].astype(float)
    print(new_covar)

    df = pd.merge(df, new_covar, on = 'IID')

#covars_new_df = covars_new_df.dropna()
covars_new_df = df.fillna('NA')

print(covars_new_df)

covars_new_df.to_csv(covar_new_file, sep = "\t", index = False)
