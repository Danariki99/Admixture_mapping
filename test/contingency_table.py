import pandas as pd 
import os
import numpy as np
import pickle

beta_path = "/private/groups/ioannidislab/smeriglio/trial/betas/output.cancer1060.glm.logistic.hybrid"
pheno_path = "/private/groups/ioannidislab/smeriglio/phe_files/cancer1060.phe"
pickle_path = "/private/groups/ioannidislab/smeriglio/pickle_files/lai_object_chr_*.pkl"
keep_path = "/private/groups/ioannidislab/smeriglio/keep_file/keep_file.txt"
 

beta = pd.read_csv(beta_path, sep = "\t")

absolute_beta = beta['BETA'].abs()
index_max = absolute_beta.idxmax()
chrom = beta.loc[index_max, '#CHROM']
index_max = beta.loc[index_max, 'ID']

pheno = pd.read_csv(pheno_path, sep = " ")
pheno.columns = ['#IID', 'cancer1060']
keep = pd.read_csv(keep_path, sep = "\t")
pheno = pd.merge(pheno, keep, on = '#IID')
pickle_path = pickle_path.replace('*', str(chrom))
with open(pickle_path, 'rb') as f:
    lai_object = pickle.load(f)
lai_object.lai = lai_object.lai.reshape(len(lai_object.lai),int(len(lai_object.lai[0])/2), 2 )

df_table = pheno.copy()
df_table['#IID'] = df_table['#IID'].astype(str)
ancestry = 3

window = lai_object.lai[index_max, :, :]
window_match = (window == ancestry).astype(int)

# Crea un set con i valori di df_table['#IID']
iid_set = set(df_table['#IID'].values)

# Inizializza i contatori
counter_yes_1 = 0
counter_yes_05 = 0
counter_yes_0 = 0
counter_no_1 = 0
counter_no_05 = 0
counter_no_0 = 0

for sample in lai_object.sample_IDs:
    if sample in iid_set:
        index = df_table.index[df_table['#IID'] == sample][0]
        mean_match = window_match[index].mean()
        
        if df_table.loc[index, 'cancer1060'] == 2:
            if mean_match == 1:
                counter_yes_1 += 1
            elif mean_match == 0.5:
                counter_yes_05 += 1
            elif mean_match == 0:
                counter_yes_0 += 1
        else:
            if mean_match == 1:
                counter_no_1 += 1
            elif mean_match == 0.5:
                counter_no_05 += 1
            elif mean_match == 0:
                counter_no_0 += 1

contingency_table = pd.DataFrame({
    '1/1': [counter_yes_1, counter_no_1],
    '1/0 or 0/1': [counter_yes_05, counter_no_05],
    '0/0': [counter_yes_0, counter_no_0]
}, index=['sì', 'no'])

print(contingency_table)

        