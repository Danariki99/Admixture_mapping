import pandas as pd
import matplotlib.pyplot as plt
import os
from qmplot import manhattanplot

folder = '/private/groups/ioannidislab/smeriglio/test-omit-ref/'
pheno = 'HC132'
base_normal_filename = 'output.*.glm.logistic.hybrid'
base_adjusted_filename = 'output.*.glm.logistic.hybrid.adjusted'

pheno_folder = folder

current_normal_filename = os.path.join(pheno_folder,base_normal_filename).replace('*',pheno)
df1 = pd.read_table(current_normal_filename, sep="\t")
df1 = df1.drop('P', axis = 1)

current_adjusted_filename = os.path.join(pheno_folder,base_adjusted_filename).replace('*',pheno)
df2 = pd.read_table(current_adjusted_filename, sep="\t")

merged_df = pd.merge(df1, df2, on=['#CHROM', 'ID'])
merged_df=merged_df.rename(columns={'UNADJ':'P'})

filename = os.path.join(pheno_folder, 'manhattan_plot.png')


# generate manhattan plot and set an output file.
ax = manhattanplot(data=merged_df, xticklabel_kws={"rotation": "vertical"})
plt.savefig(filename, bbox_inches='tight', pad_inches=0.1)
plt.close()

