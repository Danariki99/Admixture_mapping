import pandas as pd
import matplotlib.pyplot as plt
import os
from qmplot import manhattanplot

def manhattan_plot(folder):

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

    for i in range(len(ancestry_map)):

        ancestry = ancestry_map[str(i)]

        print('starting manhattan plot for ancestry: ' + ancestry)

        current_folder = folder.replace('*', ancestry)

        phenotypes = os.listdir(current_folder)
        for pheno in phenotypes:
            base_normal_filename = 'output.*.glm.logistic.hybrid'
            base_adjusted_filename = 'output.*.glm.logistic.hybrid.adjusted'

            pheno_folder = os.path.join(current_folder,pheno)

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

manhattan_plot('/private/groups/ioannidislab/smeriglio/output/no_BMI/output_ancestry_*/')


