import pandas as pd
import matplotlib.pyplot as plt

# Leggi il file
df = pd.read_csv('/private/groups/ioannidislab/smeriglio/trial/betas/output.cancer1060.glm.logistic.hybrid', sep='\t')

# Crea il barplot
plt.figure(figsize=(10,6))
plt.bar(range(len(df['BETA'].tolist())), df['BETA'])
plt.title('Betas')
plt.xlabel('ID')
plt.ylabel('BETA')
plt.savefig('/private/groups/ioannidislab/smeriglio/trial/betas/barplot.png')