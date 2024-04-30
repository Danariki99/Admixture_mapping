import pandas as pd
import matplotlib.pyplot as plt

# Leggi il file
df = pd.read_csv('/private/groups/ioannidislab/smeriglio/phe_files/cancer1060.phe', sep='\t')

print(df.columns())

