import subprocess
import pandas as pd
from io import StringIO

# Define the command
command = "awk -F'\t' '!/^##/ {print $1\"\t\"$2}' /private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_EUR.vcf"

# Run the command
process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
output, error = process.communicate()

data = StringIO(output.decode())
df = pd.read_csv(data, sep="\t")

# Split the 'POS' column into 'start_pos' and 'end_pos'
df[['POS', 'end_POS']] = df['POS'].str.split('_', expand=True)

df['#CHROM'] = df['#CHROM'].str.replace('chr', '').astype(int)

df.to_csv('/private/groups/ioannidislab/smeriglio/windows_pos/positions.csv', sep='\t', index=False)