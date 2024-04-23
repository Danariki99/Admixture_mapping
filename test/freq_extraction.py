import os
import pandas as pd

generic_filename = "/private/groups/ioannidislab/smeriglio/freq_files/*.afreq"

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
    print(f'starting ancestry {i}')
    ancestry = ancestry_map[str(i)]
    
    freq_files = os.path.join(generic_filename.replace("*", ancestry))
    
    freq_df = pd.read_csv(freq_files, sep="\t")
    print(f'mean freq for ancestry {ancestry} = {freq_df["ALT_FREQS"].mean()}')