import subprocess
import os

wind_folder = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/FUMA/ukbb/wind'

list_of_files = os.listdir(wind_folder)

for wind_filename in list_of_files:
    print('Processing file', wind_filename)
    subprocess.run(['sbatch', 'sbatch.sh', str(wind_filename)])