# to add as part of the supplementary information, we should include a table (excel)
# with mean age + sd in cases and controls per phenotypes, as well %women.
# I can help if you need!

import os
import pandas as pd
import openpyxl
import numpy as np

# Paths to important files and directories
phe_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb'
covar_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar.phe'
keep_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt'

# Read the covariate file
try:
    covar = pd.read_csv(covar_file, sep='\t')
except Exception as e:
    raise FileNotFoundError(f"Error reading the covariate file: {e}")

# Read the keep file
try:
    keep = pd.read_csv(keep_file, sep='\t')
except Exception as e:
    raise FileNotFoundError(f"Error reading the keep file: {e}")

# List of phenotype files in the directory
phe_files = os.listdir(phe_path)

# Initialize an empty DataFrame to store results
out = pd.DataFrame(columns=[
    'Phenotype', 
    'Mean_Age_Cases', 'SD_Age_Cases',
    'Mean_Age_Controls', 'SD_Age_Controls',
    '%Women_Cases', '%Men_Cases',
    '%Women_Controls', '%Men_Controls'
])
covar = covar.rename(columns={'IID': '#IID'})

# Loop through phenotype files
for phe_file in phe_files:
    try:
        # Read the phenotype file
        phe = pd.read_csv(f'{phe_path}/{phe_file}', sep=' ')
    except Exception as e:
        print(f"Error reading the file {phe_file}: {e}")
        continue
    phe_name = phe_file.replace('.phe', '')

    if pd.isna(phe[phe_name].iloc[0]):
        # Split the values in the '#IID' column by '\t' and convert to integers
        phe[['#IID', phe_name]] = phe['#IID'].str.split('\t', expand=True)
        phe['#IID'] = phe['#IID'].astype(int)
        phe[phe_name] = phe[phe_name].astype(int)

    
    # Merge phenotype data with covariate and keep files
    phe = phe.merge(covar, on='#IID', how='inner')
    phe = phe.merge(keep, on='#IID', how='inner')
    
    # Filter cases and controls based on the phenotype column
    phe_cases = phe[phe[phe_name] == 2]
    phe_controls = phe[phe[phe_name] == 1]

    # Calculate statistics for cases
    mean_age_cases = phe_cases['age'].mean()
    sd_age_cases = phe_cases['age'].std()
    woman_percentage_cases = (phe_cases[phe_cases['sex'] == 1].shape[0] / phe_cases.shape[0]) * 100 if phe_cases.shape[0] > 0 else 0
    man_percentage_cases = (phe_cases[phe_cases['sex'] == 0].shape[0] / phe_cases.shape[0]) * 100 if phe_cases.shape[0] > 0 else 0

    # Calculate statistics for controls
    mean_age_controls = phe_controls['age'].mean()
    sd_age_controls = phe_controls['age'].std()
    woman_percentage_controls = (phe_controls[phe_controls['sex'] == 1].shape[0] / phe_controls.shape[0]) * 100 if phe_controls.shape[0] > 0 else 0
    man_percentage_controls = (phe_controls[phe_controls['sex'] == 0].shape[0] / phe_controls.shape[0]) * 100 if phe_controls.shape[0] > 0 else 0

    # Append results for the current phenotype to the output DataFrame
    # Create a new row as a DataFrame
    new_row = pd.DataFrame([{
        'Phenotype': phe_name,  # Use the file name as the phenotype identifier
        'Mean_Age_Cases': mean_age_cases,
        'SD_Age_Cases': sd_age_cases,
        'Mean_Age_Controls': mean_age_controls,
        'SD_Age_Controls': sd_age_controls,
        '%Women_Cases': woman_percentage_cases,
        '%Men_Cases': man_percentage_cases,
        '%Women_Controls': woman_percentage_controls,
        '%Men_Controls': man_percentage_controls
    }])

    # Concatenate the new row with the existing DataFrame
    out = pd.concat([out, new_row], ignore_index=True)

# Save the results to an Excel file
print(out)
out.to_excel('phenotype_statistics.xlsx', index=False)



    
