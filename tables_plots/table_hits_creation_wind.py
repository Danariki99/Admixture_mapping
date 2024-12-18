# to add as part of the supplementary information, we should include a table (excel)
# with mean age + sd in cases and controls per phenotypes, as well %women.
# I can help if you need!

import os
import pandas as pd
import openpyxl
import requests
import numpy as np

import requests

def fetch_cytoband(chromosome, start, end, genome="hg19"):
    url = "https://genome.ucsc.edu/cgi-bin/hgTables"
    params = {
        "db": genome,
        "hgta_group": "allTracks",
        "hgta_track": "cytoBand",
        "hgta_table": "cytoBand",
        "hgta_regionType": "range",
        "position": f"{chromosome}:{start}-{end}",
        "hgta_outputType": "primaryTable",
        "boolshad.sendToGalaxy": "0",
        "boolshad.sendToGreat": "0",
        "hgta_doTopSubmit": "get output",
    }

    try:
        response = requests.post(url, data=params, timeout=30)
        response.raise_for_status()  # Gestisce errori HTTP
    except requests.exceptions.ReadTimeout:
        print(f"Timeout: il server non ha risposto entro il tempo previsto.")
        return "Timeout"
    except requests.exceptions.RequestException as e:
        print(f"Errore nella richiesta: {e}")
        return "Errore"
    
    # Analizza la risposta
    response_text = response.text.strip()
    if not response_text:
        return "Unknown"
    
    lines = response_text.split("\n")
    fields = lines[1].split("\t")
    return f"{fields[0].replace('chr', '')}{fields[3]}"

# Paths to important files and directories
phe_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/phe_files/ukbb'
covar_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/covar_file/ukbb/ukb24983_GWAS_covar.phe'
keep_file = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/keep_file/ukbb/keep_file.txt'
significant_path = '/private/groups/ioannidislab/smeriglio/out_cleaned_codes/vcf_files_windows/ukbb/output_more_wind_3/'

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

# Initialize an empty DataFrame to store results
out1 = pd.DataFrame(columns=[
    'Phenotype',
    'CytoBand',
    'ancestry tested',
    'OR (CI = 95%)',
    'p value',
    'START',
    'END',
    'chr'
])
covar = covar.rename(columns={'IID': '#IID'})

df_first_batch = pd.read_excel("ukbb_v1.xlsx", sheet_name="first_batch")

# Iterate over the files in the result path
succesfully_pheno_anc_list = os.listdir(significant_path)

for fold in succesfully_pheno_anc_list:
    ancestry = fold.split('/')[-1].split('_')[0]
    pheno = fold.split('/')[-1].split('_')[1]
    chrom = fold.split('/')[-1].split('_')[2]

    pheno_name = df_first_batch[df_first_batch['ID'] == pheno]['ID2']

    result_path = f'/private/groups/ioannidislab/smeriglio/out_cleaned_codes/output/ukbb/output_ancestry_{ancestry}/{pheno}/output.{pheno}.glm.logistic.hybrid'

    chr = int(chrom.replace('chr', ''))

    df = pd.read_table(result_path, sep='\t')

    filtered_df = df[df['#CHROM'] == chr]
    index_min = filtered_df['P'].idxmin()
    row_with_min_p = filtered_df.loc[index_min]

    start = row_with_min_p['POS']
    end = filtered_df.loc[index_min + 1]['POS']

    print('fetching cytoband')
    cytoband = fetch_cytoband(chrom, start, end)

    p = row_with_min_p['P']
    oddr = row_with_min_p['OR']

    out1 = pd.concat([
        out1,
        pd.DataFrame([{
            'Phenotype': pheno_name.iloc[0] if not pheno_name.empty else 'Unknown',  # Nome del fenotipo
            'CytoBand': cytoband,  
            'ancestry tested': ancestry,  
            'OR (CI = 95%)': oddr, 
            'p value': p,
            'START': start,
            'END': end,
            'chr': chr
        }])

    ], ignore_index=True)

out1.to_excel('table_hits_wind.xlsx', index=False)


    
    
    


    
