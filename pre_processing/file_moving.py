import os
import subprocess
import pandas as pd
import shutil

def excel_reading(filename, sheet_name):
    try:
        df = pd.read_excel(filename, sheet_name=sheet_name)
        return df
    except Exception as e:
        print("Errore durante la lettura del file Excel:", e)
        return None




if __name__ == '__main__':

    ukbb_dir = '/private/groups/ioannidislab/ukbb24983'

    os.chdir(ukbb_dir)
    
    excel_file = '/private/groups/ioannidislab/smeriglio/excel_file/ukbb_v1.xlsx'

    sheet_name = 'first_batch'

    excel_df = excel_reading(excel_file, sheet_name)

    phenotype_IDs = excel_df['ID'].tolist()

    suffix = '.phe'

    pheno_files_list = []

    for pheno in phenotype_IDs:
        pheno_files_list.append(pheno+suffix)

    output_folder = '/private/groups/ioannidislab/smeriglio/phe_files'

    pheno_folder = '/private/groups/ioannidislab/ukbb24983/phenotypedata/extras/highconfidenceqc/v2_2020/phe'

    for pheno in pheno_files_list:

        if pheno[:2] == 'HC':

            pheno_folder = '/private/groups/ioannidislab/ukbb24983/phenotypedata/extras/highconfidenceqc/v2_2020/phe'

            pheno_file = os.path.join(pheno_folder, pheno)

            shutil.copy(pheno_file, output_folder)

        elif pheno[:2] == 'FH':

            pheno_folder = '/private/groups/ioannidislab/ukbb24983/phenotypedata/extras/family_history/phe/'

            pheno_file = os.path.join(pheno_folder, pheno)

            shutil.copy(pheno_file, output_folder)
        
        elif pheno[:2] == 'ca':

            pheno_folder = '/private/groups/ioannidislab/ukbb24983/phenotypedata/extras/cancer/v2_2020/phe'

            pheno_file = os.path.join(pheno_folder, pheno)

            shutil.copy(pheno_file, output_folder)
        else:
            print(pheno)

    file_list = os.listdir(output_folder)

    for name in file_list:

        pheno_name = name.split('.')[0]
    
        # Construct the full file path
        filename = os.path.join(output_folder, name)

        # Remove the first column from the file using awk and overwrite the file
        command = f"awk '{{sub(/^\\S+\\s*/, \"\", $0)}} 1' {filename} > temp_file && mv temp_file {filename}"
        subprocess.run(command, shell=True)

        # Construct the awk command to add the new header line
        awk_command = f"awk -v header='#IID {pheno_name}' 'BEGIN {{print header}} {{print}}' {filename} > temp_file && mv temp_file {filename}"

        # Run the awk command to add the header directly to the file
        subprocess.run(awk_command, shell=True)