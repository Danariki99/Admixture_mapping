import os
import subprocess

filename = '/private/groups/ioannidislab/smeriglio/covar_file/ukb24983_GWAS_covar.phe'

# Remove the non-numeric column from the file using awk and overwrite the file
command = f"awk '{{printf \"%s\\t\", $2; printf \"%s\\t\", $5; for (i=10; i<=NF; i++) printf \"%s\\t\", $i; print \"\"}}' {filename} > temp_file && mv temp_file {filename}"



subprocess.run(command, shell=True)