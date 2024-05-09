import subprocess

# Define the file paths
input_file_general = '/private/groups/ioannidislab/smeriglio/merged_vcfs/ancestry_*.vcf'
output_file_general = '/private/groups/ioannidislab/smeriglio/covar_file/*_freq.txt'

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

    print('processing', ancestry)

    input_file = input_file_general.replace('*', ancestry)
    output_file = output_file_general.replace('*', ancestry)

    # Define the awk command
    awk_command = f'''awk 'BEGIN {{OFS=FS="\\t"}} $1 ~ /^##/ {{next}} $1 ~ /^#/ {{for(i=10; i<=NF; i++) id[i]=$i; next}} {{for(i=10; i<=NF; i++) {{split($i,a,":"); if(a[1]=="0/0") count[i]+=0; else if(a[1]=="1/0" || a[1]=="0/1") count[i]+=0.5; else if(a[1]=="1/1") count[i]+=1}}}} END {{print "IID", "{ancestry}"; for(i=10; i<=NF; i++) print id[i], count[i]}}' '''
    # Combine the command and the file paths
    full_command = f'{awk_command} {input_file} > {output_file}'

    # Run the command
    subprocess.run(full_command, shell=True, check=True)