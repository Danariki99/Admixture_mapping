import subprocess

n_chrom = 22

for chrom in range(1, n_chrom + 1):
    print('Processing chromosome', chrom)
    subprocess.run(['sbatch', 'sbatch.sh', str(chrom)])