import random
import subprocess

def get_sample_list_from_vcf(vcf_file):
    sample_list = []
    try:
        # Esegui il comando bcftools per ottenere la lista dei campioni dal file VCF
        result = subprocess.run(['bcftools', 'query', '-l', vcf_file], capture_output=True, text=True, check=True)
        # Dividi il risultato in righe e rimuovi eventuali spazi vuoti
        sample_list = result.stdout.strip().split('\n')
    except subprocess.CalledProcessError as e:
        print(f"Errore durante l'estrazione della lista dei campioni: {e.stderr}")
    return sample_list


def generate_fake_phenotype_file(sample_list, num_phenotypes, filename):
    with open(filename, 'w') as file:
        # Write header
        file.write('IID ')
        for i in range(1, num_phenotypes + 1):
            file.write(f'PHENO{i} ')
        file.write('\n')
        
        # Write data
        for i, sample in enumerate(sample_list, start=1):
            file.write(f'{sample} ')
            for _ in range(num_phenotypes):
                # Generate random phenotype values (0 or 1)
                phenotype_value = random.randint(1, 2)
                file.write(f'{phenotype_value} ')
            file.write('\n')


# Usage example
vcf_file = '/private/groups/ioannidislab/smeriglio/vcf_files/chr22/ancestry_EUR.vcf'  # Sostituisci con il percorso del tuo file VCF
sample_list = get_sample_list_from_vcf(vcf_file)
print("Lista dei campioni nel file VCF:")
print(len(sample_list))
num_phenotypes = 4  # Modifica a seconda del numero di fenotipi
filename = '/private/groups/ioannidislab/smeriglio/trial/fake_pheno.txt'  # Modifica con il percorso desiderato
generate_fake_phenotype_file(sample_list, num_phenotypes, filename)
print(f"Fake phenotype file '{filename}' created successfully.")
