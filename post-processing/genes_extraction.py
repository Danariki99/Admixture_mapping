import requests
import pandas as pd
import json

def fetch_genes_ucsc(chromosome, start, end):
    url = f"http://api.genome.ucsc.edu/getData/track?genome=hg38;track=refGene;chrom=chr{chromosome};start={start};end={end}"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "refGene" in data:
            genes = data["refGene"]
            return genes
        else:
            print("No genes found in the specified region.")
            return []
    else:
        print(f"Error fetching data from UCSC: {response.status_code}")
        return []

if __name__ == '__main__':
    try:
        data_file = pd.read_table('/private/groups/ioannidislab/smeriglio/output/no_BMI/output_ancestry_AFR/P_info_AFR.tsv', sep='\t')
        print("Data file loaded successfully.")
    except Exception as e:
        print("Error loading data file:", e)
        exit()

    gene_dict = {}
    gene_df = pd.DataFrame(columns=['gene_name'])
    gene_names = []

    for index, row in data_file.iterrows():
        try:
            chrom = int(row['#CHROM'])  # Assicura che il cromosoma sia un numero intero
            start = int(row['POS'])
            end = int(row['end_POS'])
            print(f"Fetching genes for: chr{chrom}, start: {start}, end: {end}")
            genes = fetch_genes_ucsc(chrom, start, end)

            if genes:
                gene_info = []

                for gene in genes:
                    gene_name = gene.get('name2', 'N/A')
                    gene_info.append({
                        'gene_name': gene_name,
                        'transcript_id': gene.get('name', 'N/A'),
                        'chromosome': gene.get('chrom', 'N/A'),
                        'tx_start': gene.get('txStart', 'N/A'),
                        'tx_end': gene.get('txEnd', 'N/A')
                        
                    })
                    gene_names.append(gene_name)

                key = f"{chrom}_{start}_{end}"  # Convert tuple to string
                gene_dict[key] = gene_info

                

                for info in gene_info:
                    print(f"Gene Name: {info['gene_name']}")
                    print(f"Transcript ID: {info['transcript_id']}")
                    print(f"Chromosome: {info['chromosome']}")
                    print(f"Start: {info['tx_start']}")
                    print(f"End: {info['tx_end']}")
                    print("---")
            else:
                print("No genes found in the specified region.")
                key = f"{chrom}_{start}_{end}"  # Convert tuple to string
                gene_dict[key] = []

        except Exception as e:
            print(f"Error processing row: {e}")

    # Stampa il dizionario finale
    print("Gene Dictionary:")
    for key, value in gene_dict.items():
        print(f"{key}: {value}")

    # save the gene dictionary to a file
    with open('/private/groups/ioannidislab/smeriglio/genes/gene_dict.json', 'w') as f:
        json.dump(gene_dict, f, indent=4)

    gene_names = list(set(gene_names))
    gene_df['gene_name'] = gene_names
    
    gene_df.to_csv('/private/groups/ioannidislab/smeriglio/genes/gene_names.tsv', sep='\t', index=False)
