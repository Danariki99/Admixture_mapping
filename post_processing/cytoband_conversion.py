import requests

def fetch_cytoband(chromosome, start, end, genome="hg38"):
    """
    Contatta il database UCSC e restituisce la banda cromosomica corrispondente
    a un intervallo genomico specificato.
    
    Args:
        chromosome (str): Nome del cromosoma (es. "chr6").
        start (int): Posizione iniziale nell'intervallo.
        end (int): Posizione finale nell'intervallo.
        genome (str): Versione del genoma (es. "hg38", "hg19", "hg18").
    
    Returns:
        str: Nome della banda cromosomica corrispondente o "Unknown" se non trovata.
    """
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
    
    # Effettua la richiesta a UCSC
    response = requests.post(url, data=params)
    if response.status_code != 200:
        raise Exception(f"Errore durante la connessione a UCSC: {response.status_code}")
    
    # Analizza la risposta
    response_text = response.text.strip()
    if not response_text:
        return "Unknown"
    
    # Converte l'output in righe e seleziona la banda più rilevante
    lines = response_text.split("\n")
    fields = lines[1].split("\t")  # Splitta la riga in una lista
    name = f"{fields[0].replace('chr', '')}{fields[3]}"
    return name


# Esempio di utilizzo
if __name__ == "__main__":
    chromosome = "chr6"
    start = 31284326
    end = 31461089
    genome = "hg38"
    
    band = fetch_cytoband(chromosome, start, end, genome)
    print(f"La banda cromosomica per {chromosome}:{start}-{end} è {band}")
