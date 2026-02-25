#!/usr/bin/python3
import argparse
import urllib.request
import re

#Genome Report für Bakterien von NCBI runterladen

parser = argparse.ArgumentParser()
#nargs='+' für mehrere Eingaben nach --organims
parser.add_argument('--organism', nargs='+', required=True)
args = parser.parse_args()
#ALle Eingaben speichern
patterns=[re.compile(p, re.IGNORECASE) for p in args.organism]


#Datei von angegebenem Ort nehmen, ftp ging nicht? also https
#-> lädt Genome Report runter
url = "https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"

response = urllib.request.urlopen(url)

#Erkenntnisse aud Datei: TSV, Name an 1. Stelle, Länge 7.Stelle, Status 16.
#für alle Ausdrücke nach Genomen + Organismus Namen suchen & Status anschauen ob Complete Genome

for line_bytes in response:
    # Bytes zu String umwandeln
    line = line_bytes.decode('utf-8', errors='ignore')

    if not line or line.startswith('#'):
        continue

    column = line.split('\t')
    matched_name = column[0].strip()
    length = column[6].strip()
    status = column[15].strip()

    #Nach Organism Name in Liste suchen
    if status == "Complete Genome":
        for pattern in patterns:
            if pattern.search(matched_name):
                print(f"{matched_name}\t{length}")
                break
