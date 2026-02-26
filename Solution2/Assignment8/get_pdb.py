#!/usr/bin/python3
import argparse
import requests


#3d struktur von Molekülen, fasta buchstaben aber pdb wo (x y z) atom im raum
#RCSB website portal zum archiv
#Protein Struktur bestimmt Aufgabe
#currently 4-character alphanumeric

parser = argparse.ArgumentParser()
parser.add_argument('--id', required=True)
parser.add_argument('--output',type=argparse.FileType('w'), required=True)
#ooptional ohne argument
parser.add_argument('--fasta', action='store_true')
args = parser.parse_args()

if len(args.id) != 4:
    print("Ungültige ID")
else:
    pdb_id = args.id.upper()

output = args.output

#rechtsklick link speichern
if args.fasta:
    url = f"https://files.rcsb.org/download/{pdb_id}.fasta"
else:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"


try:
    response = requests.get(url)
    # ausgabe auf konsole automatishx bei -
    output.write(response.text)

except Exception as e:
    print(e)








