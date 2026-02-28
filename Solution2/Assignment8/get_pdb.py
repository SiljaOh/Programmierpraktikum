#!/usr/bin/python3
import argparse
import sys
import urllib.request
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

#passt da nur eine ID eingegeben
if len(args.id) == 4 and args.id[0].isdigit():
    pdb_id = args.id
else:
    print("Ungültige ID")
    sys.exit(1)

output = args.output

if args.fasta:
    url = f"https://files.rcsb.org/download/{pdb_id}.fasta"
else:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

try:
    data = urllib.request.urlopen(url)
    response = data.read().decode('utf-8')
    output.write(response)

except Exception as e:
    print(e)
    sys.exit(1)
