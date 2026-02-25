#!/usr/bin/python3

#Swiss-Prot AC-Nummer ist ID für UniProtKB/Swiss_Prot-Datenbank, 6-10 Zeichen
#AC Nummer für jede Proteinsequenz, eine primäre AC und optional sekundäre ACs
import argparse
import requests

#Input: AC-Nummer mit --ac, muss vorhanden sein
parser = argparse.ArgumentParser()
parser.add_argument('--ac', required=True)
args = parser.parse_args()

#Nummer darf nur Buchstaben und Zahlen enthalten, args.ac wegen --ac
if args.ac.isalnum():
    ac_number = args.ac
else:
    print("Invalid AC")

#Ausgabe: passende Sequenz im FASTA Format
#Sequenz auf der Website finden: AC eingeben, zu Sequenz scrollen und runterladen
#-> daher auch passende url gefunden

response = requests.get(f"https://rest.uniprot.org/uniprotkb/{ac_number}.fasta")
print(response.text)

