#!/usr/bin/python3

#CRISPR Veränderung DNA basierend auf bakteriellem Immunsystem
#Gene ausschalten/korrigieren/neue Sequenzen einfügen
#guide-RNA führt Cas9-Protein (Schere) zur spezifischen DNA-Seq die durchtrennt wird
#Natürliches Abwehrsystem gg Viren von Bakterien
#CRISPR an sich = DNA-Abschnitte im Bakteriengenom      also allg Genschere

import argparse

#FASTA bekommen
parser = argparse.ArgumentParser()
parser.add_argument('--fasta', required=True)
args = parser.parse_args()

header = ""
genome_parts = []
crispr_length = 23
#schon ORF von AUfgabe 10
with open (args.fasta, "r") as file:
    for line in file:
        line = line.strip()

        if line.startswith(">"):
            if header != "":
                genome_sequence = "".join(genome_parts)
                for i in range(len(genome_sequence) - 22):
                    # sucht alle crispr stellen
                    crispr = genome_sequence[i:i + crispr_length]
                    if crispr[21:23] == "GG":
                        print(f"{header}\t{i + 1}")
                        print(crispr)
            #header setzen wenn noch keiner oder vorherigen orf gelesen
            header = line
            #leeren um es wieder frei für nächsten ORF
            genome_parts = []
            continue; #sonst landet header auch in der Sequenz
        #wenn nicht mit > appenden
        genome_parts.append(line)
    #letzter ORF
    if header != "":
        genome_sequence = "".join(genome_parts)
        for i in range(len(genome_sequence) - 22):
            crispr = genome_sequence[i:i + crispr_length]
            if crispr[21:23] == "GG":
                print(f"{header}\t{i + 1}")
                print(crispr)
