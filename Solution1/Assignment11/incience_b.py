#!/usr/bin/python3
import argparse
import sys

#Gegebene Seq
sequences = ["CTAG", "CG", "AACCCTGTC", "ATG"]
#Genome runtergeladen
genome_path = "C:/Users/gross/Documents/uni/ProPra/Mycobacterium_tuberculosis_CDC1551_genomic.fna"

genome_parts = []
with open(genome_path, "r") as file:
    for line in file:
        if line.startswith(">"):
            continue
        #strip() wichtig! sonst anderes Ergebnis
        genome_parts.append(line.strip())

genome_sequence = "".join(genome_parts)
gen_length = len(genome_sequence)

# tatsächliche relative Häufigkeit
rh = {
    'A': genome_sequence.count('A') / gen_length,
    'C': genome_sequence.count('C') / gen_length,
    'G': genome_sequence.count('G') / gen_length,
    'T': genome_sequence.count('T') / gen_length
}

# Erste Zeile der Ausgabe
print(f"A: {rh['A'] * 100:.2f}% C: {rh['C'] * 100:.2f}% G: {rh['G'] * 100:.2f}% T: {rh['T'] * 100:.2f}%")

for seq in sequences:
    k = len(seq)

    exp = gen_length * (0.25 ** k)
    #geg hier p = 1
    p_rh = 1.0

    #Anzahl erwartete Indizien
    for base in seq:
        #get(base, 0) für alle Zeichen, die nicht A,C,G,T sind
        p_rh *= rh.get(base, 0)
        exp_p_rh = gen_length * p_rh

    #Tatsächliche Inzidenzen
    inc_count = 0
    for i in range(gen_length - k + 1):
        # aktuelles fenster
        s = genome_sequence[i:i + k]
        if s == seq:
            inc_count += 1

    #Fold-Change
    fold_change = inc_count / exp_p_rh

    print(f"{seq}:\t{exp:.0f}\t{exp_p_rh:.0f}\t{inc_count}\t{fold_change:.1f}\t")
