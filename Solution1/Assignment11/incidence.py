#!/usr/bin/python3
import argparse
import re

#Vorkommen und Wahrscheinlichkeiten von Nukleotiden & kurzen Teilsequenzen
parser = argparse.ArgumentParser()
#nargs='+' für mehrere Eingaben nach
parser.add_argument('--sequence', nargs='+', required=True)
parser.add_argument('--genome', required=True)
args = parser.parse_args()
genome_file = args.genome
sequences = args.sequence


genome_sequence = []
#Ganze Sequenz als 1 zusammenfügen, Anfangszeile egal
with open(genome_file) as f:
    for line in f:
        if line.startswith('>'):
            continue
        genome_sequence.append(line.strip())

genome_sequence = "".join(genome_sequence)

results = []

for sequence in sequences:
    length = len(sequence)
    count = 0
    #länge der sequenz immer eins weiter nch rechst dass auch alle überlappenden gefunden werden
    for i in range(len(genome_sequence) - length + 1):
        #aktuelles fenster
        seq = genome_sequence[i:i+length]

        if seq == sequence:
            count += 1

    results.append(f"{sequence}: {count}")

print("\n".join(results))









