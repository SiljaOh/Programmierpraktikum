#!/usr/bin/python3
import argparse
import re
import sys

#Vorkommen und Wahrscheinlichkeiten von Nukleotiden & kurzen Teilsequenzen in Genom
parser = argparse.ArgumentParser()
#nargs='+' für mehrere Eingaben nach
parser.add_argument('--sequence', nargs='+', required=True)
parser.add_argument('--genome', type=argparse.FileType('r'), required=True)
args = parser.parse_args()
genome_file = args.genome
sequences = args.sequence
sequences = [s.upper() for s in sequences]

#Ungültige Sequenzen vermeiden
for s in sequences:
    if not re.match(r"^[ACGKTUN]+$", s.upper()):
        print("Ungültige Sequenz")
        sys.exit(1)

genome_parts = []
current_seq = []
#Ganze Sequenz als 1 zusammenfügen, Anfangszeile egal
for line in genome_file:
    line = line.strip()
    if line.startswith(">"):
        # Alte Sequenz von davor falls vorhaneden speichern
        if current_seq:
            genome_parts.append("".join(current_seq).upper())
            current_seq = []
        # kein header = verlängern
    else:
        current_seq.append(line)
    # letzte line
if current_seq:
    genome_parts.append("".join(current_seq).upper())


results = []

for sequence in sequences:
    length = len(sequence)
    count = 0
    for part in genome_parts:
        for i in range(len(part) - length + 1):
            # aktuelles fenster
            seq = part[i:i + length]
            if seq == sequence:
                count += 1

    results.append(f"{sequence}: {count}")

print("\n".join(results))

