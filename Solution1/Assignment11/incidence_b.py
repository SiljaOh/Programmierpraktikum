#!/usr/bin/python3
import argparse
import sys

sequences = ["CTAG", "CG", "AACCCTGTC", "ATG"]
genome_path = ["C:/Users/gross/Documents/uni/ProPra/Bakterien/Mycobacterium_tuberculosis_CDC1551_genomic.fna", "C:/Users/gross/Documents/uni/ProPra/hs_chr08.fna"]

for path in genome_path:
    genome_parts = []
    current_seq = []

    with open(path, "r") as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                # Alte Sequenz von davor falls vorhaneden speichern
                if current_seq:
                    genome_parts.append("".join(current_seq))
                    current_seq = []
            # kein header = verlängern
            else:
                current_seq.append(line)
        # letzte line
        if current_seq:
            genome_parts.append("".join(current_seq))

    full_genome_sequence = "".join(genome_parts)
    gen_length = len(full_genome_sequence)

    # tatsächliche relative Häufigkeit
    rh = {
        'A': full_genome_sequence.count('A') / gen_length,
        'C': full_genome_sequence.count('C') / gen_length,
        'G': full_genome_sequence.count('G') / gen_length,
        'T': full_genome_sequence.count('T') / gen_length
    }

    # Erste Zeile der Ausgabe
    print(f"A: {rh['A'] * 100:.2f}% C: {rh['C'] * 100:.2f}% G: {rh['G'] * 100:.2f}% T: {rh['T'] * 100:.2f}%")

    for sequence in sequences:
        length = len(sequence)
        exp = gen_length * (0.25 ** length)
        # geg hier p = 1
        p_rh = 1.0
        # Anzahl erwartete Indizien
        for base in sequence:
            # get(base, 0) für alle Zeichen, die nicht A,C,G,T sind
            p_rh *= rh.get(base, 0)
        exp_p_rh = gen_length * p_rh

        # Tatsächliche Inzidenzen
        inc_count = 0
        for part in genome_parts:
            for i in range(len(part) - length + 1):
                # aktuelles fenster
                seq = part[i:i + length]
                if seq == sequence:
                    inc_count += 1

        # Fold-Change
        fold_change = inc_count / exp_p_rh

        print(f"{sequence}:\t{exp:.0f}\t{exp_p_rh:.0f}\t{inc_count}\t{fold_change:.1f}\t")




