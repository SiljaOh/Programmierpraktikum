#!/usr/bin/python3

import sys
import os

if len(sys.argv) < 2:
    print("Usage: python db_to_fasta.py <pfad_zur_datei.db>")
    sys.exit(1)

input_path = sys.argv[1]
#output path gleiche rordner  mit .fasta
output_path = os.path.splitext(input_path)[0] + ".fasta"

with open(input_path, 'r') as f_in, open(output_path, 'w') as f_out:
    for line in f_in:
        if line.startswith(">"):
            f_out.write(line)
        elif line.startswith("AS "):
            # Schneidet "AS " ab
            f_out.write(line[3:])

