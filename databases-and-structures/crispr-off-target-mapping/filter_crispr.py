#!/usr/bin/python3
import argparse
import re

#bowtie nutzen um mögliche off targets von CRISPR
#SAM = sequence alignment/map
parser = argparse.ArgumentParser()
parser.add_argument('--sam', type=argparse.FileType('r'), required=True)
parser.add_argument('--no-off-targets', type=argparse.FileType('w'), required=True)
parser.add_argument('--with-mismatch', type=argparse.FileType('w'), required=True)
args = parser.parse_args()

sam_file = args.sam

for line in sam_file:
    #head egal
    if line.startswith('@'):
        continue

    column = line.split('\t')
    qname = column[0]
    seq = column[9]
    #Zusatz spalten nach spalte 11
    tags = column[11:]
    flag = int(column[1])

    if flag == 4:
        #FASTA format
        args.no_off_targets.write(f">{qname}\n{seq}\n")
        continue

    missmatches = 0
    in_suffix = 0

    for t in tags:
        if t.startswith("MD:Z:"):
            #value weil in 3 spalten geteilt : separated
            val = t.split(":")[2]
            # Zahlen = Matches, Buchstaben = Mismatches
            position = 0
            buff = ""

            for v in val:
                if v.isdigit():
                    buff += v
                else:
                    #wenn buchstabe
                    if buff != "":
                        position += int(buff)
                        buff = ""

                    missmatches += 1

                    if position == 21 or position == 22:
                        in_suffix += 1

                    position += 1

            if missmatches <= 3 and in_suffix >= 1:
                args.with_mismatch.write(f">{qname}\n{seq}\n")
            break



