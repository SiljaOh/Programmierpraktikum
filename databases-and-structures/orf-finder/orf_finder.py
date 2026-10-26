#!/usr/bin/env python3
import argparse
import sys

STOP = {"TAA", "TAG", "TGA"}

def read_fasta(path):
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq).upper()
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq).upper()

def revcomp(s):
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(comp.get(b,"N") for b in reversed(s))

def find_orfs(seq):
    orfs = []
    n = len(seq)
    for frame in (0,1,2):
        start = None
        for i in range(frame, n-2, 3):
            codon = seq[i:i+3]
            if start is None:
                if codon == "ATG":
                    start = i
            else:
                if codon in STOP:
                    orfs.append(seq[start:i+3])
                    start = None
    return orfs

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--output")
    args = p.parse_args()

    out = open(args.output, "w") if args.output else sys.stdout

    for header, seq in read_fasta(args.fasta):
        idx = 0
        for orf in find_orfs(seq) + find_orfs(revcomp(seq)):
            out.write(">%s_%d\n%s\n" % (header, idx, orf))
            idx += 1

    if args.output:
        out.close()

if __name__ == "__main__":
    main()