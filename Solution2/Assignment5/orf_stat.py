#!/usr/bin/env python3
import argparse, sys

def read_fasta(path):
    h, s = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(s).upper()
                h, s = line[1:].strip(), []
            else:
                s.append(line)
        if h is not None:
            yield h, "".join(s).upper()

def aa_length_from_orf_dna(dna):
    if dna is None:
        return -10**9
    n = len(dna)
    if n < 6 or n % 3 != 0:
        return -10**9
    return (n // 3) - 1

def make_bins(lower, upper, bin_size):
    bins = []
    x = lower
    while x <= upper:
        bins.append((x, min(upper, x + bin_size - 1)))
        x += bin_size
    return bins

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--histogram")
    p.add_argument("--lower", type=int, required=True)
    p.add_argument("--upper", type=int, required=True)
    p.add_argument("--bin_size", type=int, required=True)
    args = p.parse_args()

    lower, upper, bs = args.lower, args.upper, args.bin_size
    bins = make_bins(lower, upper, bs)
    counts = [0] * len(bins)

    total = 0
    for _, seq in read_fasta(args.fasta):
        aa = aa_length_from_orf_dna(seq)
        if lower <= aa <= upper:
            total += 1
            for i, (a, b) in enumerate(bins):
                if a <= aa <= b:
                    counts[i] += 1
                    break

    sys.stdout.write(str(total) + "\n")

    if args.histogram:
        with open(args.histogram, "w") as out:
            for (a, b), c in zip(bins, counts):
                out.write("%d-%d\t%d\n" % (a, b, c))

if __name__ == "__main__":
    main()