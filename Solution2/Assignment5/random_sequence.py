#!/usr/bin/env python3
import argparse, sys, random

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--length", type=int, required=True)
    p.add_argument("--gc", type=float, required=True)
    p.add_argument("--name", default="random")
    args = p.parse_args()

    L = args.length
    p_gc = args.gc
    if L < 0:
        raise SystemExit("length must be >= 0")
    if p_gc < 0.0 or p_gc > 1.0:
        raise SystemExit("gc must be between 0 and 1")

    seq = []
    for _ in range(L):
        if random.random() < p_gc:
            seq.append("G" if random.random() < 0.5 else "C")
        else:
            seq.append("A" if random.random() < 0.5 else "T")

    sys.stdout.write(">%s\n%s\n" % (args.name, "".join(seq)))

if __name__ == "__main__":
    main()