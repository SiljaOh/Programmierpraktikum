#!/usr/bin/env python3
from collections import defaultdict
import argparse, sys

#fasta file einlesen
def read_fasta(path):
    rec, h, buf = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    rec.append((h, "".join(buf)))
                h, buf = line[1:], []
            else:
                buf.append(line)
    if h is not None:
        rec.append((h, "".join(buf)))
    return rec

def unique_genes(records, k, start):
    if k <= 0 or start < -1:
        return 0
    m = defaultdict(set)

    for gid, (_, seq) in enumerate(records):
        s = seq.upper()
        if start == -1:
            if len(s) < k:
                continue
            for i in range(len(s) - k + 1):
                m[s[i:i+k]].add(gid)
        else:
            if len(s) >= start + k:
                m[s[start:start+k]].add(gid)

    hit = [0] * len(records)
    for genes in m.values():
        if len(genes) == 1:
            hit[next(iter(genes))] = 1
    return sum(hit)

def main():
    p = argparse.ArgumentParser(add_help=False)
    p.add_argument("--fasta", required=True)
    p.add_argument("--k", nargs="+", type=int, required=True)
    p.add_argument("--start", type=int, default=-1)  # -1 = Teil A
    p.add_argument("--plot", action="store_true")
    p.add_argument("--plot-out", default="unique_plot.png")
    a, _ = p.parse_known_args()

    rec = read_fasta(a.fasta)
    n=len(rec)
    ks=[]
    percents=[]
    for k in a.k:
        u = unique_genes(rec, k, a.start)
        ks.append(k)
        percents.append((u / n * 100.0) if n else 0.0)
        sys.stdout.write(f"{k}\t{u}\n")

    if a.plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(ks, percents, marker="o")
        plt.xlabel("k" if a.start == -1 else "Basen vom Genanfang (k)")
        plt.ylabel("Prozent eindeutig identifizierbarer Gene")
        title = "Teil a: alle k-mers" if a.start == -1 else f"Teil b: Startposition s={a.start}"
        plt.title(title)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(a.plot_out, dpi=200)

if __name__ == "__main__":
    main()