#!/usr/bin/python3

import re
import sys

def parse_summary(filepath):
    metrics = {}
    with open(filepath, 'r') as f:
        for line in f:
            # Zeilen wie: "q3      : Mean:   60.1   Dev: ..."
            match = re.match(r'\s*(\S+)\s*:\s*Mean:\s*([\d.]+)', line)
            if match:
                name = match.group(1)
                mean = float(match.group(2))
                metrics[name] = mean
    return metrics

files = {'GOR1': 'summary_casp_gor1.txt', 'GOR3': 'summary_casp_gor3.txt', 'GOR4': 'summary_casp_gor4.txt'}

prefix = sys.argv[1] if len(sys.argv) > 1 else ''

data = {}
for name, fname in files.items():
    path = prefix + fname
    try:
        data[name] = parse_summary(path)
    except FileNotFoundError:
        print(f"Datei nicht gefunden: {path}")
        sys.exit(1)

metrics_order = ['q3', 'qObs_H', 'qObs_E', 'qObs_C', 'SOV', 'SOV_H', 'SOV_E', 'SOV_C']

# Ausgabe
print(f"{'Metric':<10} {'GOR1':>8} {'GOR3':>8} {'GOR4':>8}")

for m in metrics_order:
    row = f"{m:<10}"
    for gor in ['GOR1', 'GOR3', 'GOR4']:
        val = data[gor].get(m, float('nan'))
        row += f" {val:>8.1f}"
    print(row)
print("=" * 42)