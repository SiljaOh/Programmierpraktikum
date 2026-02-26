#!/usr/bin/python3

def parse_fasta(file):
    sequences = {}
    current_header = None
    current_seq = []

    for line in file:
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            if current_header:
                sequences[current_header] = "".join(current_seq)
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)

    if current_header:
        sequences[current_header] = "".join(current_seq)

    return sequences


def parse_fasta_content(content):
    """Für Web: content aus upload"""
    sequences = {}
    current_header = None
    current_seq = []

    for line in content.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                sequences[current_header] = "".join(current_seq)
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)

    if current_header:
        sequences[current_header] = "".join(current_seq)
    return sequences