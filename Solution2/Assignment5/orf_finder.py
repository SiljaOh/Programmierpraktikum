#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


@dataclass(frozen=True)
class FastaRecord:
    header: str
    seq: str


def read_fasta(path: str | Path) -> Iterator[FastaRecord]:

    header: Optional[str] = None
    chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                # flush previous record
                if header is not None:
                    seq = "".join(chunks).replace(" ", "").upper()
                    yield FastaRecord(header=header, seq=seq)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line)

    # flush last record
    if header is not None:
        seq = "".join(chunks).replace(" ", "").upper()
        yield FastaRecord(header=header, seq=seq)


def revcomp(seq: str) -> str:

    comp: Dict[str, str] = {
        "A": "T",
        "T": "A",
        "C": "G",
        "G": "C",
        "N": "N",
    }
    # map unknowns to N
    return "".join(comp.get(b, "N") for b in reversed(seq.upper()))


def find_orfs_in_one_strand(seq: str) -> List[str]:

    orfs: List[str] = []
    n = len(seq)

    for frame in (0, 1, 2):
        current_start: Optional[int] = None  # base index of the ATG within seq

        # step codon-wise in this frame
        for i in range(frame, n - 2, 3):
            codon = seq[i : i + 3]
            if current_start is None:
                if codon == START_CODON:
                    current_start = i
            else:
                if codon in STOP_CODONS:
                    # ORF includes the stop codon
                    orf_seq = seq[current_start : i + 3]
                    orfs.append(orf_seq)
                    current_start = None
                # else: keep scanning; ignore internal ATGs by design

    return orfs


def format_fasta(header: str, seq: str, width: int = 80) -> str:
    """
    FASTA formatting with optional line wrapping.
    """
    if width <= 0:
        return f">{header}\n{seq}\n"
    lines = [seq[i : i + width] for i in range(0, len(seq), width)]
    return f">{header}\n" + "\n".join(lines) + "\n"


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Find ORFs (ATG..Stop) in 6 frames and output as FASTA.")
    p.add_argument("--fasta", required=True, help="Input FASTA file (can contain multiple records).")
    p.add_argument("--output", default=None, help="Output FASTA file (default: stdout).")
    # --db is intentionally not implemented per user request
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    out_handle = None
    try:
        if args.output:
            out_handle = open(args.output, "w", encoding="utf-8")
            out = out_handle
        else:
            out = sys.stdout

        orf_counter_global = 0

        for rec in read_fasta(args.fasta):
            # Find ORFs on forward strand
            forward_orfs = find_orfs_in_one_strand(rec.seq)

            # Find ORFs on reverse-complement strand
            rc = revcomp(rec.seq)
            reverse_orfs = find_orfs_in_one_strand(rc)

            # Output all ORFs for this record, numbered sequentially per input record
            idx = 0
            for orf_seq in forward_orfs + reverse_orfs:
                header = f"{rec.header}_{idx}"
                out.write(format_fasta(header, orf_seq, width=80))
                idx += 1
                orf_counter_global += 1

        return 0

    except FileNotFoundError:
        print(f"ERROR: Input FASTA not found: {args.fasta}", file=sys.stderr)
        return 2
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    finally:
        if out_handle is not None:
            out_handle.close()


if __name__ == "__main__":
    raise SystemExit(main())