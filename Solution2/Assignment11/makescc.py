#!/usr/bin/env python3

import argparse
import urllib.request
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import math

# makescc.py (Aufgabe 11)


AA3_TO_AA1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


@dataclass
class ResidueAtom:
    chain: str
    resseq: int
    serial: int
    aa1: str
    x: float
    y: float
    z: float
    ss: str = "C"   # H (Helix), E (Extended for beta-sheet), C (Coil)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate sscc table for a PDB id or local PDB file")
    p.add_argument("--id", required=True, help="PDB id (e.g. 1CRN) or local path to .pdb file")
    p.add_argument("--distance", required=True, type=float, help="contact cutoff distance in Angstrom")
    p.add_argument("--type", required=True, dest="atom_type", help="atom name, e.g. CA or CB")
    p.add_argument("--length", required=True, type=int, help="local contact max sequence distance")
    p.add_argument("--contactmatrix", default=None, help="optional: write contact matrix to this file")
    return p.parse_args()


# Funktion 1) download_pdb
def download_pdb(pdb_id_or_path: str) -> List[str]:
    import os
    if os.path.isfile(pdb_id_or_path):
        with open(pdb_id_or_path, "r", encoding="utf-8", errors="replace") as f:
            return f.readlines()

    pdb_id = pdb_id_or_path.strip().upper()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

    with urllib.request.urlopen(url) as response:
        text = response.read().decode("utf-8", errors="replace")

    return text.splitlines(keepends=True)


# Funktion 2) extract_ss_ranges
def extract_ss_ranges(
    pdb_lines: List[str]
) -> Tuple[List[Tuple[str, int, int]], List[Tuple[str, int, int]]]:
    # HELIX und SHEET Bereiche aus Header lesen
    helix: List[Tuple[str, int, int]] = []
    sheet: List[Tuple[str, int, int]] = []

    for line in pdb_lines:
        rec = line[0:6].strip()

        if rec == "HELIX":
            chain = line[19:20].strip()
            start_s = line[21:25].strip()
            end_s = line[33:37].strip()
            if chain and start_s and end_s:
                helix.append((chain, int(start_s), int(end_s)))

        elif rec == "SHEET":
            chain = line[21:22].strip() or "?"
            start_s = line[22:26].strip()
            end_s = line[33:37].strip()
            if chain and start_s and end_s:
                sheet.append((chain, int(start_s), int(end_s)))

    return helix, sheet


# Funktion 3) keep_only_model_one
def keep_only_model_one(pdb_lines: List[str]) -> List[str]:
    # wenn MODEL existiert, nur Modell 1 behalten, sonst alles
    if not any(line.startswith("MODEL") for line in pdb_lines):
        return pdb_lines

    out: List[str] = []
    in_first = False

    for line in pdb_lines:
        if line.startswith("MODEL"):
            in_first = (line[10:14].strip() == "1")
            continue

        if line.startswith("ENDMDL"):
            if in_first:
                break
            continue

        if in_first:
            out.append(line)

    return out


# Funktion 4) parse_residue_atoms
def parse_residue_atoms(
    pdb_lines: List[str],
    atom_type: str,
    aa3_to_aa1: dict
) -> List[ResidueAtom]:
    # pro Residue genau ein Atom sammeln, zB CA
    want = atom_type.strip()
    residues: List[ResidueAtom] = []
    seen = set()

    for line in pdb_lines:
        if not line.startswith("ATOM"):
            continue

        if line[12:16].strip() != want:
            continue

        alt_loc = line[16:17]
        if alt_loc not in (" ", "A"):
            continue

        resname3 = line[17:20].strip()
        chain = line[21:22].strip()
        resseq_s = line[22:26].strip()
        serial_s = line[6:11].strip()

        if not resseq_s or not serial_s:
            continue

        key = (chain, int(resseq_s))
        if key in seen:
            continue
        seen.add(key)

        aa1 = aa3_to_aa1.get(resname3, "X")

        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())

        residues.append(
            ResidueAtom(chain, int(resseq_s), int(serial_s), aa1, x, y, z)
        )

    return residues


# Funktion 5) ss_for_residue
def ss_for_residue(
    chain: str,
    resseq: int,
    helix_ranges,
    sheet_ranges
) -> str:
    # H falls in HELIX Bereich, sonst E falls in SHEET Bereich, sonst C
    for ch, start, end in helix_ranges:
        if ch == chain and start <= resseq <= end:
            return "H"

    for ch, start, end in sheet_ranges:
        if ch == chain and start <= resseq <= end:
            return "E"

    return "C"


# Funktion 6) compute_contacts
def compute_contacts(residues, distance_cutoff, local_length):
    # Kontakt wenn Distanz kleiner cutoff
    n = len(residues)
    global_counts = [0] * n
    local_counts = [0] * n

    contact_matrix = [[0] * n for _ in range(n)]
    distance_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        ai = residues[i]
        for j in range(i + 1, n):
            aj = residues[j]

            dx = ai.x - aj.x
            dy = ai.y - aj.y
            dz = ai.z - aj.z
            dist = math.sqrt(dx * dx + dy * dy + dz * dz)

            distance_matrix[i][j] = dist
            distance_matrix[j][i] = dist

            if dist < distance_cutoff:
                contact_matrix[i][j] = 1
                contact_matrix[j][i] = 1

                same_chain = ai.chain == aj.chain
                seq_gap = abs(ai.resseq - aj.resseq)

                if same_chain and 0 < seq_gap < local_length:
                    local_counts[i] += 1
                    local_counts[j] += 1
                else:
                    global_counts[i] += 1
                    global_counts[j] += 1

    return global_counts, local_counts, contact_matrix, distance_matrix


# Funktion 7) write_contact_matrix (nur wenn --contactmatrix benutzt wird)
def write_contact_matrix(path, residues, contact_matrix):
    labels = [f"{r.chain}{r.resseq}" for r in residues]

    with open(path, "w", encoding="utf-8") as f:
        f.write("\t" + "\t".join(labels) + "\n")
        for i, row in enumerate(contact_matrix):
            f.write(labels[i] + "\t" + "\t".join(str(v) for v in row) + "\n")


def main() -> int:

    # 1) Argumente/Variable lesen
    args = parse_args()

    # 2) Werte in Variable schreiben
    pdb_id = args.id
    distance_cutoff = args.distance
    atom_type = args.atom_type
    local_length = args.length
    contactmatrix_path = args.contactmatrix

    # 3) PDB laden
    pdb_lines_all = download_pdb(pdb_id)

    # 4) Sekundaerstruktur Bereiche auslesen
    helix_ranges, sheet_ranges = extract_ss_ranges(pdb_lines_all)

    # 5) Modell 1 ablesen
    pdb_lines_model_one = keep_only_model_one(pdb_lines_all)

    # 6) Residues sammeln
    residues = parse_residue_atoms(
        pdb_lines_model_one,
        atom_type,
        AA3_TO_AA1
    )

    # 7) ss pro Residue setzen
    for r in residues:
        r.ss = ss_for_residue(r.chain, r.resseq, helix_ranges, sheet_ranges)

    # 8) Kontakte berechnen
    global_counts, local_counts, contact_matrix, _ = compute_contacts(
        residues,
        distance_cutoff,
        local_length
    )

    # 9) sscc ausgeben
    print("chain\tpos\tserial\taa\tss\tglobal\tlocal")
    for i, r in enumerate(residues):
        print(
            f"{r.chain}\t{r.resseq}\t{r.serial}\t"
            f"{r.aa1}\t{r.ss}\t"
            f"{global_counts[i]}\t{local_counts[i]}"
        )

    # 10) Dateien schreiben
    if contactmatrix_path is not None:
        write_contact_matrix(contactmatrix_path, residues, contact_matrix)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())