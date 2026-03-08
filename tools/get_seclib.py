#!/usr/bin/python3
import argparse
import os
import urllib.request

parser = argparse.ArgumentParser()
parser.add_argument('--DSSP', type=argparse.FileType('r', encoding='utf-8'), nargs='+')
parser.add_argument('--PDB', type=argparse.FileType('r', encoding='utf-8'), nargs='+')
parser.add_argument('--id', nargs='+')
# WICHTIG: 'a' Modus für das Sammeln aller CASP-Ordner
parser.add_argument('--output', default='casp_all.seclib', type=argparse.FileType('a', encoding='utf-8'))
args = parser.parse_args()

if not args.DSSP and not args.PDB and not args.id:
    parser.error("Mindestens --DSSP, --PDB oder --id angeben!")

seen_entries = set()

def write_entry(pdb_id, chain_id, full_aa, full_ss):
    # Eindeutiger Key (ID + Kette)
    key = (pdb_id.upper(), chain_id.strip())
    if key in seen_entries:
        return
    seen_entries.add(key)

    # GOR braucht Kontext: Wir filtern alles unter 30 Aminosäuren
    if len(full_aa) >= 30:
        args.output.write(f"> {pdb_id.upper()}{chain_id.strip()}00\n")
        args.output.write(f"AS {full_aa}\n")
        args.output.write(f"SS {full_ss}\n\n")

# --- DSSP VERARBEITUNG ---
if args.DSSP:
    for dssp in args.DSSP:
        # ID aus Dateiname (z.B. 4P7C aus 4P7C.cif.dssp)
        pdb_id = os.path.basename(dssp.name).split('.')[0].upper()[:4]
        chains = {}
        start = False

        for line in dssp:
            if '#  RESIDUE AA' in line:
                start = True
                continue
            if start and len(line) > 16:
                aa = line[13]
                if aa in ['!', '*']: continue
                chain_id = line[11].strip() or 'A'
                ss_raw = line[16]

                # 3-State Mapping
                if ss_raw in ['H', 'G', 'I']: ss = 'H'
                elif ss_raw in ['E', 'B']: ss = 'E'
                else: ss = 'C'

                if chain_id not in chains:
                    chains[chain_id] = {'aa': [], 'ss': []}
                chains[chain_id]['aa'].append(aa.upper())
                chains[chain_id]['ss'].append(ss)

        for chain_id, seqs in chains.items():
            write_entry(pdb_id, chain_id, ''.join(seqs['aa']), ''.join(seqs['ss']))
        print(f"Done DSSP: {pdb_id}")

# --- PDB VERARBEITUNG ---
aa_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLU': 'E',
    'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'MSE': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T',
    'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

if args.PDB:
    for pdb in args.PDB:
        pdb_id = os.path.basename(pdb.name).split('.')[0].upper()[:4]
        chains = {}
        helix_or_sheet = []

        for line in pdb:
            if line.startswith("ENDMDL"): break # Nur erstes Modell

            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                chain_id = line[21].strip() or 'A'
                res_id = line[22:26].strip()
                if chain_id not in chains:
                    chains[chain_id] = {'aa': [], 'res_ids': []}
                if res_name in aa_map:
                    if not chains[chain_id]['res_ids'] or chains[chain_id]['res_ids'][-1] != res_id:
                        chains[chain_id]['aa'].append(aa_map[res_name])
                        chains[chain_id]['res_ids'].append(res_id)

            if line.startswith("HELIX"):
                helix_or_sheet.append(('H', line[19], line[21:25].strip(), line[33:37].strip()))
            elif line.startswith("SHEET"):
                helix_or_sheet.append(('E', line[21], line[22:26].strip(), line[33:37].strip()))

        for chain_id, seqs in chains.items():
            ss_seq = ['C'] * len(seqs['aa'])
            for typ, chain, start, end in helix_or_sheet:
                if chain != chain_id: continue
                for i, res_id in enumerate(seqs['res_ids']):
                    try:
                        r_val = int(''.join(filter(str.isdigit, res_id)))
                        s_val = int(''.join(filter(str.isdigit, start)))
                        e_val = int(''.join(filter(str.isdigit, end)))
                        if s_val <= r_val <= e_val:
                            ss_seq[i] = typ
                    except ValueError: pass

            write_entry(pdb_id, chain_id, ''.join(seqs['aa']), ''.join(ss_seq))
        print(f"Done PDB: {pdb_id}")