#!/usr/bin/python3
import argparse
import os
import urllib.request

parser = argparse.ArgumentParser()
parser.add_argument('--DSSP', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--PDB', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--id', nargs='+')
parser.add_argument('--output', default='output.seclib', type=argparse.FileType('w'))
args = parser.parse_args()

if not args.DSSP and not args.PDB and not args.id:
    parser.error("Mindestens entweder --DSSP oder --PDB angeben!")

def find_local_file(id, extensions):
    search_file = [f"{id.lower()}.{ext}" for ext in extensions] + [f"{id.upper()}.{ext}" for ext in extensions]
    for root, dirs, files in os.walk('.'):
        for fname in files:
            for name in search_file:
                if fname.lower() == name.lower():
                    return os.path.join(root, fname)
    return None

def write_entry(pdb_id, chain_id, full_aa, full_ss):
    key = (pdb_id.upper(), chain_id)
    if key in seen_entries:
        return
    seen_entries.add(key)
    args.output.write(f"> {pdb_id}{chain_id}00\nAS {full_aa}\nSS {full_ss}\n\n")

seen_entries = set()
downloaded_files = []

if args.id:
    for id in args.id:
        local_dssp = find_local_file(id, ['dssp'])
        if local_dssp:
            if args.DSSP is None:
                args.DSSP = []
            args.DSSP.append(open(local_dssp, 'r'))
            continue

        local_pdb = find_local_file(id, ['pdb'])
        if local_pdb:
            if args.PDB is None:
                args.PDB = []
            args.PDB.append(open(local_pdb, 'r'))
            continue

        url = f"https://files.rcsb.org/download/{id.upper()}.pdb"
        filename = f"{id.lower()}.pdb"
        urllib.request.urlretrieve(url, filename)
        downloaded_files.append(open(filename, 'r'))

    if downloaded_files:
        args.PDB = (args.PDB or []) + downloaded_files

if args.DSSP:
    for dssp in args.DSSP:
        chains = {}
        start = False

        for line in dssp:
            line = line.strip()
            if line.startswith("HEADER"):
                pdb_id = line.split()[-1]
            if line.startswith('  #  RESIDUE AA'):
                start = True
                continue
            if start:
                aa = line[13]
                ss_before = line[16]
                chain_id = line[11]
                if ss_before in ['H', 'G', 'I']:
                    ss = 'H'
                elif ss_before in ['E', 'B']:
                    ss = 'E'
                else:
                    ss = 'C'

                if chain_id not in chains:
                    chains[chain_id] = {'aa': [], 'ss': []}
                chains[chain_id]['aa'].append(aa)
                chains[chain_id]['ss'].append(ss)

        for chain_id, seqs in chains.items():
            write_entry(pdb_id, chain_id, ''.join(seqs['aa']), ''.join(seqs['ss']))

aa_map = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M',  # MSE -> MET minimal Änderung
    'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

if args.PDB:
    for pdb in args.PDB:

        if 'pdb_id' not in locals():
            pdb_id = os.path.basename(pdb.name).split('.')[0].upper()

        chains = {}
        helix_or_sheet = []

        for line in pdb:
            line = line.strip()
            if line.startswith("ENDMDL"):
                break
            if line.startswith("MODEL") and line != "MODEL 1":
                break

            if line.startswith("HEADER"):
                pdb_id = line.split()[-1]

            if line.startswith("ATOM"):
                res_name = line[17:20].strip()
                chain_id = line[21].strip()
                res_id = line[22:26].strip()
                current_id = (chain_id, res_id)

                if chain_id not in chains:
                    chains[chain_id] = {'aa': [], 'residue_ids': [], 'last_id': None}

                if current_id != chains[chain_id]['last_id']:
                    if res_name in aa_map:
                        chains[chain_id]['aa'].append(aa_map[res_name])
                        chains[chain_id]['residue_ids'].append(current_id)
                        chains[chain_id]['last_id'] = current_id

            if line.startswith("HELIX"):
                helix_or_sheet.append(('H', line[19], line[21:25].strip(), line[33:37].strip()))
            elif line.startswith("SHEET"):
                helix_or_sheet.append(('E', line[21], line[22:26].strip(), line[33:37].strip()))

        for chain_id, seqs in chains.items():
            ss_seq = ['C'] * len(seqs['aa'])
            for typ, chain, start, end in helix_or_sheet:
                if chain != chain_id:
                    continue
                for i, (c, res_id) in enumerate(seqs['residue_ids']):
                    if c == chain:
                        try:
                            # wir behandeln Insertionscodes minimal korrekt
                            r_id_int = int(''.join(filter(str.isdigit, res_id)))
                            s_id_int = int(''.join(filter(str.isdigit, start)))
                            e_id_int = int(''.join(filter(str.isdigit, end)))
                            if s_id_int <= r_id_int <= e_id_int:
                                ss_seq[i] = typ
                        except ValueError:
                            pass

            write_entry(pdb_id, chain_id, ''.join(seqs['aa']), ''.join(ss_seq))

for f in downloaded_files:
    f.close()
    os.remove(f.name)