#!/usr/bin/python3
import argparse
import os
import urllib.request

parser = argparse.ArgumentParser()
#entweder oder
parser.add_argument('--DSSP', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--PDB', type=argparse.FileType('r'), nargs='+')
parser.add_argument('--id', nargs='+')
parser.add_argument('--output', default='output.seclib', type=argparse.FileType('w'))
args = parser.parse_args()

#code duplis vermeidne

if not args.DSSP and not args.PDB and not args.id:
    parser.error("Mindestens entweder --DSSP oder --PDB angeben!")

#extensions= dateiendungen
def find_local_file(id, extensions):
    #für gross un dklein testen
    search_file = [f"{id.lower()}.{ext}" for ext in extensions] + [f"{id.upper()}.{ext}" for ext in extensions]
    #aktueller ordnerpfad, unterordner, dateien,
    # os.walk ale ordner ab aktuelles verzeichnis, rekursiv
    for root, dirs, files in os.walk('.'):
        for fname in files:
            for name in search_file:
                if fname.lower() == name.lower():
                    return os.path.join(root, fname)
    return None

result = ""
seen_entries = set()
downloaded_files = []

if args.id:
    for id in args.id:
        #upper falls klein geschrieben
        local_dssp = find_local_file(id, ['dssp'])
        if local_dssp:
            if args.DSSP is None:
                #kien eingelesen aber lokal vorhanden -> jetzt hinzufügen
                #erst leer erstellen
                args.DSSP = []
            args.DSSP.append(open(local_dssp, 'r'))
            continue  # kein PDB nötig, wenn DSSP vorhanden

        local_pdb = find_local_file(id, ['pdb'])
        if local_pdb:
            if args.PDB is None:
                args.PDB = []
            args.PDB.append(open(local_pdb, 'r'))
            continue

        #sonst runterlsden
        url = f"https://files.rcsb.org/download/{id.upper()}.pdb"
        filename = f"{id.lower()}.pdb"
        urllib.request.urlretrieve(url, filename)
        downloaded_files.append(open(filename, 'r'))

    if args.PDB:
        args.PDB = args.PDB + downloaded_files
    else:
        if downloaded_files:
            args.PDB = downloaded_files
        else: args.PDB = None



if args.DSSP:
    dssp_files = args.DSSP
    for dssp in dssp_files:

        # chain_id -> {'aa': [], 'ss': []}
        chains = {}
        start = False

        for line in dssp:
            line = line.strip()

            if line.startswith("HEADER"):
                parts = line.split()
                #lteztes Element ist ID
                pdb_id = parts[-1]

            if line.startswith('  #  RESIDUE AA'):
                start = True
                #passt, danach starten
                continue

            if start:
                #aa sequenz nacheinander in zeilen, ss auch
                aa = line[13]
                ss_before = line[16]
                chain_id = line[11]
                #Helix: H G I
                #Sheet: E B
                #C = Rest
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

            entry_key = (pdb_id.upper(), chain_id)
            #duplikate vermeiden
            if entry_key in seen_entries:
                continue
            seen_entries.add(entry_key)

            full_aa = ''.join(seqs['aa'])
            full_ss = ''.join(seqs['ss'])
            result += f"> {pdb_id}{chain_id}00\nAS {full_aa}\nSS {full_ss}\n\n"

aa_map = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
          'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
          'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
          'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'}

if args.PDB:
    pdb_files = args.PDB
    for pdb in pdb_files:

        # chain_id -> {'aa': [], 'ss': [], 'residue_ids': []
        chains = {}
        helix_or_sheet = []

        for line in pdb:
            line = line.strip()
            if line.startswith("HEADER"):
                parts = line.split()
                # lteztes Element ist ID
                pdb_id = parts[-1]

            if line.startswith("ATOM"):
                res_name = line[17:20].strip()  # zb GLU
                chain_id = line[21].strip() #A, B, C, oder so
                res_id = line[22:26].strip()  # danach zb 2

                #buchstabe plus nummer
                current_id = (chain_id, res_id)

                if chain_id not in chains:
                    chains[chain_id] = {'aa': [], 'residue_ids': [], 'last_id': None}

                if current_id != chains[chain_id]['last_id']:
                    #es kann sowas wie MSE (modifiziert MET) deswegen absicherung
                    if res_name in aa_map:
                        # wenn neue ID -> bei nächster AS, AS hinzufügen
                        chains[chain_id]['aa'].append(aa_map[res_name])
                        chains[chain_id]['residue_ids'].append(current_id)
                        chains[chain_id]['last_id'] = current_id

            if line.startswith("HELIX"):
                #zb von A (chain)  44 bis 55
                helix_or_sheet.append(('H', line[19], line[21:25].strip(), line[33:37].strip()))
            elif line.startswith("SHEET"):
                helix_or_sheet.append(('E', line[21], line[22:26].strip(), line[33:37].strip()))

        #weil rest ist C
        for chain_id, seqs in chains.items():
            #dupllikate vermieden
            entry_key = (pdb_id.upper(), chain_id)
            if entry_key in seen_entries:
                continue
            seen_entries.add(entry_key)

            ss_seq = ['C'] * len(seqs['aa'])
            for typ, chain, start, end in helix_or_sheet:
                if chain != chain_id:
                    continue
                #sucht in der Liste residue_ids nach Tupel (chain, start) und gibt Position zurück
                start_idx = seqs['residue_ids'].index((chain, start))
                end_idx = seqs['residue_ids'].index((chain, end))
                for i in range(start_idx, end_idx + 1):
                    ss_seq[i] = typ

            full_aa = ''.join(seqs['aa'])
            full_ss = ''.join(ss_seq)
            result += f"> {pdb_id}{chain_id}00\nAS {full_aa}\nSS {full_ss}\n\n"


args.output.write(result)
#runtergeladene pdb löschen
for f in downloaded_files:
    f.close()
    os.remove(f.name)

