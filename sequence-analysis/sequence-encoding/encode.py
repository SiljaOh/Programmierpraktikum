#!/usr/bin/env python3
import os
import re
import argparse


#read csv, convert to dict [{collumn1: value1, collumn2: value2, ...} , {...} , ...]
def parse_csv(filepath):
    data_list = []

    with open(filepath, 'r') as f:
        header = f.readline().strip().split(',')

        for line in f:
            line = line.strip()
            if not line:
                continue
            values = line.split(',')
            row_dict = dict(zip(header, values))
            data_list.append(row_dict)

    return data_list


#produces tsv: Data_Type \t Count of Data type
#No header
#example: Integ Cluster   1
def write_exptypes(data,output_file):
    counts = {} #for all Data types

    for row in data:
        dt = row['Data_Type']
        counts[dt] = counts.get(dt, 0) + 1

    with open(output_file, 'w') as f:
        for dt in sorted(counts.keys()):
            f.write(f"{dt}\t{counts[dt]}\n")


#Cell_Type = Zelllinie
#produces tsv: Cell_Type \t Count of Antibodies for ChIP-Seq Experiments
#No header
#inlcude Cell_Type s without Chip Seq Experiment
#example: 8988T   0 (--> No experiment with ChIP-seq)
#interpretation of generated file: a celltype has ChIp-seq experiments with n different Antibodies
def write_antibodies(data,output_file):

    all_cells = sorted(list(set(row['Cell_Type'] for row in data))) #set eliminates duplicates
    # { cell_type : {antibody1, antibody2, ...} } dictionary with set als values
    cell_antibodies = {cell: set() for cell in all_cells} #set!!!

    for row in data:
        if row['Data_Type'] == 'ChIP-seq':
            cell = row['Cell_Type']
            factors = row['Experimental_Factors']

            # Extract value after "Antibody="
            match = re.search(r'Antibody=(\S+)', factors)
            if match:
                cell_antibodies[cell].add(match.group(1))

    with open(output_file, 'w') as f:
        for cell in all_cells:
            count = len(cell_antibodies[cell])
            f.write(f"{cell}\t{count}\n")



#produces tsv: cell line \t RNAseq Accession \t ChIPseq Accession
#DCC number = DCC_Accession = eindeutige ID für jedes Expriment
#with header
#sammeln von cell lines und der DCC_Accession von RNAseq experiment und von ChIP-seq experiment für H3K27me3
#überpürufen ob jede cell line mind eine RNA Accession hat und mind eine ChIP Accession für H3K27me3
#sortieren, wenn es mehrere IDs gibt für mehrere Experimente

def write_chip_rna_seq(data, output_file):
    # Dictionary with list: { "cell type": {"rna": [], "chip": []} } #use list instead of set to avoid losing duplicates
    res = {}

    for row in data:
        cell = row['Cell_Type']
        dcc = row['DCC_Accession'].strip()
        dtype = row['Data_Type'].strip()
        factors = row['Experimental_Factors']

        if not dcc:  # skip empty accessions!!!!!
            continue

        #DCC-Accession Duplicates are allowed
        #Zelllinien cant be duplicates:
        if cell not in res:
            res[cell] = {'rna': [], 'chip': []}

        if dtype == 'RNA-seq':
            res[cell]['rna'].append(dcc)

        if dtype == 'ChIP-seq':
            match = re.search(r'Antibody=(\S+)', factors)
            if match and match.group(1) == 'H3K27me3':
                res[cell]['chip'].append(dcc)

    with open(output_file, 'w') as f:
        #header
        f.write("cell line\tRNAseq Accession\tChIPseq Accession\n")

        for cell in sorted(res.keys()):
            if res[cell]['rna'] and res[cell]['chip']: #Sicherstellung dass IDS für RNA seq und ChIP seq vorhanden sind
                rna_ids = ",".join(sorted(res[cell]['rna']))
                chip_ids = ",".join(sorted(res[cell]['chip']))
                f.write(f"{cell}\t{rna_ids}\t{chip_ids}\n")


#Header encode.csv:
#Data_Type,Cell_Type,Experimental_Factors,Treatment,Lab,PI,Assembly,Status,GEO_Accession,DCC_Accession,Date_Unrestricted
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    data = parse_csv(args.input)
    write_exptypes(data, os.path.join(args.output, "exptypes.tsv"))
    write_antibodies(data, os.path.join(args.output, "antibodies.tsv"))
    write_chip_rna_seq(data, os.path.join(args.output, "chip_rna_seq.tsv"))

if __name__ == "__main__":
    main()

#Ausgabe: encode.py --input <path_to_file> --output <path_to_dir>
#Creates files:exptypes.tsv  antibodies.tsv  chip_rna_seq.tsv in output dir