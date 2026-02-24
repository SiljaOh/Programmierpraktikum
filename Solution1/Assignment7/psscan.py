#!/usr/bin/env python3
import sys
import argparse
import re
import urllib.request #for --web

#read fasta, converting to Dic {Header : Seq, Header : Seq, ... }
def parse_fasta(file):
    sequences = {}
    current_header = None
    current_seq=[]

    for line in file:
        line = line.strip() #remove any invisible characters
        if not line: continue  #remove empty lines

        if line.startswith(">"):
            if current_header:
                sequences[current_header] = "".join(current_seq)

            #new header without >
            current_header = line[1:]
            current_seq=[]

        else:
            current_seq.append(line)

    #save last sequence after iteration
    if current_header:
        sequences[current_header] = "".join(current_seq)

    return sequences

#get arguments
#argument looks like this:
# cat FASTA.seq | psscan.py --fasta - [--pattern <pattern>|--web ←-,→ <ProSiteID>] [--extern]<ID1>\t<startPosition>\t<matchedSequence><ID2>\t<startPosition>\t<matchedSequence>
def main():

    parser = argparse.ArgumentParser(description="prosite pattern scanner")

    #--fasta: file or - Stdin, required
    parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                        help="file path or '-' for Pipe")

    #--pattern: String, optional
    parser.add_argument('--pattern', type=str, help="PROSITE Pattern")

    #--web: Erweiterung: instead of string pattern, use PROSITE ID, optional
    parser.add_argument('--web', type=str, help="PROSITE-ID")

    #BONUS TO DO
    # --extern: remote scan with prosite website, optional
    # action='store_true' : boolean switch: flag is not there: false , flag is there: true
    parser.add_argument('--extern', action='store_true', help="Scan via PROSITE-Website")

    args = parser.parse_args()


    seq=parse_fasta(args.fasta)
    current_pattern= None
    if args.pattern:
        current_pattern = args.pattern
    elif args.web:
        current_pattern = pattern_from_web(args.web)

    if current_pattern is None:
        print("No pattern", file=sys.stderr)
        sys.exit(1)

    regex_str = prosite_to_regex(current_pattern)
    scan_seq(seq, regex_str)


#translate prosite pattern to regex string
#prosite pattern is a AminoAcidSequence
def prosite_to_regex(pattern):

    regex = pattern.replace("<","^").replace(">","$")
    #{ABC} -> [^ABC]
    regex = re.sub(r'\{([A-Za-z]+)\}', r'[^\1]', regex)
    # 'x' --> '.'
    regex = regex.replace("x", ".")
    #remove -
    regex = regex.replace("-", "")
    # (n) or (n, m) -> {n} or {n, m}
    regex = re.sub(r'\((\d+(?:,\d+)?)\)', r'{\1}', regex)

    return regex


#scan seq ( of fasta ) for regex ( of translated prosite pattern)
#Lookahead-Anker (?=...) find all positions
#python starts counting with 0 , we start with 1
def scan_seq(seq, regex):
    search_regex = f"(?=({regex}))"
    #Lookahead trick, find all pos even though the last one is not over yet

    for header, seq in seq.items():
        for match in re.finditer(search_regex, seq):
            start_pos = match.start() + 1
            matched_substring = match.group(1)
            print(f"{header}\t{start_pos}\t{matched_substring}")


#find Prosite-Pattern from --web with prosite id
def pattern_from_web(prosite_id):
    url = "https://ftp.expasy.org/databases/prosite/prosite.dat"
    try:
        with urllib.request.urlopen(url) as response:
            found=False
            for line in response:
                line = line.decode('utf-8').rstrip('\n')
                if line.startswith("AC") and prosite_id in line: #AC = prosite id
                    found = True
                if found and line.startswith("PA"): #PA = pattern
                    return line[2:].strip().rstrip('.') #remove . and PA from pattern
    except Exception as e:
        print(f"Error fetching pattern: {e}", file=sys.stderr)
        sys.exit(1)
    return None

#example output from bioclient2:
#gi|15595239|ref|NP_248731.1| probable hemagglutinin [Pseudomonas aeruginosa]    451     LHVDTLGSA
#meaning of the example output: header + pos where seq that matches pattern starts  + the seq that matched pattern

if __name__ == "__main__":
    main()