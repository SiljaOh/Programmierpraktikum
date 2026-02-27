#!/usr/bin/env python3

import sys
import argparse
import psscan #use psscan.prosite_to_regex() , psscan.pattern_from_web(), psscan.parse_fasta()
import mysql.connector


def main():

    parser = argparse.ArgumentParser(description="database prosite pattern scanner")

    #fasta oder db
    group_mode = parser.add_mutually_exclusive_group(required=True)
    group_mode.add_argument('--db', action='store_true', help="Suche in DB")
    group_mode.add_argument('--fasta', type=str, help="Pfad zur FASTA Datei")

    #pattern oder web
    group_pat = parser.add_mutually_exclusive_group(required=True)
    group_pat.add_argument('--pattern', type=str, help="PROSITE Pattern")
    group_pat.add_argument('--web', type=str, help="PROSITE ID")

    args = parser.parse_args()

    # 1. Pattern holen from web  , nutze ppscan
    if args.web:
        current_pattern = psscan.pattern_from_web(args.web)
    else:
        current_pattern = args.pattern

    if not current_pattern:
        print("Error: No pattern found.", file=sys.stderr)
        sys.exit(1)

    # 2. Regex übersetzen (Nutzt psscan Modul)
    regex_str = psscan.prosite_to_regex(current_pattern)

    #3. DB suche, neu
    if args.db:
        try:
            conn = mysql.connector.connect(
                host="mysql2-ext.bio.ifi.lmu.de",
                user="bioprakt17",
                password="$1$18yyYqpW$y//O6M/6c6ke8EGazjJTx/",
                database="bioprakt17"
            )
            cursor = conn.cursor()

            # Nur Header ausgeben, SQL filtert die Sequenzen mit REGEXP! Platzhalter %s
            query = "SELECT description FROM sequence WHERE seq REGEXP %s"
            cursor.execute(query, (regex_str,))

            for (header,) in cursor:
                print(header)

            cursor.close()
            conn.close()
        except mysql.connector.Error as err:
            print(f"DB Error: {err}", file=sys.stderr)
            sys.exit(1)

    #oder 3. fasta einlesen
    elif args.fasta:
        if args.fasta == '-' :
            seq_dict = psscan.parse_fasta(sys.stdin)
        else:
            try:
                with open(args.fasta, "r") as fasta_file:
                    seq_dict = psscan.parse_fasta(fasta_file)
            except FileNotFoundError:
                print("Error: File not found", file=sys.stderr)
                sys.exit(1)

        psscan.scan_seq(seq_dict, regex_str)

if __name__ == "__main__":
    main()