#!/usr/bin/env python3
import argparse, sys, getpass

STOP = {"TAA", "TAG", "TGA"}

def read_fasta(path):
    h, s = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if h is not None:
                    yield h, "".join(s).upper()
                h, s = line[1:].strip(), []
            else:
                s.append(line)
        if h is not None:
            yield h, "".join(s).upper()

def revcomp(s):
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(comp.get(b, "N") for b in reversed(s.upper()))

def find_orfs(seq):
    res, n = [], len(seq)
    for frame in (0,1,2):
        start = None
        for i in range(frame, n-2, 3):
            codon = seq[i:i+3]
            if start is None:
                if codon == "ATG":
                    start = i
            else:
                if codon in STOP:
                    res.append((seq[start:i+3], frame, start, i+2))
                    start = None
    return res

def connect_mysql(host, user, dbname, password):
    try:
        import mysql.connector
        return mysql.connector.connect(host=host, user=user, password=password, database=dbname)
    except Exception:
        import pymysql
        return pymysql.connect(host=host, user=user, password=password, database=dbname, autocommit=False)

def get_or_create_organism_id(cur, name):
    cur.execute("SELECT organism_id FROM organism WHERE name=%s", (name,))
    row = cur.fetchone()
    if row:
        return row[0]
    cur.execute("INSERT INTO organism (name) VALUES (%s)", (name,))
    return cur.lastrowid

def insert_orf(cur, header, idx, strand, frame, start, end, seq, organism_id):
    desc = "ORF strand=%s frame=%d start=%d end=%d" % (strand, frame, start, end)
    cur.execute(
        "INSERT INTO sequence (accession, description, seq, length, organism_id) VALUES (%s,%s,%s,%s,%s)",
        (None, desc, seq, len(seq), organism_id)
    )
    sid = cur.lastrowid
    cur.execute(
        "INSERT INTO source_entry (source, source_key, sequence_id) VALUES (%s,%s,%s)",
        ("orf_finder", "%s_%d" % (header, idx), sid)
    )

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--output")
    p.add_argument("--db", action="store_true")
    p.add_argument("--db-host", default="mysql2-ext.bio.ifi.lmu.de")
    p.add_argument("--db-user", default="bioprakt17")
    p.add_argument("--db-name", default="bioprakt17")
    p.add_argument("--db-pass")
    p.add_argument("--organism")
    args = p.parse_args()

    out = open(args.output, "w") if args.output else sys.stdout

    conn = None
    cur = None
    organism_id = None

    if args.db:
        pw = args.db_pass if args.db_pass is not None else getpass.getpass("DB password: ")
        conn = connect_mysql(args.db_host, args.db_user, args.db_name, pw)
        cur = conn.cursor()
        if args.organism:
            organism_id = get_or_create_organism_id(cur, args.organism)

    for header, seq in read_fasta(args.fasta):
        fwd = find_orfs(seq)
        rev = find_orfs(revcomp(seq))

        idx = 0
        for (orfseq, frame, start, end) in fwd:
            out.write(">%s_%d\n%s\n" % (header, idx, orfseq))
            if cur is not None:
                insert_orf(cur, header, idx, "+", frame, start, end, orfseq, organism_id)
            idx += 1

        for (orfseq, frame, start, end) in rev:
            out.write(">%s_%d\n%s\n" % (header, idx, orfseq))
            if cur is not None:
                insert_orf(cur, header, idx, "-", frame, start, end, orfseq, organism_id)
            idx += 1

        if conn is not None:
            conn.commit()

    if args.output:
        out.close()
    if cur is not None:
        cur.close()
    if conn is not None:
        conn.close()

if __name__ == "__main__":
    main()



