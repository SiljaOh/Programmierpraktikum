#!/usr/bin/env python3
import argparse
import sys
import getpass

STOP = {"TAA", "TAG", "TGA"}

def read_fasta(path):
    header = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq).upper()
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq).upper()

def revcomp(s):
    comp = {"A":"T","T":"A","C":"G","G":"C"}
    return "".join(comp.get(b,"N") for b in reversed(s))

def find_orfs(seq):
    orfs = []
    n = len(seq)
    for frame in (0,1,2):
        start = None
        for i in range(frame, n-2, 3):
            codon = seq[i:i+3]
            if start is None:
                if codon == "ATG":
                    start = i
            else:
                if codon in STOP:
                    orfs.append(seq[start:i+3])
                    start = None
    return orfs




def connect_mysql(host, user, dbname, password):
    """
    Tries mysql-connector first, falls back to pymysql.
    """
    try:
        import mysql.connector  # type: ignore
        return mysql.connector.connect(
            host=host, user=user, password=password, database=dbname
        )
    except Exception:
        import pymysql  # type: ignore
        return pymysql.connect(
            host=host, user=user, password=password, database=dbname, autocommit=False
        )

def ensure_orf_table(cur):

    cur.execute("""
        CREATE TABLE IF NOT EXISTS orf (
            orf_id      INT AUTO_INCREMENT PRIMARY KEY,
            source      VARCHAR(60) NOT NULL,
            source_key  VARCHAR(120) NOT NULL,
            header      VARCHAR(255) NOT NULL,
            idx         INT NOT NULL,
            seq         LONGTEXT NOT NULL,
            length      INT NOT NULL,
            created_at  TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            UNIQUE KEY uniq_source_key (source, source_key),
            INDEX idx_header (header)
        ) DEFAULT CHARSET=utf8mb4
    """)

def insert_orf(cur, header, idx, orf_seq):
    source = "orf_finder"
    source_key = f"{header}_{idx}"
    cur.execute(
        """
        INSERT INTO orf (source, source_key, header, idx, seq, length)
        VALUES (%s, %s, %s, %s, %s, %s)
        ON DUPLICATE KEY UPDATE
            seq=VALUES(seq),
            length=VALUES(length)
        """,
        (source, source_key, header, idx, orf_seq, len(orf_seq))
    )



def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fasta", required=True)
    p.add_argument("--output")

    # DB flags (added)
    p.add_argument("--db", action="store_true", help="Insert ORFs into database table 'orf'")
    p.add_argument("--db-host", default="mysql2-ext.bio.ifi.lmu.de")
    p.add_argument("--db-user", default="bioprakt17")
    p.add_argument("--db-name", default="bioprakt17")
    p.add_argument("--db-pass")

    args = p.parse_args()

    out = open(args.output, "w") if args.output else sys.stdout

    conn = None
    cur = None

    if args.db:
        pw = args.db_pass if args.db_pass is not None else getpass.getpass("DB password: ")
        conn = connect_mysql(args.db_host, args.db_user, args.db_name, pw)
        cur = conn.cursor()
        ensure_orf_table(cur)

    for header, seq in read_fasta(args.fasta):
        idx = 0
        for orf in find_orfs(seq) + find_orfs(revcomp(seq)):
            out.write(">%s_%d\n%s\n" % (header, idx, orf))
            if cur is not None:
                insert_orf(cur, header, idx, orf)
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