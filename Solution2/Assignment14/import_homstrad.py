#!/usr/bin/env python3
import os, argparse, re
from getpass import getpass
import mysql.connector

#Parsen der .ali-Dateien
def parse_ali(path):
    with open(path, "r", errors="replace") as f:
        content = f.read()
    # Gruppe 1: Der Key direkt hinter >P1;
    # Gruppe 3: Die Sequenz bis zum *
    pattern = re.compile(r">P1;([^\r\n]+)[\r\n]+([^\r\n]+)[\r\n]+(.*?)\*", re.DOTALL)
    for match in pattern.finditer(content):
        final_key = match.group(1).strip() # PDB ID aus der Datei
        aligned = "".join(match.group(3).split()) # Sequenz ohne Leerzeichen
        if aligned:
            yield final_key, aligned

#Hilfsfunktion für einzelne DB-Queries
def one(cur, sql, params=()):
    cur.execute(sql, params)
    res = cur.fetchone()
    return res if res else None


def import_family(conn, base, fam):
    folder_path = os.path.join(base, fam)
    if not os.path.isdir(folder_path): return

    files = os.listdir(folder_path)
    final_path = None
    for f in files:
        if f.lower().endswith(".ali") or f.lower().endswith(".tem"):
            final_path = os.path.join(folder_path, f)
            break
    if not final_path: return

    with conn.cursor(buffered=True) as cur:

        cur.execute("INSERT IGNORE INTO family(name) VALUES(%s)", (fam,))
        fam_row = one(cur, "SELECT family_id FROM family WHERE name=%s", (fam,))
        fam_id = fam_row[0]

        cur.execute("INSERT IGNORE INTO alignment(name, type) VALUES(%s, 'msa')", (fam,))
        align_row = one(cur, "SELECT alignment_id FROM alignment WHERE name=%s AND type='msa'", (fam,))
        align_id = align_row[0]

        n = 0
        for key, aligned in parse_ali(final_path):
            aligned_db = aligned.replace(".", "-")
            ungapped = aligned_db.replace("-", "")

            cur.execute("SELECT sequence_id FROM sequence WHERE seq=%s", (ungapped,))
            seq_row = cur.fetchone()
            if seq_row:
                seq_id = seq_row[0]
            else:
                cur.execute("INSERT INTO sequence(seq, length) VALUES(%s,%s)", (ungapped, len(ungapped)))
                seq_id = cur.lastrowid

            cur.execute("INSERT IGNORE INTO source_entry(source, source_key, sequence_id) VALUES(%s,%s,%s)",
                        ("pdb", key, seq_id))

            cur.execute("INSERT IGNORE INTO family_member(family_id, source, source_key) VALUES(%s,%s,%s)",
                        (fam_id, "pdb", key))

            cur.execute("INSERT INTO alignment_row(alignment_id, sequence_id, aligned_seq) VALUES(%s,%s,%s)",
                        (align_id, seq_id, aligned_db))
            n += 1

        conn.commit()
        if n > 0:
            print(f"[OK] {fam}: {n} sequences")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--user", required=True)
    ap.add_argument("--homstrad_dir", required=True)
    ap.add_argument("--password", default=None)
    args = ap.parse_args()

    pw = args.password or os.environ.get("MYSQL_PW") or getpass("Password: ")

    try:
        conn = mysql.connector.connect(
            host=args.host,
            user=args.user,
            password=pw,
            database=args.db,
            autocommit=False
        )
    except Exception as e:
        print(f"Verbindungsfehler: {e}")
        return

    try:
        dirs = sorted([d for d in os.listdir(args.homstrad_dir)
                       if os.path.isdir(os.path.join(args.homstrad_dir, d))])
        for fam in dirs:
            import_family(conn, args.homstrad_dir, fam)
    finally:
        conn.close()

    print("\nFertig! Prüfe jetzt 'SELECT COUNT(*), COUNT(DISTINCT source_key) FROM family_member;'")


if __name__ == "__main__":
    main()






