#!/usr/bin/env python3
# Importiert HOMSTRAD .ali Dateien in unser MySQL/MariaDB Datenbank.

import os, argparse #os fuer Pfade/ORdner
from getpass import getpass #passworteingabe
import mysql.connector #conne


def parse_ali(path):
    """Liest eine .ali Datei und gibt Einträge (pdb, chain, aligned_seq) zurück."""
    pdb = chain = None
    seq = []
    in_seq = False

    for line in open(path, "r", errors="replace"):
        line = line.strip()
        if not line or line.startswith("C;"):   # Kommentare/leere Zeilen überspringen
            continue

        if line.startswith(">P1;"):            # neuer Eintrag beginnt
            pdb = chain = None
            seq = []
            in_seq = False
            continue

        if line.lower().startswith("structurex:"):   # hier stehen PDB + Chain
            parts = line.split(":")
            pdb = (parts[1].strip().lower() if len(parts) > 1 else "")
            chain = (parts[3].strip() if len(parts) > 3 and parts[3].strip() else "-")
            in_seq = True
            continue

        if in_seq:                              # Sequenzzeilen sammeln
            line = line.replace(" ", "")
            if "*" in line:                     # '*' beendet die Sequenz
                seq.append(line.split("*", 1)[0])
                aligned = "".join(seq)
                if pdb and aligned:
                    yield pdb, chain, aligned
                in_seq = False
            else:
                seq.append(line)


def one(cur, sql, params=()):
    """Führt SQL aus und gibt genau eine Zeile (oder None) zurück."""
    cur.execute(sql, params)
    return cur.fetchone()


def import_family(conn, base, fam):
    ali = os.path.join(base, fam, fam + ".ali")
    if not os.path.isfile(ali):
        print("[SKIP]", fam, "kein .ali")
        return

    with conn.cursor(buffered=True) as cur:
        # 1) family holen/erzeugen
        cur.execute("INSERT IGNORE INTO family(name) VALUES(%s)", (fam,))
        fam_id = one(cur, "SELECT family_id FROM family WHERE name=%s", (fam,))[0]

        # 2) alignment holen/erzeugen
        cur.execute("INSERT IGNORE INTO alignment(name,type) VALUES(%s,'msa')", (fam,))
        align_id = one(cur, "SELECT alignment_id FROM alignment WHERE name=%s AND type='msa'", (fam,))[0]

        n = 0
        for pdb, chain, aligned in parse_ali(ali):
            source = "pdb"
            key = f"{pdb}_{chain}"
            ungapped = aligned.replace("-", "")  # echte Sequenz ohne '-'

            # 3) sequence_id finden oder neu anlegen
            row = one(cur, "SELECT sequence_id FROM source_entry WHERE source=%s AND source_key=%s",
                      (source, key))
            if row:
                seq_id = row[0]
            else:
                cur.execute("INSERT INTO sequence(seq,length,organism_id) VALUES(%s,%s,NULL)",
                            (ungapped, len(ungapped)))
                seq_id = cur.lastrowid
                cur.execute("INSERT IGNORE INTO source_entry(source,source_key,sequence_id) VALUES(%s,%s,%s)",
                            (source, key, seq_id))

            # 4) Mitgliedschaft speichern
            cur.execute("INSERT IGNORE INTO family_member(family_id,source,source_key) VALUES(%s,%s,%s)",
                        (fam_id, source, key))

            # 5) alignment row speichern (aligned_seq mit '-')
            cur.execute("INSERT IGNORE INTO alignment_row(alignment_id,sequence_id,aligned_seq) VALUES(%s,%s,%s)",
                        (align_id, seq_id, aligned))

            n += 1

        conn.commit()
        print("[OK] Imported", fam, "with", n, "sequences")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--user", required=True)
    ap.add_argument("--homstrad_dir", required=True)
    ap.add_argument("--family", default=None)  # wenn leer: alle Familien
    ap.add_argument("--password", default=None)
    args = ap.parse_args()

    pw = args.password or os.environ.get("MYSQL_PW") or getpass("MySQL password: ")

    conn = mysql.connector.connect(host=args.host, user=args.user, password=pw, database=args.db)
    try:
        if args.family:
            import_family(conn, args.homstrad_dir, args.family)
        else:
            for fam in sorted(d for d in os.listdir(args.homstrad_dir)
                              if os.path.isdir(os.path.join(args.homstrad_dir, d))):
                import_family(conn, args.homstrad_dir, fam)
    finally:
        conn.close()


if __name__ == "__main__":
    main()