#!/usr/bin/env python3

import argparse
import getpass
import pymysql
import pymysql.cursors

# get_seq.py -> holt Sequenzen aus der DB und gibt FASTA aus
# DB Struktur (MariaDB) Gruppe 17:
#   source_entry(source, source_key, sequence_id)
#   sequence(sequence_id, seq, accession, ...)

# ArgumentParser erstellen
parser = argparse.ArgumentParser(description="Fetch sequences from DB and print FASTA")

# metavar: Anzeige-Name des Arguments im Usage/Help (nur Darstellung, keine Auswirkung auf die Programmlogik)
parser.add_argument('-i', '--id', dest='id', nargs="+", required=True, metavar='ids',
                    help='Input Id(s) (source_entry.source_key)')
parser.add_argument('-s', '--source', dest='source', required=True, metavar='source',
                    help='Input Source (source_entry.source)')

args = parser.parse_args()

# Passwort eingeben (nicht im Code speichern)
dbpass = getpass.getpass("DB password: ")

# Verbindung an DB
connection = pymysql.connect(
    host="mysql2-ext.bio.ifi.lmu.de",
    port=3306,
    user="bioprakt17",
    password=dbpass,
    database="bioprakt17",
    cursorclass=pymysql.cursors.DictCursor,
)

# SQL
# JOIN: source_entry wird mit sequence über sequence_id verknüpft
#   se.source = gewünschte Quelle
#   se.source_key IN = Liste der übergebenen IDs
# ORDER BY -> reproduzierbare/geordnete Ausgabe
placeholders = ",".join(["%s"] * len(args.id))
sql = f"""
    SELECT se.source_key, s.accession, s.seq
    FROM source_entry se
    JOIN sequence s ON s.sequence_id = se.sequence_id
    WHERE se.source = %s
      AND se.source_key IN ({placeholders})
    ORDER BY se.source_key
"""

# Parameterliste für die vorbereitete SQL-Abfrage
# erstes Element: source
# danach: alle übergebenen IDs
params = [args.source] + args.id

# SQL-Query ausführen
# parametrisierte Ausführung: SQL-Struktur und Werte werden getrennt übergeben
# verhindert Manipulation der Abfrage (SQL-Injection)
try:
    with connection.cursor() as cur:
        cur.execute(sql, params)        # alle Treffer als Liste von Dictionaries
        rows = cur.fetchall()
finally:
    connection.close()                    # DB-Verbindung schließen (immer)

# FASTA Ausgabe
for row in rows:
    header = row["accession"] if row["accession"] else row["source_key"]
    print(f">{header}")
    print(row["seq"])

# __name__ == "__main__" nur beim direkten Start der Datei
# d.h. erhindert Ausführung beim Import als Modul
if __name__ == "__main__":
    main()