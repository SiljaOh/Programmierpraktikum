#!/usr/bin/env python3

import os
import argparse
import mysql.connector

# 1. Argumente parsen
parser = argparse.ArgumentParser()
parser.add_argument('--input', required=True, help="Pfad zur SwissProt Datei")
args = parser.parse_args()

# 2. Verbindung herstellen
db = mysql.connector.connect(
    host="mysql2-ext.bio.ifi.lmu.de",
    user="bioprakt17",
    password=os.environ.get("DB_PASS", ""),
    database="bioprakt17"
)
cursor = db.cursor()

#LAUFZEIT VERBESSERUNG
# Wir speichern IDs im RAM, um SELECT-Anfragen zu vermeiden
#Idee:
org_cache = {}  # { "Homo sapiens": 1 }
#weniger SELECT Anfragen, da zb. viele Organismen Homo sapiens als OS haben
kw_cache = {}   # { "Glycolysis": 5 }
#gleiche idee für keywords


# 3. Einlesen von SwissProt-Datei:
#ID = Identifier (brauchen wir nicht)
#AC = Accession (in sequence)
#DE = Description  (in sequence)
#OS = Organism Species (in organism)
#KW = Keywords (in keyword)
#SQ = Die exakte Abfolge der Aminosäuren (in sequence)
#Jeder Eintrag endert mit //
# Idee: einen Eintrag einlesen --> commiten in db --> nächster Eintrag
#variablen danach reseten

with open(args.input, 'r') as file:
    accession = "" #string zb. P04406; UNIQUE
    description = "" #string zb. Glyceraldehyde-3-phosphate dehydrogenase OS=Homo sapiens...
    organism = "" #string zb. Homo sapiens (Human).
    keywords = [] #list! zb: [Glycolysis, Oxidoreductase, NAD]
    sequence_parts = []
    #list! zb: SEQUENCE   335 AA;  36053 MW;  D8A1E1E6E1F1F1F1 CRC64; header
     #MGKVKVGVNG FGRIGRLVTR AAFNSGKVDI VAINDPFIDL NYMVYMFQYD STHGKFHGTV
    #zu: [MGKVKVGVNG, FGRIGRLVTR , ...] am ende zu einer sequence joinen!

    in_sequence_block = False #schalter

    for line in file:
        line = line.strip() # Entfernt Leerzeichen am Anfang/Ende und das \n

        #PARSING LOGIK

        #sequenz einlesen:
        #1. wenn schalter AN ist
        if in_sequence_block and not line.startswith("//"):
            # Entferne alle Leerzeichen und Zahlen aus seq
            clean_seq = "".join([c for c in line if c.isalpha()])
            sequence_parts.append(clean_seq)
            continue


        # 2. Schalter AUS: suche nach AC, DE, OS, SQ, KW
        if line.startswith("AC"):
            accession = line.replace("AC", "").strip().split(";")[0]

        elif line.startswith("DE"):
            description += line.replace("DE", "").strip() + " "

        elif line.startswith("OS"):
            organism = line.replace("OS", "").strip().rstrip(".")

        elif line.startswith("SQ"):
            #schlater umlegen und die header zeile wird geskippt!
            in_sequence_block = True
            continue

        elif line.startswith("KW"):
            # Punkt am Ende weg, splitten und Leerzeichen weg
            kw_line = line.replace("KW", "").strip().rstrip(".")
            if kw_line: #KW darf nicht leer sein
                new_kws = [k.strip() for k in kw_line.split(";") if k.strip()]
                keywords.extend(new_kws)


        # 3. Wenn der Eintrag zu Ende ist (//)
        elif line.startswith("//"):

            if not accession:  # Falls die Datei mit // startet oder ein leerer Eintrag kommt
                continue

            final_seq = "".join(sequence_parts) #finales joinen
            length = len(final_seq) #für length in sequence

            #UPLOAD LOGIK mit LAUFZEITVERBESSERUNG
            #duplikate verhindern: sequence nur einfügen wenn accession noch nich exsistiert
            """
            OS: 
            eingelesener organism wird zu name von organism
            organism_id von organism = organism_id von sequence
            """

            if organism not in org_cache:
                # Da kein Foreign Key existiert: Erst schauen, ob er da ist
                cursor.execute("SELECT organism_id FROM organism WHERE name = %s", (organism,))
                res = cursor.fetchone()
                if res:
                    org_cache[organism] = res[0]
                else:
                    cursor.execute("INSERT INTO organism (name) VALUES (%s)", (organism,))
                    org_cache[organism] = cursor.lastrowid

            curr_org_id = org_cache[organism]

            """
            SQ:
            eingelesene sequence wird zu seq von sequence
            length wird zu length von sequence
            
            AC:
            eingelesene accession wird zu accession in sequence
            
            DE:
            eigelesene description wird zu description in sequence
            """
            sql_seq = """INSERT IGNORE INTO sequence (accession, description, seq, length, organism_id) VALUES (%s, %s, %s, %s, %s)"""
            cursor.execute(sql_seq, (accession, description.strip(), final_seq, length, curr_org_id))

            if cursor.lastrowid != 0:
                curr_seq_id = cursor.lastrowid
            else:
                # Falls Accession schon da war ID kurz abrufen
                cursor.execute("SELECT sequence_id FROM sequence WHERE accession = %s", (accession,))
                curr_seq_id = cursor.fetchone()[0]

            """
            KW:
            eingelesene keywords wird zu word von keyword
            -->keyword_id von keyword wird zu keyword_id von sequence_keyword
            -->sequence_id von sequence_keyword wird zu sequence_id von sequence
            """
            for kw in keywords:
                if kw not in kw_cache:
                    cursor.execute("SELECT keyword_id FROM keyword WHERE word = %s", (kw,))
                    kw_res = cursor.fetchone()
                    if kw_res:
                        kw_cache[kw] = kw_res[0]
                    else:
                        cursor.execute("INSERT INTO keyword (word) VALUES (%s)", (kw,))
                        kw_cache[kw] = cursor.lastrowid

                curr_kw_id = kw_cache[kw]

                # Verknüpfung INSERT IGNORE ist schneller als SELECT
                cursor.execute("INSERT IGNORE INTO sequence_keyword (sequence_id, keyword_id) VALUES (%s, %s)",
                               (curr_seq_id, curr_kw_id))

            #Reset
            accession = ""
            description = ""
            organism = ""
            keywords = []
            sequence_parts = []
            in_sequence_block = False
            continue

db.commit()
cursor.close()
db.close()