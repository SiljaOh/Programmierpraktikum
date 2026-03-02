#!/usr/bin/env python3
import sys
#Überlappung : start neuesGEn< ende vorheriges gen
if len(sys.argv) != 2: #Programm erwartet geneau 1x Datei als argument
    sys.exit("Usage: python overlaps.py saccharomyces.gff")

groups = {} #hier Gene sortiert speichern(Chromosom, Strand)

for line in open(sys.argv[1]): #Datei einlesen
    if line.startswith("#"): #info zeilen ignorieren
        continue
    p = line.rstrip().split("\t") #9 SPalten getrennt durch tabs
    if len(p) != 9 or p[2] != "gene" or p[6] not in "+-": # nur vollständige Zeilen, 'gene', strand +/-
        continue
    seqid, start, end, strand, attrs = p[0], int(p[3]), int(p[4]), p[6], p[8] # wichitge WErte holen
    gid = next((x.split("=")[1] for x in attrs.split(";") if x.startswith("ID=")), "?") #gen ID holen
    groups.setdefault((seqid, strand), []).append((start, end, gid)) #gene gruppieren

count = 0

for (seqid, strand), genes in groups.items(): #überlappung suchen
    genes.sort() #Gene werden nach start pos soritert
    ps, pe, pid = genes[0] #erstes Gen markieren
    for s, e, gid in genes[1:]: #jedes neue Gen mit vorherigem vergleichen
        if s <= pe: #neuer start leigt vor altem ende
            count += 1
            print(seqid, strand, f"{pid}:{ps}-{pe}", f"{gid}:{s}-{e}", sep="\t") #ausgabe
            pe = max(pe, e)
        else:
            ps, pe, pid = s, e, gid #kein overlap-> neues referenz GEn setzen

print(f"#Overlaps_same_strand={count}", file=sys.stderr) #totale anzahl der Overlaps ausgeben