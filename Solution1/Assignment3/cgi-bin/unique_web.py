#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import cgi
import html
import os
import tempfile
import traceback

# --------- meine Logik (leicht angepasst) ---------

def read_fasta(path):
    listgenes, aktheader, bufferseq = [], None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if aktheader is not None:
                    listgenes.append((aktheader, "".join(bufferseq)))
                aktheader, bufferseq = line[1:], []
            else:
                bufferseq.append(line)
    if aktheader is not None:
        listgenes.append((aktheader, "".join(bufferseq)))
    return listgenes

def unique_genes(records, k, start):
    if k <= 0 or start < -1:
        return 0
    m = defaultdict(set)

    for gid, (_, seq) in enumerate(records):
        s = seq.upper()
        if start == -1:
            if len(s) < k:
                continue
            for i in range(len(s) - k + 1):
                m[s[i:i+k]].add(gid)
        else:
            if len(s) >= start + k:
                m[s[start:start+k]].add(gid)

    hit = [0] * len(records)
    for genes in m.values():
        if len(genes) == 1:
            hit[next(iter(genes))] = 1
    return sum(hit)

def parse_k_list(s):
    # erlaubt "5 7 9" oder "5,7,9"
    s = (s or "").replace(",", " ").strip()
    if not s:
        return []
    return [int(x) for x in s.split()]

# --------- CGI Teil ---------

print("Content-Type: text/html; charset=utf-8\n")

form = cgi.FieldStorage()

def page(title, body):
    print(f"""<!doctype html>
<html lang="de">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 24px; }}
    .box {{ max-width: 1000px; }}
    pre {{ background:#f6f6f6; padding:12px; border-radius:8px; overflow:auto; }}
    input[type=text], input[type=number] {{ padding: 6px; width: 340px; }}
    input[type=file] {{ margin: 6px 0 14px; }}
    button {{ padding: 8px 12px; }}
    table {{ border-collapse: collapse; margin-top: 12px; }}
    td, th {{ border: 1px solid #ddd; padding: 6px 10px; }}
  </style>
</head>
<body>
  <div class="box">
    {body}
  </div>
</body>
</html>""")

# Wenn noch keine Datei hochgeladen wurde -> Formular anzeigen
if "fasta" not in form:
    page(
        "Aufgabe 12 – unique",
        """
        <h1>Aufgabe 12 – unique</h1>
        <p>FASTA hochladen, k-Werte eingeben (z.B. <code>5 7 9</code>), start wählen.</p>

        <form method="post" enctype="multipart/form-data">
          <label><b>FASTA Datei</b></label><br>
          <input type="file" name="fasta" required><br>

          <label><b>k Werte</b> (z.B. 5 7 9 oder 5,7,9)</label><br>
          <input type="text" name="k" value="5 7 9" required><br><br>

          <label><b>start</b> (-1 = Teil A)</label><br>
          <input type="number" name="start" value="-1" required><br><br>

          <button type="submit">Ausführen</button>
        </form>
        """
    )
    raise SystemExit

# Upload + Parameter auslesen
try:
    fasta_item = form["fasta"]
    if not getattr(fasta_item, "file", None):
        raise ValueError("FASTA Upload wurde nicht korrekt übertragen.")

    ks = parse_k_list(form.getfirst("k", ""))
    start = int(form.getfirst("start", "-1"))
    if not ks:
        raise ValueError("Keine k-Werte angegeben.")
    if start < -1:
        raise ValueError("start muss -1 oder >= 0 sein.")

    # FASTA temporär speichern
    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(fasta_item.file.read())
    tmp.close()

    recs = read_fasta(tmp.name)
    n = len(recs)

    rows = []
    for k in ks:
        u = unique_genes(recs, k, start)
        pct = (u / n * 100.0) if n else 0.0
        rows.append((k, u, pct))

    tsv = "".join(f"{k}\t{u}\n" for (k, u, _) in rows)
    table_rows = "\n".join(f"<tr><td>{k}</td><td>{u}</td><td>{pct:.2f}%</td></tr>" for (k, u, pct) in rows)

    page(
        "Aufgabe 12 – Ergebnis",
        f"""
        <h1>Aufgabe 12 – Ergebnis</h1>
        <p><a href="./">Zurück zur Startseite</a></p>
        <p>Anzahl Gene: <b>{n}</b> | start: <b>{start}</b></p>

        <table>
          <tr><th>k</th><th>#unique</th><th>% eindeutig</th></tr>
          {table_rows}
        </table>

        <h3>TSV Output</h3>
        <pre>{html.escape(tsv)}</pre>
        """
    )

    os.unlink(tmp.name)

except Exception:
    page("Fehler", f"<h1>Fehler</h1><pre>{html.escape(traceback.format_exc())}</pre>")