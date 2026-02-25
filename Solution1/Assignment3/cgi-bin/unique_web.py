#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import defaultdict
import cgi
import html
import os
import tempfile
import traceback

# ----------------------------
# Konfiguration
# ----------------------------
MAX_UPLOAD_BYTES = 15 * 1024 * 1024   # 15 MB (soft limit)
MAX_K_COUNT = 50                      # max. wie viele k-Werte akzeptiert werden
MAX_K_VALUE = 1000000                 # Schutz gegen Quatschwerte (z.B. 10^9)

# ----------------------------
# Fachlogik (
# ----------------------------

def read_fasta(path):
    """Liest FASTA ein -> list[(header, seq)]. Wirft ValueError bei kaputter FASTA."""
    records, header, buf = [], None, []
    with open(path) as f:
        for line_no, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # Neue Sequenz beginnt
                if header is not None:
                    records.append((header, "".join(buf)))
                header = line[1:].strip()
                buf = []
            else:
                # Sequenzzeile
                if header is None:
                    raise ValueError(f"FASTA ungültig (Sequenz vor erstem Header) in Zeile {line_no}.")
                buf.append(line)

    if header is not None:
        records.append((header, "".join(buf)))

    if not records:
        raise ValueError("FASTA enthält keine Records (keine '>' Header gefunden).")
    return records

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

# ----------------------------
# Parsing & Validierung
# ----------------------------

def parse_k_list(k_raw: str):
    """Erlaubt '5 7 9' oder '5,7,9'. Wirft ValueError bei Ungültigem."""
    if k_raw is None:
        raise ValueError("Bitte k-Werte angeben.")
    s = k_raw.replace(",", " ").strip()
    if not s:
        raise ValueError("Bitte k-Werte angeben (z.B. '5 7 9').")

    toks = s.split()
    if len(toks) > MAX_K_COUNT:
        raise ValueError(f"Zu viele k-Werte ({len(toks)}). Maximal erlaubt: {MAX_K_COUNT}.")

    ks = []
    for t in toks:
        try:
            k = int(t)
        except ValueError:
            raise ValueError(f"Ungültiger k-Wert: '{t}' (muss Integer sein).")
        if k <= 0:
            raise ValueError(f"Ungültiger k-Wert: {k} (muss > 0 sein).")
        if k > MAX_K_VALUE:
            raise ValueError(f"k ist zu groß ({k}). Bitte <= {MAX_K_VALUE}.")
        ks.append(k)

    # Duplikate entfernen, aber Reihenfolge erhalten
    seen = set()
    out = []
    for k in ks:
        if k not in seen:
            seen.add(k)
            out.append(k)
    return out

def parse_start(start_raw: str):
    if start_raw is None or start_raw == "":
        return -1
    try:
        s = int(start_raw)
    except ValueError:
        raise ValueError(f"Ungültiger start-Wert: '{start_raw}' (muss Integer sein).")
    if s < -1:
        raise ValueError("start muss -1 oder >= 0 sein.")
    return s

def read_upload_to_temp(upload_item):
    """Speichert Upload in tempfile und gibt Pfad zurück. Wirft ValueError bei Problemen."""
    if upload_item is None or not getattr(upload_item, "file", None):
        raise ValueError("FASTA-Datei wurde nicht korrekt übertragen.")

    # Optional: Dateiname anzeigen/prüfen
    filename = getattr(upload_item, "filename", "") or ""
    # Soft extension check (nicht hart, weil manche FASTA keine Endung haben)
    if filename and not any(filename.lower().endswith(ext) for ext in [".fa", ".fasta", ".fna", ".ffn", ".txt", ".faa"]):
        # Nicht blocken, nur Hinweis (wir blocken nicht, weil sonst legitime Dateien nerven)
        pass

    data = upload_item.file.read()
    if data is None:
        raise ValueError("Upload konnte nicht gelesen werden.")
    if len(data) == 0:
        raise ValueError("Upload ist leer.")
    if len(data) > MAX_UPLOAD_BYTES:
        raise ValueError(f"Upload zu groß ({len(data)} Bytes). Max: {MAX_UPLOAD_BYTES} Bytes.")

    tmp = tempfile.NamedTemporaryFile(delete=False)
    tmp.write(data)
    tmp.close()
    return tmp.name, filename

# ----------------------------
# HTML Helper
# ----------------------------

def render_page(title, body_html):
    print("Content-Type: text/html; charset=utf-8\n")
    print(f"""<!doctype html>
<html lang="de">
<head>
  <meta charset="utf-8">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      --tum-blue:#1f4f82;
      --tum-blue-dark:#173b60;
      --bg:#ffffff;
      --card:#f6f7f9;
      --border:#e2e6ea;
      --text:#111;
      --muted:#555;
      --danger:#b00020;
    }}
    body {{
      margin: 0;
      font-family: Arial, Helvetica, sans-serif;
      background: var(--bg);
      color: var(--text);
    }}
    .wrap {{
      max-width: 1000px;
      margin: 28px auto;
      padding: 0 16px;
    }}
    h1 {{
      margin: 0 0 12px;
      color: var(--tum-blue-dark);
      font-size: 28px;
    }}
    .card {{
      background: var(--card);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 16px;
      margin: 12px 0;
    }}
    .hint {{
      color: var(--muted);
      margin: 6px 0 14px;
    }}
    label {{ font-weight: 700; display:block; margin-top: 10px; }}
    input[type=text], input[type=number] {{
      width: min(480px, 100%);
      padding: 8px;
      border-radius: 10px;
      border: 1px solid var(--border);
      margin-top: 6px;
      background: #fff;
    }}
    input[type=file] {{
      margin-top: 6px;
    }}
    button {{
      margin-top: 14px;
      padding: 10px 14px;
      border-radius: 10px;
      border: 0;
      background: var(--tum-blue);
      color: #fff;
      cursor: pointer;
      font-weight: 700;
    }}
    button:hover {{ filter: brightness(0.95); }}
    .danger {{
      color: var(--danger);
      font-weight: 700;
      margin: 0 0 10px;
    }}
    pre {{
      background: #fff;
      border: 1px solid var(--border);
      border-radius: 10px;
      padding: 12px;
      overflow: auto;
      margin: 10px 0 0;
    }}
    table {{
      border-collapse: collapse;
      width: 100%;
      background: #fff;
      border: 1px solid var(--border);
      border-radius: 10px;
      overflow: hidden;
    }}
    th, td {{
      border-bottom: 1px solid var(--border);
      padding: 10px 10px;
      text-align: left;
    }}
    th {{
      background: #f0f3f7;
      color: var(--tum-blue-dark);
    }}
    .toplinks a {{
      color: var(--tum-blue);
      text-decoration: none;
      font-weight: 700;
    }}
    .toplinks a:hover {{ text-decoration: underline; }}
  </style>
</head>
<body>
  <div class="wrap">
    {body_html}
  </div>
</body>
</html>""")

def form_html(default_k="5 7 9", default_start="-1", error_msg=None):
    err = f'<div class="card"><div class="danger">{html.escape(error_msg)}</div></div>' if error_msg else ""
    return f"""
    <div class="toplinks"><a href="../">← Zur Startseite</a></div>
    <h1>Aufgabe 12 – unique</h1>
    <div class="card">
      <div class="hint">
        Lade eine FASTA-Datei hoch. Gib k-Werte ein (z.B. <code>5 7 9</code>).
        <br><b>start = -1</b> → Teil a (alle k-mers), <b>start ≥ 0</b> → Teil b (fester Start).
      </div>

      <form method="post" enctype="multipart/form-data">
        <label>FASTA Datei</label>
        <input type="file" name="fasta" required>

        <label>k Werte</label>
        <input type="text" name="k" value="{html.escape(default_k)}" required>

        <label>start</label>
        <input type="number" name="start" value="{html.escape(default_start)}" required>

        <button type="submit">Ausführen</button>
      </form>
    </div>
    {err}
    """

# ----------------------------
# CGI Entry
# ----------------------------

try:
    form = cgi.FieldStorage()

    # Erstaufruf: Formular
    if "fasta" not in form:
        render_page("Aufgabe 12 – unique", form_html())
        raise SystemExit

    # Parameter lesen (mit Defaults für Re-Render)
    k_raw = form.getfirst("k", "")
    start_raw = form.getfirst("start", "-1")

    # Upload -> temp
    fasta_path = None
    fasta_filename = ""
    fasta_path, fasta_filename = read_upload_to_temp(form["fasta"])

    # parse/validate k & start
    ks = parse_k_list(k_raw)
    start = parse_start(start_raw)

    # FASTA lesen/validieren
    recs = read_fasta(fasta_path)
    n = len(recs)

    # k zusätzlich gegen Sequenzen checken (freundlicher Hinweis)
    max_len = max(len(seq) for _, seq in recs) if recs else 0
    if start == -1:
        too_big = [k for k in ks if k > max_len]
        if too_big:
            # nicht fatal, aber Hinweis
            pass
    else:
        # fester start: mindestens ein Gen muss lang genug sein
        if all(len(seq) < start + min(ks) for _, seq in recs):
            raise ValueError("Kein Gen ist lang genug für start + k. Bitte start/k anpassen.")

    # Ergebnisse berechnen
    rows = []
    for k in ks:
        u = unique_genes(recs, k, start)
        pct = (u / n * 100.0) if n else 0.0
        rows.append((k, u, pct))

    # HTML zusammenbauen
    table_rows = "\n".join(
        f"<tr><td>{k}</td><td>{u}</td><td>{pct:.2f}%</td></tr>"
        for (k, u, pct) in rows
    )
    tsv = "".join(f"{k}\t{u}\n" for (k, u, _) in rows)

    meta = f"Datei: <b>{html.escape(fasta_filename) if fasta_filename else '(Upload)'}</b> | Gene: <b>{n}</b> | start: <b>{start}</b>"
    body = f"""
      <div class="toplinks"><a href="../">← Zur Startseite</a></div>
      <h1>Aufgabe 12 – Ergebnis</h1>

      <div class="card">
        <div class="hint">{meta}</div>
        <table>
          <tr><th>k</th><th>#unique</th><th>% eindeutig</th></tr>
          {table_rows}
        </table>

        <h3 style="margin:14px 0 6px;color:var(--tum-blue-dark);">TSV Output</h3>
        <pre>{html.escape(tsv)}</pre>
      </div>

      <div class="card">
        <h3 style="margin:0 0 8px;color:var(--tum-blue-dark);">Neue Rechnung</h3>
        {form_html(default_k=k_raw or "5 7 9", default_start=str(start))}
      </div>
    """
    render_page("Aufgabe 12 – Ergebnis", body)

except ValueError as e:
    # „saubere“ Nutzerfehler: Formular erneut anzeigen + Fehlermeldung
    try:
        k_raw = locals().get("k_raw", "5 7 9")
        start_raw = locals().get("start_raw", "-1")
        render_page("Aufgabe 12 – unique (Fehler)", form_html(default_k=k_raw, default_start=start_raw, error_msg=str(e)))
    except Exception:
        print("Content-Type: text/plain; charset=utf-8\n")
        print("Fehler:", e)

except Exception:
    # unerwartete Fehler: Debug-Trace anzeigen (besser als 500 ohne Infos)
    print("Content-Type: text/html; charset=utf-8\n")
    print("<h1>Server-Fehler</h1>")
    print("<pre>" + html.escape(traceback.format_exc()) + "</pre>")

finally:
    # Tempfile aufräumen
    try:
        if "fasta_path" in locals() and locals()["fasta_path"]:
            os.unlink(locals()["fasta_path"])
    except Exception:
        pass