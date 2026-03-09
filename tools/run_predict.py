#!/usr/bin/python3

import subprocess
import html
import base64
import os
import cgi
import tempfile
import urllib.request

CGI_DIR = os.path.dirname(os.path.abspath(__file__))
JAR_PATH          = os.path.join(CGI_DIR, "predict.jar")
VALIDATE_JAR_PATH = os.path.join(CGI_DIR, "validateGor.jar")


def run_gor_predict(model_path, seq_path, fmt, probabilities):
    cmd = [
        "java", "-cp", JAR_PATH, "de.lmu.ifi.bio.gor.cli.Predict",
        "--model", model_path,
        "--seq",   seq_path,
        "--format", fmt,
    ]
    if probabilities:
        cmd.append("--probabilities")
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=300)
    except Exception as e:
        return e


def run_validate(pred_path, ref_path, summary_path, detailed_path):
    cmd = [
        "java", "-jar", VALIDATE_JAR_PATH,
        "-p", pred_path,
        "-r", ref_path,
        "-f", "txt",
        "-s", summary_path,
        "-d", detailed_path,
    ]
    try:
        return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              universal_newlines=True, timeout=300)
    except Exception as e:
        return e


def fetch_fasta_from_pdb(pdb_id):
    pdb_id = pdb_id.strip()
    if not pdb_id or len(pdb_id) != 4 or not pdb_id.isalnum():
        return None, f"Ungültige PDB ID: '{pdb_id}'."
    tmp = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode='w')
    tmp.close()
    try:
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}/display"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            content = resp.read().decode('utf-8')
            with open(tmp.name, 'w') as f:
                f.write(content)
        if os.path.getsize(tmp.name) < 10:
            return None, f"PDB ID '{pdb_id.upper()}' lieferte keine FASTA-Daten."
        return tmp.name, None
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code} für PDB '{pdb_id.upper()}'."
    except Exception as e:
        return None, str(e)


def parse_summary_scores(summary_text):
    """Liest Mean-Werte aus Summary-Text."""
    scores = {}
    for line in summary_text.splitlines():
        for key in ("q3", "SOV", "qObs_H", "qObs_E", "qObs_C", "SOV_H", "SOV_E", "SOV_C"):
            stripped = line.strip()
            if stripped.startswith(key + ":") or stripped.startswith(key + " "):
                parts = line.split("Mean:")
                if len(parts) > 1:
                    try:
                        scores[key] = float(parts[1].strip().split()[0])
                    except:
                        pass
    return scores


#html

HTML_HEAD = """<!DOCTYPE html>
<html lang="de">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GOR Sekundärstruktur-Vorhersage</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <script src="https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js"></script>
    <style>
        :root {
            --bio-green: #0d6b5c; --bio-green-light: #d4ede8; --bio-green-pale: #e8f6f3;
            --bio-teal: #2a9d8f; --bio-dark: #0f2926; --bio-deep: #0a1f1c;
            --text: #1a2e2a; --text-muted: #4a6560; --card-bg: #ffffff; --border: #c5e0db;
        }
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body { font-family: "DM Sans", system-ui, sans-serif; background: var(--bio-deep);
               color: var(--text); line-height: 1.6; min-height: 100vh; }
        body::before {
            content: ""; position: fixed; inset: 0;
            background:
                linear-gradient(180deg, rgba(13,107,92,0.06) 0%, transparent 40%),
                radial-gradient(ellipse 80% 50% at 50% -20%, rgba(42,157,143,0.15), transparent),
                repeating-linear-gradient(0deg, transparent, transparent 32px, rgba(13,107,92,0.03) 32px, rgba(13,107,92,0.03) 33px),
                repeating-linear-gradient(90deg, transparent, transparent 32px, rgba(13,107,92,0.03) 32px, rgba(13,107,92,0.03) 33px);
            pointer-events: none; z-index: 0;
        }
        .header {
            position: relative;
            background: linear-gradient(160deg, var(--bio-deep) 0%, var(--bio-dark) 40%, var(--bio-green) 100%);
            color: #fff; padding: 40px 24px 48px; text-align: center; overflow: hidden;
        }
        .header::after { content:""; position:absolute; bottom:0; left:0; right:0; height:60px;
            background: linear-gradient(to top, var(--bio-deep), transparent); pointer-events:none; }
        .header-inner { position: relative; z-index: 1; }
        .header h1 { margin: 0 0 6px; font-size: 1.6rem; font-weight: 700; letter-spacing: 0.02em; }
        .header p  { margin: 0; opacity: 0.85; font-size: 0.95rem; }
        .hero-dots { display:flex; gap:8px; justify-content:center; margin-top:14px; }
        .hero-dots span { width:8px; height:8px; border-radius:50%; background:rgba(255,255,255,0.45); }
        .main { position:relative; z-index:1; max-width:900px; margin:-24px auto 0; padding:0 20px 56px; }
        .main-panel {
            background: linear-gradient(180deg, rgba(255,255,255,0.97) 0%, rgba(248,252,251,0.98) 100%);
            border-radius:20px; padding:36px 32px 40px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.25), 0 0 0 1px rgba(255,255,255,0.1) inset;
        }
        .card { background:var(--card-bg); border:1px solid var(--border); border-radius:14px;
                padding:24px 24px 20px; margin-bottom:20px; }
        .card-title { font-size:0.78em; font-weight:700; text-transform:uppercase;
                      letter-spacing:0.08em; color:var(--bio-green); margin-bottom:16px; }
        .tab-buttons { display:flex; gap:8px; margin-bottom:18px; }
        .tab-btn {
            flex:1; padding:9px 12px; border:2px solid var(--border); border-radius:10px;
            background:white; cursor:pointer; font-weight:600; font-size:0.88em;
            color:var(--bio-green); font-family:"DM Sans",sans-serif; transition:all 0.2s;
        }
        .tab-btn:hover  { border-color:var(--bio-teal); background:var(--bio-green-pale); }
        .tab-btn.active { border-color:var(--bio-green); background:var(--bio-green); color:white; }
        .tab-content { display:none; } .tab-content.active { display:block; }
        label { display:block; font-weight:600; font-size:0.9em; color:var(--text); margin-bottom:7px; }
        select, input[type="text"] {
            width:100%; padding:10px 13px; border:1.5px solid var(--border); border-radius:10px;
            font-size:0.93em; font-family:"DM Sans",sans-serif; color:var(--text);
            background:var(--bio-green-pale); outline:none; transition:border-color 0.2s;
        }
        select:focus, input[type="text"]:focus { border-color:var(--bio-teal); background:white; }
        input[type="file"] {
            width:100%; padding:10px 13px; border:1.5px dashed var(--border); border-radius:10px;
            font-size:0.88em; font-family:"DM Sans",sans-serif; background:var(--bio-green-pale); cursor:pointer;
        }
        textarea {
            width:100%; padding:10px 13px; border:1.5px solid var(--border); border-radius:10px;
            font-size:0.88em; font-family:"JetBrains Mono",monospace; color:var(--text);
            background:var(--bio-green-pale); outline:none; resize:vertical; transition:border-color 0.2s;
        }
        textarea:focus { border-color:var(--bio-teal); background:white; }
        small { display:block; margin-top:6px; font-size:0.82em; color:var(--text-muted); }
        .checkbox-row { display:flex; align-items:center; gap:10px; }
        .checkbox-row input[type="checkbox"] { width:17px; height:17px; accent-color:var(--bio-green); cursor:pointer; }
        .checkbox-row label { margin:0; font-weight:500; cursor:pointer; font-size:0.93em; }
        .submit-btn {
            width:100%; padding:14px;
            background:linear-gradient(135deg, var(--bio-dark), var(--bio-green));
            color:white; border:none; border-radius:12px; font-size:1em; font-weight:700;
            font-family:"DM Sans",sans-serif; cursor:pointer; letter-spacing:0.03em;
            transition:opacity 0.2s; margin-top:5px;
        }
        .submit-btn:hover { opacity:0.88; }
        .status-msg { padding:13px 17px; border-radius:12px; margin:16px 0; font-size:0.93em; border-left:4px solid; }
        .error   { background:#fdf0ef; color:#8b2a22; border-color:#c0392b; }
        .success { background:var(--bio-green-pale); color:var(--bio-dark); border-color:var(--bio-green); }
        .info    { background:#eef5fb; color:#1a4d7a; border-color:#2980b9; }
        .result-html-box {
            background:#ffffff; color:#1a2e2a; padding:24px 28px; border-radius:14px;
            font-family:"JetBrains Mono",monospace; font-size:0.95em; max-height:600px;
            overflow-y:auto; overflow-x:auto; margin-top:14px; border:1.5px solid var(--border);
            line-height:1.8; white-space:pre; box-shadow:0 2px 12px rgba(13,107,92,0.08);
        }
        .result-title { font-weight:700; font-size:0.9em; color:var(--bio-green);
                        margin:20px 0 4px; text-transform:uppercase; letter-spacing:0.06em; }
        .download-btn {
            display:inline-block; background:linear-gradient(135deg, var(--bio-dark), var(--bio-green));
            color:white; padding:11px 22px; text-decoration:none; border-radius:10px;
            font-weight:700; font-size:0.92em; margin-top:14px; margin-right:8px; transition:opacity 0.2s;
        }
        .download-btn:hover { opacity:0.88; }
        /* Score table */
        .score-table { width:100%; border-collapse:collapse; margin-top:14px;
                       font-family:"JetBrains Mono",monospace; font-size:0.9em; }
        .score-table th { background:var(--bio-green); color:white; padding:10px 16px;
                          text-align:left; font-size:0.82em; letter-spacing:0.05em; }
        .score-table td { padding:9px 16px; border-bottom:1px solid var(--bio-green-light); }
        .score-table tr:last-child td { border-bottom:none; }
        .score-table tr:hover td { background:var(--bio-green-pale); }
        .score-good { color:#0d6b5c; font-weight:700; }
        .score-mid  { color:#e07b00; font-weight:600; }
        .score-low  { color:#c0392b; font-weight:600; }
        /* Jmol */
        .jmol-wrap { background:#fff; border:1.5px solid var(--border); border-radius:14px;
                     padding:20px; margin-top:14px; box-shadow:0 2px 12px rgba(13,107,92,0.08); }
        .jmol-wrap h4 { font-size:0.82em; font-weight:700; text-transform:uppercase;
                        letter-spacing:0.07em; color:var(--bio-green); margin-bottom:12px; }
        .jmol-legend { display:flex; gap:14px; flex-wrap:wrap; margin-top:12px; font-size:0.82em; color:var(--text-muted); }
        .jmol-legend span { display:flex; align-items:center; gap:6px; }
        .leg-dot { width:11px; height:11px; border-radius:50%; display:inline-block; }
        footer { position:relative; z-index:1;
            background:linear-gradient(180deg, var(--bio-dark) 0%, var(--bio-deep) 100%);
            color:rgba(255,255,255,0.88); padding:22px 24px; font-size:0.88rem;
            text-align:center; border-top:1px solid rgba(13,107,92,0.4); margin-top:40px;
        }
        footer a { color:#9ae6d9; text-decoration:none; }
        footer a:hover { text-decoration:underline; }
    </style>
    <script>
        function switchTab(tabId) {
            document.querySelectorAll('.seq-tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.seq-tab-content').forEach(c => c.classList.remove('active'));
            document.getElementById('tab-' + tabId).classList.add('active');
            document.getElementById('content-' + tabId).classList.add('active');
            document.getElementById('seq_mode').value = tabId;
        }
        function switchRefTab(tabId) {
            document.querySelectorAll('.ref-tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.ref-tab-content').forEach(c => c.classList.remove('active'));
            document.getElementById('ref-tab-' + tabId).classList.add('active');
            document.getElementById('ref-content-' + tabId).classList.add('active');
            document.getElementById('ref_mode').value = tabId;
        }
    </script>
</head>
<body>
<header class="header">
    <div class="header-inner">
        <h1>GOR Sekundärstruktur-Vorhersage</h1>
        <p>WS 2025/2026 <span>&middot;</span> Gruppe 17</p>
        <div class="hero-dots"><span></span><span></span><span></span></div>
    </div>
</header>
<main class="main">
<div class="main-panel">
"""

HTML_FORM = """
<form method="post" enctype="multipart/form-data">
    <input type="hidden" name="seq_mode" id="seq_mode" value="text">
    <input type="hidden" name="ref_mode" id="ref_mode" value="none">

    <!-- 1. Sequenz -->
    <div class="card">
        <div class="card-title">1. Sequenz eingeben</div>
        <div class="tab-buttons">
            <button type="button" class="seq-tab-btn tab-btn active" id="tab-text"
                    onclick="switchTab('text')">Sequenz eingeben</button>
            <button type="button" class="seq-tab-btn tab-btn" id="tab-pdbid"
                    onclick="switchTab('pdbid')">PDB ID</button>
            <button type="button" class="seq-tab-btn tab-btn" id="tab-upload"
                    onclick="switchTab('upload')">FASTA hochladen</button>
        </div>
        <div class="seq-tab-content tab-content active" id="content-text">
            <label>Aminosäuresequenz (FASTA-Format)</label>
            <textarea name="seq_text" rows="5" placeholder=">MeinProtein&#10;ACDEFGHIKLMNPQRSTVWY..."></textarea>
            <small>Erste Zeile mit &gt;ID, dann Sequenz</small>
        </div>
        <div class="seq-tab-content tab-content" id="content-pdbid">
            <label>PDB ID (z.B. 1TIM)</label>
            <input type="text" name="pdb_id" placeholder="z.B. 1TIM, 2HHB...">
            <small>Sequenz + 3D-Struktur werden automatisch von rcsb.org geladen.</small>
        </div>
        <div class="seq-tab-content tab-content" id="content-upload">
            <label>FASTA-Datei hochladen</label>
            <input type="file" name="seq_file" accept=".fasta,.fa,.txt">
            <small>Erlaubte Formate: .fasta, .fa</small>
        </div>
    </div>

    <!-- 2. Modell -->
    <div class="card">
        <div class="card-title">2. Modell hochladen</div>
        <label>Trainiertes GOR-Modell</label>
        <input type="file" name="model_file" accept=".txt,.model">
        <small>Modell-Datei aus dem Training-Tool (GOR I, III oder IV)</small>
    </div>

    <!-- 3. Referenzstruktur (optional) -->
    <div class="card">
        <div class="card-title">3. Referenzstruktur (optional) &ndash; für Q3 &amp; SOV Scores</div>
        <div class="tab-buttons">
            <button type="button" class="ref-tab-btn tab-btn active" id="ref-tab-none"
                    onclick="switchRefTab('none')">Keine</button>
            <button type="button" class="ref-tab-btn tab-btn" id="ref-tab-upload"
                    onclick="switchRefTab('upload')">Datei hochladen</button>
            <button type="button" class="ref-tab-btn tab-btn" id="ref-tab-text"
                    onclick="switchRefTab('text')">Text einfügen</button>
        </div>
        <div class="ref-tab-content tab-content active" id="ref-content-none">
            <small style="margin-top:0">Kein Vergleich – nur Vorhersage anzeigen.</small>
        </div>
        <div class="ref-tab-content tab-content" id="ref-content-upload">
            <label>Referenz-Datei (.seclib)</label>
            <input type="file" name="ref_file" accept=".seclib,.txt">
            <small>Bekannte Sekundärstruktur im seclib-Format</small>
        </div>
        <div class="ref-tab-content tab-content" id="ref-content-text">
            <label>Referenz einfügen</label>
            <textarea name="ref_text" rows="5" placeholder="&gt; ProteinID&#10;AS ACDEFG...&#10;SS HHHCCC..."></textarea>
        </div>
    </div>

    <!-- 4. Optionen -->
    <div class="card">
        <div class="card-title">4. Optionen</div>
        <div class="checkbox-row" style="margin-bottom:12px">
            <input type="checkbox" name="probabilities" id="probabilities" value="1">
            <label for="probabilities">Wahrscheinlichkeiten ausgeben (PH, PE, PC)</label>
        </div>
        <div class="checkbox-row" style="margin-bottom:12px">
            <input type="checkbox" name="show_jmol" id="show_jmol" value="1" checked>
            <label for="show_jmol">3D-Struktur mit JSmol anzeigen (nur bei PDB ID)</label>
        </div>
    </div>

    <button type="submit" name="predict_action" class="submit-btn">Vorhersage starten</button>
</form>
"""

HTML_FOOTER = """
</div>
</main>
<footer>
    ProgPrakt WS 2025/2026 &middot; Gruppe 17 &nbsp;|&nbsp;
    <a href="../index.html">&larr; Zurück zur Hauptseite</a>
</footer>
</body>
</html>
"""


#Main

def main():
    fs = cgi.FieldStorage()
    print("Content-Type: text/html; charset=utf-8\n")
    print(HTML_HEAD)
    print(HTML_FORM)

    if "predict_action" not in fs:
        print(HTML_FOOTER)
        return

    seq_mode      = fs.getvalue("seq_mode")      or "text"
    ref_mode      = fs.getvalue("ref_mode")      or "none"
    probabilities = fs.getvalue("probabilities") == "1"
    show_jmol     = fs.getvalue("show_jmol")     == "1"

    tmp_seq = tmp_model = tmp_ref = None
    seq_path = model_path = ref_path = None
    pdb_id_used = None
    error_msg = None

    # Modell
    if "model_file" in fs and fs["model_file"].filename:
        tmp_m = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
        tmp_m.write(fs["model_file"].file.read())
        tmp_m.close()
        model_path = tmp_m.name
        tmp_model  = model_path
    else:
        error_msg = "Bitte ein Modell hochladen."

    # Sequenz
    if not error_msg:
        if seq_mode == "text":
            seq_text = (fs.getvalue("seq_text") or "").strip()
            if not seq_text:
                error_msg = "Bitte eine Sequenz eingeben."
            else:
                tmp_s = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode='w')
                tmp_s.write(seq_text); tmp_s.close()
                seq_path = tmp_s.name; tmp_seq = seq_path

        elif seq_mode == "pdbid":
            pdb_id_used = (fs.getvalue("pdb_id") or "").strip().upper()
            if not pdb_id_used:
                error_msg = "Bitte eine PDB ID eingeben."
            else:
                print(f"<div class='status-msg info'>Lade PDB <strong>{pdb_id_used}</strong> von rcsb.org…</div>")
                seq_path, error_msg = fetch_fasta_from_pdb(pdb_id_used)
                tmp_seq = seq_path

        elif seq_mode == "upload":
            if "seq_file" in fs and fs["seq_file"].filename:
                allowed = {'.fasta', '.fa', '.txt'}
                suffix  = os.path.splitext(fs["seq_file"].filename)[1].lower()
                if suffix not in allowed:
                    error_msg = f"Ungültiges Format '{suffix}'. Erlaubt: .fasta, .fa"
                else:
                    tmp_s = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                    tmp_s.write(fs["seq_file"].file.read()); tmp_s.close()
                    seq_path = tmp_s.name; tmp_seq = seq_path
            else:
                error_msg = "Bitte eine FASTA-Datei hochladen."

    #Referenz
    if not error_msg and ref_mode != "none":
        if ref_mode == "upload":
            if "ref_file" in fs and fs["ref_file"].filename:
                tmp_r = tempfile.NamedTemporaryFile(suffix=".seclib", delete=False)
                tmp_r.write(fs["ref_file"].file.read()); tmp_r.close()
                ref_path = tmp_r.name; tmp_ref = ref_path
            else:
                error_msg = "Bitte eine Referenz-Datei hochladen."
        elif ref_mode == "text":
            ref_text = (fs.getvalue("ref_text") or "").strip()
            if not ref_text:
                error_msg = "Bitte Referenz-Inhalt einfügen."
            else:
                tmp_r = tempfile.NamedTemporaryFile(suffix=".seclib", delete=False, mode='w')
                tmp_r.write(ref_text); tmp_r.close()
                ref_path = tmp_r.name; tmp_ref = ref_path

    #Fehler
    if error_msg:
        print(f"<div class='status-msg error'><strong>Fehler:</strong> {html.escape(error_msg)}</div>")
        for t in [tmp_seq, tmp_model, tmp_ref]:
            if t and os.path.exists(t): os.unlink(t)
        print(HTML_FOOTER)
        return

    # Vorhersage laufen lassen
    res = run_gor_predict(model_path, seq_path, "html", probabilities)

    if not isinstance(res, subprocess.CompletedProcess):
        print(f"<div class='status-msg error'><strong>Systemfehler:</strong> {html.escape(str(res))}</div>")
        print(HTML_FOOTER)
        return

    if res.returncode != 0 or not res.stdout.strip():
        stderr = res.stderr
        if "keinen gültigen FASTA-Header" in stderr or "startet nicht mit" in stderr:
            msg = "Ungültige Sequenz: FASTA-Format erwartet (&gt;ID, dann Sequenz)."
        elif "loadModel" in stderr or "NumberFormatException" in stderr or "NoSuchElementException" in stderr:
            msg = "Ungültiges Modell."
        elif "IOException" in stderr:
            msg = "Fehler beim Lesen der Eingabedatei."
        elif not res.stdout.strip():
            msg = "Keine Ausgabe – Sequenz muss mindestens 17 Aminosäuren lang sein."
        else:
            msg = f"Fehler bei der Vorhersage.<br><small style='font-family:monospace'>{html.escape(stderr[:300])}</small>"
        print(f"<div class='status-msg error'><strong>Fehler:</strong> {msg}</div>")
        print(HTML_FOOTER)
        return

    # Vorhersage anzeigen
    print("<div class='status-msg success'><strong>Vorhersage erfolgreich!</strong></div>")
    print("<div class='result-title'>Vorhersage:</div>")
    print(f"<div class='result-html-box'>{res.stdout}</div>")

    # TXT Download
    res_txt = run_gor_predict(model_path, seq_path, "txt", probabilities)
    if isinstance(res_txt, subprocess.CompletedProcess) and res_txt.stdout.strip():
        b64 = base64.b64encode(res_txt.stdout.encode()).decode()
        print(f'<a href="data:text/plain;base64,{b64}" download="prediction.txt" class="download-btn">&#11015; prediction.txt herunterladen</a>')

    #Validierungsscores
    if ref_path and isinstance(res_txt, subprocess.CompletedProcess) and res_txt.stdout.strip():
        # Header normalisieren: "> 1TIM_1|Chains A, B|..." -> "> 1TIM_1"
        def normalize_prd(text):
            lines = []
            for line in text.splitlines():
                if line.startswith(">"):
                    lines.append("> " + line[1:].strip().split()[0])
                elif line.startswith("PS "):
                    lines.append("PS " + line[3:].replace("-", "C"))
                else:
                    lines.append(line)
            return "\n".join(lines) + "\n"

        tmp_prd = tempfile.NamedTemporaryFile(suffix=".prd", delete=False, mode='w')
        tmp_prd.write(normalize_prd(res_txt.stdout)); tmp_prd.close()
        tmp_sum = tempfile.NamedTemporaryFile(suffix=".sum", delete=False); tmp_sum.close()
        tmp_det = tempfile.NamedTemporaryFile(suffix=".evl", delete=False); tmp_det.close()

        val_res = run_validate(tmp_prd.name, ref_path, tmp_sum.name, tmp_det.name)

        if isinstance(val_res, subprocess.CompletedProcess) and val_res.returncode == 0:
            print("<div class='result-title'>Validierungsscores:</div>")

            if os.path.exists(tmp_sum.name) and os.path.getsize(tmp_sum.name) > 0:
                with open(tmp_sum.name) as f:
                    summary_txt = f.read()
                scores = parse_summary_scores(summary_txt)

                def score_class(v):
                    if v >= 65: return "score-good"
                    if v >= 50: return "score-mid"
                    return "score-low"

                metrics = [
                    ("Q3 (gesamt)", "q3"), ("SOV (gesamt)", "SOV"),
                    ("Q_H (Helix)", "qObs_H"), ("Q_E (Strang)", "qObs_E"), ("Q_C (Coil)", "qObs_C"),
                    ("SOV_H", "SOV_H"), ("SOV_E", "SOV_E"), ("SOV_C", "SOV_C"),
                ]
                print("<table class='score-table'>")
                print("<tr><th>Metrik</th><th>Wert</th></tr>")
                for name, key in metrics:
                    val = scores.get(key)
                    if val is not None:
                        print(f"<tr><td>{name}</td><td class='{score_class(val)}'>{val:.1f}%</td></tr>")
                print("</table>")

                b64s = base64.b64encode(summary_txt.encode()).decode()
                print(f'<a href="data:text/plain;base64,{b64s}" download="validation_summary.txt" class="download-btn">&#11015; Summary</a>')

            if os.path.exists(tmp_det.name) and os.path.getsize(tmp_det.name) > 0:
                with open(tmp_det.name) as f:
                    det_txt = f.read()
                b64d = base64.b64encode(det_txt.encode()).decode()
                print(f'<a href="data:text/plain;base64,{b64d}" download="validation_detailed.txt" class="download-btn">&#11015; Details</a>')
        else:
            print("<div class='status-msg info'>Validierung nicht möglich – Formate prüfen.</div>")

        for t in [tmp_prd.name, tmp_sum.name, tmp_det.name]:
            if os.path.exists(t): os.unlink(t)

    #JSmol
    if show_jmol and pdb_id_used:
        print("<div class='result-title'>3D-Struktur (JSmol):</div>")
        jmol_script = (
            f"load ={pdb_id_used}; background white; cartoons on; wireframe off; rockets off; "
            "select helix; color [230,57,70]; "
            "select sheet; color [244,162,97]; "
            "select coil;  color [173,181,189]; "
            "spin off; zoom 80;"
        )
        pdb_escaped = html.escape(pdb_id_used)
        print(f"""<div class="jmol-wrap">
      <h4>Sekundärstruktur &ndash; {pdb_escaped}</h4>
      <div id="jmol_div"></div>
      <div class="jmol-legend">
        <span><span class="leg-dot" style="background:#e63946"></span> Helix (H)</span>
        <span><span class="leg-dot" style="background:#f4a261"></span> Beta-Strang (E)</span>
        <span><span class="leg-dot" style="background:#adb5bd"></span> Coil / Loop (C)</span>
      </div>
    </div>
    <script>
      var jmolInfo = {{
        width: 500, height: 440,
        use: "HTML5",
        color: "#f8fcfb",
        j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
        script: "{jmol_script}"
      }};
      var myJmol = Jmol.getApplet("myJmol", jmolInfo);
      document.getElementById("jmol_div").innerHTML = Jmol.getAppletHtml(myJmol);
    </script>""")


    #afräumen
    for t in [tmp_seq, tmp_model, tmp_ref]:
        if t and os.path.exists(t): os.unlink(t)

    print(HTML_FOOTER)


if __name__ == "__main__":
    main()
