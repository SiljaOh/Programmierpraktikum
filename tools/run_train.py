#!/usr/bin/python3

import subprocess
import html
import base64
import os
import cgi
import tempfile
import urllib.request
from mycgi import Form

CGI_DIR = os.path.dirname(os.path.abspath(__file__))
JAR_PATH = os.path.join(CGI_DIR, "train.jar")
CB513_PATH = "/mnt/extstud/praktikum/bioprakt/Data/GOR/CB513/CB513DSSP.db"
GET_SECLIB = os.path.join(CGI_DIR, "get_seclib.py")


def run_gor_train(db_path, method, model_out):
    cmd = [
        "java", "-cp", JAR_PATH, "de.lmu.ifi.bio.gor.cli.Train",
        "--db", db_path,
        "--method", str(method),
        "--model", model_out
    ]
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=300
        )
        return result
    except Exception as e:
        return e


def is_model_empty(model_path):
    with open(model_path, "r") as f:
        content = f.read()
    numbers = [x for x in content.split() if x.lstrip('-').isdigit()]
    return bool(numbers) and all(x == '0' for x in numbers)


def fetch_pdb_as_seclib(pdb_id):
    pdb_id = pdb_id.strip()
    if not pdb_id or len(pdb_id) != 4 or not pdb_id.isalnum():
        return None, f"Ungültige PDB ID: '{pdb_id}'. Eine PDB ID besteht aus genau 4 alphanumerischen Zeichen (z.B. 1TIM, 2HHB)."

    tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp_seclib = tempfile.NamedTemporaryFile(suffix=".seclib", delete=False, mode='w')
    tmp_pdb.close()
    tmp_seclib.close()

    try:
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=30) as response:
                if response.status != 200:
                    return None, f"PDB ID '{pdb_id.upper()}' nicht gefunden. Bitte prüfen Sie die ID auf rcsb.org."
                with open(tmp_pdb.name, 'wb') as f:
                    f.write(response.read())
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, f"Ungültige PDB ID '{pdb_id.upper()}': Eintrag existiert nicht in der PDB."
            return None, f"HTTP Fehler {e.code} beim Laden von PDB '{pdb_id.upper()}'."
        except urllib.error.URLError as e:
            return None, f"Netzwerkfehler: {str(e.reason)}"

        if os.path.getsize(tmp_pdb.name) < 100:
            return None, f"PDB ID '{pdb_id.upper()}' lieferte keine gültige Datei."

        result = subprocess.run(
            ["python3", GET_SECLIB, "--PDB", tmp_pdb.name, "--output", tmp_seclib.name],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            timeout=60
        )
        if result.returncode != 0:
            return None, f"Fehler bei der Konvertierung: {result.stderr}"

        return tmp_seclib.name, None
    except Exception as e:
        return None, f"Unerwarteter Fehler: {str(e)}"
    finally:
        if os.path.exists(tmp_pdb.name):
            os.unlink(tmp_pdb.name)


def main():
    fs = cgi.FieldStorage()
    form = Form()

    print("Content-Type: text/html; charset=utf-8\n")

    print("""
<!DOCTYPE html>
<html lang="de">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GOR Model Trainer</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }

        body {
            font-family: 'Segoe UI', system-ui, sans-serif;
            background: #f0f4f3;
            color: #1a2e2a;
            min-height: 100vh;
        }

        .site-header {
            background: linear-gradient(135deg, #1a3a30 0%, #2d5a45 60%, #3a7a5a 100%);
            color: white;
            text-align: center;
            padding: 40px 20px 35px;
        }

        .site-header h1 {
            font-size: 1.9em;
            font-weight: 700;
            letter-spacing: 0.01em;
            margin-bottom: 6px;
        }

        .site-header .subtitle {
            font-size: 0.95em;
            opacity: 0.75;
            letter-spacing: 0.03em;
        }

        .page-content {
            max-width: 800px;
            margin: 0 auto;
            padding: 35px 20px 60px;
        }

        .card {
            background: white;
            border-radius: 16px;
            padding: 30px;
            margin-bottom: 20px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.06);
        }

        .card-title {
            font-size: 0.8em;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: #3a7a5a;
            margin-bottom: 18px;
        }

        .tab-buttons {
            display: flex;
            gap: 8px;
            margin-bottom: 20px;
        }

        .tab-btn {
            flex: 1;
            padding: 10px 14px;
            border: 2px solid #d4e8df;
            border-radius: 10px;
            background: white;
            cursor: pointer;
            font-weight: 600;
            font-size: 0.9em;
            color: #4a7a65;
            transition: all 0.2s;
        }

        .tab-btn:hover {
            border-color: #3a7a5a;
            background: #f0f8f4;
        }

        .tab-btn.active {
            border-color: #3a7a5a;
            background: #3a7a5a;
            color: white;
        }

        .tab-content { display: none; }
        .tab-content.active { display: block; }

        label {
            display: block;
            font-weight: 600;
            font-size: 0.9em;
            color: #2d4a40;
            margin-bottom: 8px;
        }

        select, input[type="text"] {
            width: 100%;
            padding: 11px 14px;
            border: 1.5px solid #c8e0d4;
            border-radius: 10px;
            font-size: 0.95em;
            color: #1a2e2a;
            background: #f8fcfa;
            outline: none;
            transition: border-color 0.2s;
        }

        select:focus, input[type="text"]:focus {
            border-color: #3a7a5a;
            background: white;
        }

        input[type="file"] {
            width: 100%;
            padding: 10px 14px;
            border: 1.5px dashed #c8e0d4;
            border-radius: 10px;
            font-size: 0.9em;
            background: #f8fcfa;
            cursor: pointer;
        }

        small {
            display: block;
            margin-top: 6px;
            font-size: 0.82em;
            color: #7a9e90;
        }

        .submit-btn {
            width: 100%;
            padding: 14px;
            background: linear-gradient(135deg, #2d5a45, #3a7a5a);
            color: white;
            border: none;
            border-radius: 12px;
            font-size: 1em;
            font-weight: 700;
            cursor: pointer;
            letter-spacing: 0.03em;
            transition: opacity 0.2s;
            margin-top: 5px;
        }

        .submit-btn:hover { opacity: 0.9; }

        .status-msg {
            padding: 14px 18px;
            border-radius: 12px;
            margin: 16px 0;
            font-size: 0.95em;
            border-left: 4px solid;
        }

        .error {
            background: #fdf0ef;
            color: #8b2a22;
            border-color: #c0392b;
        }

        .success {
            background: #edf7f2;
            color: #1e6644;
            border-color: #27ae60;
        }

        .info {
            background: #eef5fb;
            color: #1a4d7a;
            border-color: #2980b9;
        }

        .preview-box {
            background: #1a2e2a;
            color: #a8d5bc;
            padding: 20px;
            border-radius: 12px;
            font-family: 'Courier New', monospace;
            font-size: 0.85em;
            white-space: pre;
            max-height: 280px;
            overflow-y: auto;
            overflow-x: auto;
            margin-top: 16px;
        }

        .preview-title {
            font-weight: 700;
            font-size: 0.95em;
            color: #2d4a40;
            margin: 20px 0 6px;
        }

        .download-btn {
            display: inline-block;
            background: linear-gradient(135deg, #2d5a45, #3a7a5a);
            color: white;
            padding: 12px 24px;
            text-decoration: none;
            border-radius: 10px;
            font-weight: 700;
            font-size: 0.95em;
            margin-top: 16px;
            transition: opacity 0.2s;
        }

        .download-btn:hover { opacity: 0.9; }

        .back-link {
            text-align: center;
            margin-top: 30px;
        }

        .back-link a {
            color: #3a7a5a;
            text-decoration: none;
            font-weight: 600;
            font-size: 0.9em;
        }

        .back-link a:hover { text-decoration: underline; }
    </style>
    <script>
        function switchTab(tabId) {
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
            document.getElementById('tab-' + tabId).classList.add('active');
            document.getElementById('content-' + tabId).classList.add('active');
            document.getElementById('input_mode').value = tabId;
        }
    </script>
</head>
<body>

<div class="site-header">
    <h1>GOR Model Training</h1>
    <div class="subtitle">WS 2025/2026 &middot; Gruppe 17</div>
</div>

<div class="page-content">
    <form method="post" enctype="multipart/form-data">
        <input type="hidden" name="input_mode" id="input_mode" value="internal">

        <div class="card">
            <div class="card-title">1. Trainingsdaten wählen</div>

            <div class="tab-buttons">
                <button type="button" class="tab-btn active" id="tab-internal" onclick="switchTab('internal')">Interne DB</button>
                <button type="button" class="tab-btn" id="tab-pdbid" onclick="switchTab('pdbid')">PDB ID</button>
                <button type="button" class="tab-btn" id="tab-upload" onclick="switchTab('upload')">Datei hochladen</button>
            </div>

            <div class="tab-content active" id="content-internal">
                <label>Interne Datenbank auswählen</label>
                <select name="internal_db">
                    <option value="cb513">CB513 (513 Proteine, Benchmark)</option>
                </select>
            </div>

            <div class="tab-content" id="content-pdbid">
                <label>PDB ID eingeben (z.B. 1TIM)</label>
                <input type="text" name="pdb_id" placeholder="z.B. 1TIM, 2HHB...">
                <small>Sequenz wird automatisch von rcsb.org geladen. Nur Protein-Strukturen werden unterstützt.</small>
            </div>

            <div class="tab-content" id="content-upload">
                <label>Datei hochladen (.seclib, .db, .pdb, .dssp)</label>
                <input type="file" name="upload_file" accept=".seclib,.db,.pdb,.dssp">
                <small>Erlaubte Formate: seclib, db, pdb, dssp</small>
            </div>
        </div>

        <div class="card">
            <div class="card-title">2. GOR Algorithmus wählen</div>
            <label>Methode</label>
            <select name="method">
                <option value="gor1">GOR I</option>
                <option value="gor3">GOR III</option>
                <option value="gor4">GOR IV</option>
            </select>
        </div>

        <div class="card">
            <div class="card-title">3. Modell-Ausgabedatei</div>
            <label>Dateiname</label>
            <input type="text" name="model_name" value="my_gor_model.txt">
        </div>

        <button type="submit" name="train_action" class="submit-btn">Training starten</button>
    </form>
    """)

    if "train_action" in fs:
        input_mode = fs.getvalue("input_mode") or "internal"
        method = fs.getvalue("method") or "gor1"
        model_name = (fs.getvalue("model_name") or "my_gor_model.txt").strip()
        tmp_file = None
        final_db = None
        error_msg = None

        if input_mode == "internal":
            final_db = CB513_PATH

        elif input_mode == "pdbid":
            pdb_id = (fs.getvalue("pdb_id") or "").strip()
            if not pdb_id:
                error_msg = "Bitte eine PDB ID eingeben."
            else:
                print(f"<div class='status-msg info'>Lade PDB <strong>{pdb_id.upper()}</strong> von rcsb.org...</div>")
                final_db, error_msg = fetch_pdb_as_seclib(pdb_id)
                tmp_file = final_db

        elif input_mode == "upload":
            if "upload_file" in fs and fs["upload_file"].filename:
                uploaded = fs["upload_file"]
                allowed = {'.seclib', '.db', '.pdb', '.dssp'}
                suffix = os.path.splitext(uploaded.filename)[1].lower()
                if suffix not in allowed:
                    error_msg = f"Ungültiges Dateiformat '{suffix}'. Erlaubte Formate: .seclib, .db, .pdb, .dssp"
                else:
                    tmp = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
                    tmp.write(uploaded.file.read())
                    tmp.close()
                    final_db = tmp.name
                    tmp_file = final_db
            else:
                error_msg = "Bitte eine Datei auswählen."

        if error_msg:
            print(f"<div class='status-msg error'><strong>Fehler:</strong> {html.escape(error_msg)}</div>")
        elif final_db:
            res = run_gor_train(final_db, method, model_name)

            if isinstance(res, subprocess.CompletedProcess):
                if res.returncode == 0:
                    if os.path.exists(model_name) and is_model_empty(model_name):
                        print("<div class='status-msg error'><strong>Fehler:</strong> Die PDB ID enthält keine Protein-Sequenz (z.B. RNA oder DNA). GOR funktioniert nur mit Proteinen — bitte eine andere PDB ID verwenden.</div>")
                    elif os.path.exists(model_name):
                        print("<div class='status-msg success'><strong>Training erfolgreich!</strong> Das Modell wurde erstellt.</div>")

                        with open(model_name, "r") as f:
                            preview = "".join([f.readline() for _ in range(20)])
                        print("<div class='preview-title'>Modell-Vorschau (Header):</div>")
                        print(f"<div class='preview-box'>{html.escape(preview)}</div>")

                        with open(model_name, "rb") as f:
                            b64 = base64.b64encode(f.read()).decode()
                        print(f'<a href="data:text/plain;base64,{b64}" download="{model_name}" class="download-btn">&#11015; {model_name} herunterladen</a>')
                else:
                    print(f"<div class='status-msg error'><strong>Java-Fehler:</strong><br><pre>{html.escape(res.stderr)}</pre></div>")
            else:
                print(f"<div class='status-msg error'><strong>Systemfehler:</strong> {html.escape(str(res))}</div>")

            if tmp_file and os.path.exists(tmp_file):
                os.unlink(tmp_file)

    print("""
    <div class="back-link">
        <a href="../index.html">&larr; Zurück zur Hauptseite</a>
    </div>
</div>
</body>
</html>
    """)


if __name__ == "__main__":
    main()