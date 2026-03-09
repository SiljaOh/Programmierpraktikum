#!/usr/bin/python3

import subprocess
import html
import base64
import os
import cgi
import tempfile
import zipfile
import shutil

CGI_DIR = os.path.dirname(os.path.abspath(__file__))
JAR_PATH = os.path.join(CGI_DIR, "predict.jar")


def run_gor5_predict(model_path, maf_dir, fmt, probabilities):
    cmd = [
        "java", "-cp", JAR_PATH, "de.lmu.ifi.bio.gor.cli.Predict",
        "--model", model_path,
        "--maf", maf_dir,
        "--format", fmt,
    ]
    if probabilities:
        cmd.append("--probabilities")
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


def main():
    fs = cgi.FieldStorage()

    print("Content-Type: text/html; charset=utf-8\n")

    print("""
<!DOCTYPE html>
<html lang="de">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GOR5 – Multiple Alignment Prediction</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
    <style>
        :root {
            --bio-green: #0d6b5c;
            --bio-green-light: #d4ede8;
            --bio-green-pale: #e8f6f3;
            --bio-teal: #2a9d8f;
            --bio-dark: #0f2926;
            --bio-deep: #0a1f1c;
            --text: #1a2e2a;
            --text-muted: #4a6560;
            --card-bg: #ffffff;
            --border: #c5e0db;
        }

        * { box-sizing: border-box; margin: 0; padding: 0; }

        body {
            font-family: "DM Sans", system-ui, sans-serif;
            background: var(--bio-deep);
            color: var(--text);
            line-height: 1.6;
            min-height: 100vh;
        }

        body::before {
            content: "";
            position: fixed;
            inset: 0;
            background:
                linear-gradient(180deg, rgba(13,107,92,0.06) 0%, transparent 40%),
                radial-gradient(ellipse 80% 50% at 50% -20%, rgba(42,157,143,0.15), transparent),
                repeating-linear-gradient(0deg, transparent, transparent 32px, rgba(13,107,92,0.03) 32px, rgba(13,107,92,0.03) 33px),
                repeating-linear-gradient(90deg, transparent, transparent 32px, rgba(13,107,92,0.03) 32px, rgba(13,107,92,0.03) 33px);
            pointer-events: none;
            z-index: 0;
        }

        .header {
            position: relative;
            background: linear-gradient(160deg, var(--bio-deep) 0%, var(--bio-dark) 40%, var(--bio-green) 100%);
            color: #fff;
            padding: 40px 24px 48px;
            text-align: center;
            overflow: hidden;
        }

        .header::after {
            content: "";
            position: absolute;
            bottom: 0; left: 0; right: 0;
            height: 60px;
            background: linear-gradient(to top, var(--bio-deep), transparent);
            pointer-events: none;
        }

        .header-inner { position: relative; z-index: 1; }

        .header h1 {
            margin: 0 0 6px 0;
            font-size: 1.6rem;
            font-weight: 700;
            letter-spacing: 0.02em;
        }

        .header p {
            margin: 0;
            opacity: 0.85;
            font-size: 0.95rem;
        }

        .badge {
            display: inline-block;
            background: rgba(255,255,255,0.15);
            border: 1px solid rgba(255,255,255,0.3);
            border-radius: 20px;
            padding: 3px 14px;
            font-size: 0.78rem;
            font-weight: 600;
            letter-spacing: 0.06em;
            margin-top: 10px;
        }

        .hero-dots {
            display: flex;
            gap: 8px;
            justify-content: center;
            margin-top: 14px;
        }

        .hero-dots span {
            width: 8px; height: 8px;
            border-radius: 50%;
            background: rgba(255,255,255,0.45);
        }

        .main {
            position: relative;
            z-index: 1;
            max-width: 900px;
            margin: -24px auto 0;
            padding: 0 20px 56px;
        }

        .main-panel {
            background: linear-gradient(180deg, rgba(255,255,255,0.97) 0%, rgba(248,252,251,0.98) 100%);
            border-radius: 20px;
            padding: 36px 32px 40px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.25), 0 0 0 1px rgba(255,255,255,0.1) inset;
        }

        .card {
            background: var(--card-bg);
            border: 1px solid var(--border);
            border-radius: 14px;
            padding: 24px 24px 20px;
            margin-bottom: 20px;
        }

        .card-title {
            font-size: 0.78em;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.08em;
            color: var(--bio-green);
            margin-bottom: 16px;
        }

        .info-box {
            background: var(--bio-green-pale);
            border: 1px solid var(--bio-green-light);
            border-radius: 10px;
            padding: 11px 15px;
            font-size: 0.88em;
            color: var(--bio-green);
            margin-bottom: 16px;
        }

        label {
            display: block;
            font-weight: 600;
            font-size: 0.9em;
            color: var(--text);
            margin-bottom: 7px;
        }

        input[type="file"] {
            width: 100%;
            padding: 10px 13px;
            border: 1.5px dashed var(--border);
            border-radius: 10px;
            font-size: 0.88em;
            font-family: "DM Sans", sans-serif;
            background: var(--bio-green-pale);
            cursor: pointer;
        }

        small {
            display: block;
            margin-top: 6px;
            font-size: 0.82em;
            color: var(--text-muted);
        }

        .checkbox-row {
            display: flex;
            align-items: center;
            gap: 10px;
        }

        .checkbox-row input[type="checkbox"] {
            width: 17px;
            height: 17px;
            accent-color: var(--bio-green);
            cursor: pointer;
        }

        .checkbox-row label {
            margin: 0;
            font-weight: 500;
            cursor: pointer;
            font-size: 0.93em;
        }

        .submit-btn {
            width: 100%;
            padding: 14px;
            background: linear-gradient(135deg, var(--bio-dark), var(--bio-green));
            color: white;
            border: none;
            border-radius: 12px;
            font-size: 1em;
            font-weight: 700;
            font-family: "DM Sans", sans-serif;
            cursor: pointer;
            letter-spacing: 0.03em;
            transition: opacity 0.2s;
            margin-top: 5px;
        }

        .submit-btn:hover { opacity: 0.88; }

        .status-msg {
            padding: 13px 17px;
            border-radius: 12px;
            margin: 16px 0;
            font-size: 0.93em;
            border-left: 4px solid;
        }

        .error { background: #fdf0ef; color: #8b2a22; border-color: #c0392b; }
        .success { background: var(--bio-green-pale); color: var(--bio-dark); border-color: var(--bio-green); }
        .info { background: #eef5fb; color: #1a4d7a; border-color: #2980b9; }

        .result-html-box {
            background: #ffffff;
            color: #1a2e2a;
            padding: 24px 28px;
            border-radius: 14px;
            font-family: "JetBrains Mono", monospace;
            font-size: 0.95em;
            max-height: 600px;
            overflow-y: auto;
            overflow-x: auto;
            margin-top: 14px;
            border: 1.5px solid var(--border);
            line-height: 1.8;
            white-space: pre;
            box-shadow: 0 2px 12px rgba(13,107,92,0.08);
        }

        .result-title {
            font-weight: 700;
            font-size: 0.9em;
            color: var(--bio-green);
            margin: 20px 0 4px;
            text-transform: uppercase;
            letter-spacing: 0.06em;
        }

        .download-btn {
            display: inline-block;
            background: linear-gradient(135deg, var(--bio-dark), var(--bio-green));
            color: white;
            padding: 11px 22px;
            text-decoration: none;
            border-radius: 10px;
            font-weight: 700;
            font-size: 0.92em;
            margin-top: 14px;
            transition: opacity 0.2s;
        }

        .download-btn:hover { opacity: 0.88; }

        footer {
            position: relative;
            z-index: 1;
            background: linear-gradient(180deg, var(--bio-dark) 0%, var(--bio-deep) 100%);
            color: rgba(255,255,255,0.88);
            padding: 22px 24px;
            font-size: 0.88rem;
            text-align: center;
            border-top: 1px solid rgba(13,107,92,0.4);
            margin-top: 40px;
        }

        footer a { color: #9ae6d9; text-decoration: none; }
        footer a:hover { text-decoration: underline; }
    </style>
</head>
<body>

<header class="header">
    <div class="header-inner">
        <h1>GOR5 – Multiple Alignment Prediction</h1>
        <p>WS 2025/2026 <span>&middot;</span> Gruppe 17</p>
        <div class="badge">MSA-basierte Vorhersage</div>
        <div class="hero-dots"><span></span><span></span><span></span></div>
    </div>
</header>

<main class="main">
<div class="main-panel">
<form method="post" enctype="multipart/form-data">

    <div class="card">
        <div class="card-title">1. Alignment-Dateien hochladen</div>
        <div class="info-box">
            <strong>GOR5</strong> nutzt evolutionäre Information aus multiplen Sequenzalignments (MSA).
            Bitte ein ZIP-Archiv mit <code>.aln</code>-Dateien hochladen.
        </div>
        <label>ZIP-Archiv mit .aln-Dateien</label>
        <input type="file" name="maf_zip" accept=".zip">
        <small>Alle .aln-Dateien im ZIP werden als MSA-Eingabe verarbeitet.</small>
    </div>

    <div class="card">
        <div class="card-title">2. Modell hochladen</div>
        <label>Trainiertes GOR-Modell</label>
        <input type="file" name="model_file" accept=".txt,.model">
        <small>Modell-Datei aus dem Training-Tool (GOR I, III oder IV)</small>
    </div>

    <div class="card">
        <div class="card-title">3. Optionen</div>
        <div class="checkbox-row">
            <input type="checkbox" name="probabilities" id="probabilities" value="1">
            <label for="probabilities">Wahrscheinlichkeiten ausgeben (PH, PE, PC)</label>
        </div>
    </div>

    <button type="submit" name="predict_action" class="submit-btn">GOR5 Vorhersage starten</button>
</form>
    """)

    if "predict_action" in fs:
        probabilities = fs.getvalue("probabilities") == "1"
        tmp_model = None
        tmp_maf_dir = None
        model_path = None
        error_msg = None


        if "model_file" in fs and fs["model_file"].filename:
            uploaded_model = fs["model_file"]
            tmp_m = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
            tmp_m.write(uploaded_model.file.read())
            tmp_m.close()
            model_path = tmp_m.name
            tmp_model = model_path
        else:
            error_msg = "Bitte ein Modell hochladen."


        if not error_msg:
            if "maf_zip" in fs and fs["maf_zip"].filename:
                if not fs["maf_zip"].filename.endswith(".zip"):
                    error_msg = "Bitte eine ZIP-Datei hochladen."
                else:
                    tmp_zip = tempfile.NamedTemporaryFile(suffix=".zip", delete=False)
                    tmp_zip.write(fs["maf_zip"].file.read())
                    tmp_zip.close()
                    tmp_maf_dir = tempfile.mkdtemp()
                    try:
                        with zipfile.ZipFile(tmp_zip.name, 'r') as z:
                            aln_files = [f for f in z.namelist() if f.endswith(".aln")]
                            if not aln_files:
                                error_msg = "Keine .aln-Dateien im ZIP gefunden."
                            else:
                                z.extractall(tmp_maf_dir)
                    except zipfile.BadZipFile:
                        error_msg = "Ungültige ZIP-Datei."
                    finally:
                        os.unlink(tmp_zip.name)
            else:
                error_msg = "Bitte ein ZIP-Archiv mit .aln-Dateien hochladen."

        if error_msg:
            print(f"<div class='status-msg error'><strong>Fehler:</strong> {html.escape(error_msg)}</div>")
        elif model_path and tmp_maf_dir:
            # Ordner mit .aln-Dateien finden
            maf_dir = tmp_maf_dir
            for root, dirs, files in os.walk(tmp_maf_dir):
                if any(f.endswith(".aln") for f in files):
                    maf_dir = root
                    break

            res = run_gor5_predict(model_path, maf_dir, "html", probabilities)

            if isinstance(res, subprocess.CompletedProcess):
                if res.returncode == 0 and res.stdout.strip():
                    print("<div class='status-msg success'><strong>Vorhersage erfolgreich!</strong></div>")
                    print("<div class='result-title'>Ergebnis:</div>")
                    print(f"<div class='result-html-box'>{res.stdout}</div>")

                    # TXT für Download
                    res_txt = run_gor5_predict(model_path, maf_dir, "txt", probabilities)
                    if isinstance(res_txt, subprocess.CompletedProcess) and res_txt.stdout.strip():
                        b64 = base64.b64encode(res_txt.stdout.encode()).decode()
                        print(f'<a href="data:text/plain;base64,{b64}" download="gor5_prediction.txt" class="download-btn">&#11015; gor5_prediction.txt herunterladen</a>')
                elif res.returncode != 0:
                    stderr = res.stderr
                    if "loadModel" in stderr or "NumberFormatException" in stderr or "NoSuchElementException" in stderr:
                        error_msg = "Ungültiges Modell: Die hochgeladene Datei ist kein gültiges GOR-Modell. Bitte ein trainiertes Modell aus dem Training-Tool verwenden."
                    elif "IOException" in stderr:
                        error_msg = "Fehler beim Lesen der Eingabedatei. Bitte Dateiformat prüfen."
                    else:
                        error_msg = f"Unbekannter Fehler bei der Vorhersage. Bitte Eingaben prüfen.<br><small style='font-family:monospace'>{html.escape(stderr[:300])}</small>"
                    print(f"<div class='status-msg error'><strong>Fehler:</strong> {error_msg}</div>")
                else:
                    print("<div class='status-msg error'><strong>Fehler:</strong> Keine Ausgabe erzeugt. Bitte .aln-Dateien und Modell prüfen.</div>")
            else:
                print(f"<div class='status-msg error'><strong>Systemfehler:</strong> {html.escape(str(res))}</div>")

        # Aufräumen
        if tmp_model and os.path.exists(tmp_model):
            os.unlink(tmp_model)
        if tmp_maf_dir and os.path.exists(tmp_maf_dir):
            shutil.rmtree(tmp_maf_dir, ignore_errors=True)

    print("""
</div>
</main>

<footer>
    ProgPrakt WS 2025/2026 &middot; Gruppe 17 &nbsp;|&nbsp;
    <a href="../index.html">&larr; Zurück zur Hauptseite</a>
</footer>

</body>
</html>
    """)


if __name__ == "__main__":
    main()