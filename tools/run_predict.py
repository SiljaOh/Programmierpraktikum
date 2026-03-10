#!/usr/bin/python3

import subprocess
import html
import base64
import os
import cgi
import tempfile
import urllib.request
import re
import json

CGI_DIR           = os.path.dirname(os.path.abspath(__file__))
JAR_PATH          = os.path.join(CGI_DIR, "predict.jar")
VALIDATE_JAR_PATH = os.path.join(CGI_DIR, "validateGor.jar")


def get_cartoon_style(seq_str):
    mapping = {'H': '=', 'E': '>', 'C': '.', '-': ' '}
    return "".join([mapping.get(c.upper(), '.') for c in seq_str])


def get_cartoon_aligned(ps_line):
    if len(ps_line) >= 3:
        return "   " + get_cartoon_style(ps_line[3:])
    return "   " + get_cartoon_style(ps_line)


def run_gor_predict(model_path, seq_path, fmt, probabilities, show_scores):
    cmd = [
        "java", "-cp", JAR_PATH, "de.lmu.ifi.bio.gor.cli.Predict",
        "--model", model_path,
        "--seq",   seq_path,
        "--format", fmt,
    ]
    if probabilities:
        cmd.append("--probabilities")

    if show_scores:
        cmd.append("--show-scores")

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
    pdb_id = pdb_id.strip().upper()
    if not pdb_id or len(pdb_id) != 4 or not pdb_id.isalnum():
        return None, f"Ungültige PDB ID: '{pdb_id}'."
    tmp = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode='w')
    tmp.close()
    try:
        url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            content = resp.read().decode('utf-8')
            with open(tmp.name, 'w') as f:
                f.write(content)
        if os.path.getsize(tmp.name) < 10:
            return None, f"PDB ID '{pdb_id}' lieferte keine FASTA-Daten."
        return tmp.name, None
    except urllib.error.HTTPError as e:
        return None, f"HTTP {e.code} für PDB '{pdb_id}'."
    except Exception as e:
        return None, str(e)


def fetch_seclib_from_pdb(pdb_id):
    """
    Holt die mmCIF-Datei von RCSB und extrahiert Sekundärstruktur (H/E/C)
    per Kette+Residue aus _struct_conf (Helices) und _struct_sheet_range (Sheets).
    Baut daraus eine .seclib-Datei im Format:
        > PDBID_CHAIN
        AS <Aminosäuresequenz>
        PS <H/E/C-String>
    """
    pdb_id = pdb_id.strip().upper()
    if not pdb_id or len(pdb_id) != 4 or not pdb_id.isalnum():
        return None, None, f"Ungültige PDB ID: '{pdb_id}'."
    try:
        # mmCIF holen
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            cif_text = resp.read().decode('utf-8', errors='replace')
    except urllib.error.HTTPError as e:
        return None, None, f"HTTP {e.code} für PDB '{pdb_id}'."
    except Exception as e:
        return None, None, str(e)

    # FASTA auch holen für Sequenz und Kettenlängen
    fasta_path, err = fetch_fasta_from_pdb(pdb_id)
    if err:
        return None, err

    # FASTA einlesen: chain_id -> Sequenz
    chain_seqs = {}
    with open(fasta_path) as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and seq_parts:
                    # Header-Format: >PDBID_CHAIN|...
                    chain_id = header.split("|")[0].split("_")[-1] if "_" in header else "A"
                    chain_seqs[chain_id] = "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            elif line:
                seq_parts.append(line)
        if header and seq_parts:
            chain_id = header.split("|")[0].split("_")[-1] if "_" in header else "A"
            chain_seqs[chain_id] = "".join(seq_parts)
    os.unlink(fasta_path)

    if not chain_seqs:
        return None, None, f"Keine Sequenzdaten für '{pdb_id}' gefunden."

    # Helix-Ranges aus mmCIF parsen: chain -> [(start, end), ...]
    helix_ranges = {}
    sheet_ranges = {}

    # _struct_conf block (Helices: HELX_P)
    in_struct_conf = False
    conf_headers = []
    conf_rows = []
    lines = cif_text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "loop_":
            # Check ob nächster Block _struct_conf ist
            block_headers = []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                block_headers.append(lines[i].strip())
                i += 1
            if any("_struct_conf." in h for h in block_headers):
                # Daten lesen bis nächster loop_ oder #
                block_data = []
                while i < len(lines):
                    l = lines[i].strip()
                    if l.startswith("loop_") or l.startswith("#"):
                        break
                    if l and not l.startswith("_"):
                        block_data.append(l)
                    i += 1
                # Spalten-Indices ermitteln
                def col(name):
                    for idx, h in enumerate(block_headers):
                        if h == name:
                            return idx
                    return -1
                ci_type   = col("_struct_conf.conf_type_id")
                ci_chain  = col("_struct_conf.beg_label_asym_id")
                ci_start  = col("_struct_conf.beg_label_seq_id")
                ci_end    = col("_struct_conf.end_label_seq_id")
                if ci_type >= 0 and ci_chain >= 0 and ci_start >= 0 and ci_end >= 0:
                    for row in block_data:
                        parts = row.split()
                        if len(parts) > max(ci_type, ci_chain, ci_start, ci_end):
                            conf_type = parts[ci_type]
                            chain     = parts[ci_chain]
                            try:
                                s = int(parts[ci_start])
                                e = int(parts[ci_end])
                            except:
                                continue
                            if conf_type.startswith("HELX"):
                                helix_ranges.setdefault(chain, []).append((s, e))
            elif any("_struct_sheet_range." in h for h in block_headers):
                block_data = []
                while i < len(lines):
                    l = lines[i].strip()
                    if l.startswith("loop_") or l.startswith("#"):
                        break
                    if l and not l.startswith("_"):
                        block_data.append(l)
                    i += 1
                def col2(name):
                    for idx, h in enumerate(block_headers):
                        if h == name:
                            return idx
                    return -1
                ci_chain = col2("_struct_sheet_range.beg_label_asym_id")
                ci_start = col2("_struct_sheet_range.beg_label_seq_id")
                ci_end   = col2("_struct_sheet_range.end_label_seq_id")
                if ci_chain >= 0 and ci_start >= 0 and ci_end >= 0:
                    for row in block_data:
                        parts = row.split()
                        if len(parts) > max(ci_chain, ci_start, ci_end):
                            chain = parts[ci_chain]
                            try:
                                s = int(parts[ci_start])
                                e = int(parts[ci_end])
                            except:
                                continue
                            sheet_ranges.setdefault(chain, []).append((s, e))
            continue
        i += 1

    # Residue-Nummerierung: label_seq_id -> position in Sequenz (1-basiert)
    # Wir brauchen die Mapping-Tabelle chain -> [label_seq_ids in Reihenfolge]
    # Aus _atom_site: label_asym_id, label_seq_id, label_comp_id (nur CA-Atome)
    chain_resnum_order = {}
    in_atom = False
    atom_headers = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == "loop_":
            block_headers = []
            i += 1
            while i < len(lines) and lines[i].strip().startswith("_"):
                block_headers.append(lines[i].strip())
                i += 1
            if any("_atom_site." in h for h in block_headers):
                def acol(name):
                    for idx, h in enumerate(block_headers):
                        if h == name:
                            return idx
                    return -1
                ci_group  = acol("_atom_site.group_PDB")
                ci_name   = acol("_atom_site.label_atom_id")
                ci_chain  = acol("_atom_site.label_asym_id")
                ci_seqid  = acol("_atom_site.label_seq_id")
                ci_altloc = acol("_atom_site.label_alt_id")
                while i < len(lines):
                    l = lines[i].strip()
                    if l.startswith("loop_") or l.startswith("#"):
                        break
                    if l and not l.startswith("_"):
                        parts = l.split()
                        needed = max(ci_group, ci_name, ci_chain, ci_seqid) if all(x >= 0 for x in [ci_group, ci_name, ci_chain, ci_seqid]) else -1
                        if needed >= 0 and len(parts) > needed:
                            if parts[ci_group] == "ATOM" and parts[ci_name] == "CA":
                                alt = parts[ci_altloc] if ci_altloc >= 0 else "."
                                if alt in (".", "A"):
                                    chain  = parts[ci_chain]
                                    seqid  = parts[ci_seqid]
                                    try:
                                        seqnum = int(seqid)
                                    except:
                                        i += 1
                                        continue
                                    if chain not in chain_resnum_order:
                                        chain_resnum_order[chain] = []
                                    if seqnum not in chain_resnum_order[chain]:
                                        chain_resnum_order[chain].append(seqnum)
                    i += 1
                break
            continue
        i += 1

    # seclib aufbauen
    seclib_lines = []
    for chain, seq in sorted(chain_seqs.items()):
        res_nums = chain_resnum_order.get(chain, [])
        ss = ['C'] * len(seq)
        helices = helix_ranges.get(chain, [])
        sheets  = sheet_ranges.get(chain, [])
        for pos_idx, resnum in enumerate(res_nums):
            if pos_idx >= len(ss):
                break
            for (s, e) in helices:
                if s <= resnum <= e:
                    ss[pos_idx] = 'H'
                    break
            else:
                for (s, e) in sheets:
                    if s <= resnum <= e:
                        ss[pos_idx] = 'E'
                        break
        ss_str = "".join(ss)
        entry_id = f"{pdb_id}_{chain}"
        seclib_lines.append(f"> {entry_id}")
        seclib_lines.append(f"AS {seq}")
        seclib_lines.append(f"PS {ss_str}")
        seclib_lines.append("")

    if not seclib_lines:
        return None, None, f"Keine Sekundärstruktur-Daten für '{pdb_id}' extrahierbar."

    tmp = tempfile.NamedTemporaryFile(suffix=".seclib", delete=False, mode='w')
    tmp.write("\n".join(seclib_lines))
    tmp.close()
    return tmp.name, chain_resnum_order, None


def parse_summary_scores(summary_text):
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


def normalize_seq_input(text):
    lines = text.splitlines()
    blocks = []
    current_header = None
    current_seq = []
    auto_index = 1

    def flush():
        nonlocal current_header, current_seq, auto_index
        if current_seq:
            if current_header is None:
                current_header = f">Seq{auto_index}"
                auto_index += 1
            blocks.append((current_header, "".join(current_seq)))
            current_header = None
            current_seq = []

    for line in lines:
        line = line.strip()
        if not line:
            flush()
        elif line.startswith(">"):
            flush()
            current_header = line
        else:
            seq_clean = "".join(c for c in line if c.isalpha())
            if seq_clean:
                if current_header is None and current_seq:
                    flush()
                current_seq.append(seq_clean)
                if current_header is None:
                    flush()

    flush()
    return "\n\n".join(f"{h}\n{s}" for h, s in blocks) + "\n"


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
        .result-panel {
            background: linear-gradient(180deg, rgba(255,255,255,0.97) 0%, rgba(248,252,251,0.98) 100%);
            border-radius:20px; padding:36px 32px 40px; margin-top:24px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.25), 0 0 0 1px rgba(255,255,255,0.1) inset;
        }
        .result-panel-title {
            font-size:1.1em; font-weight:700; color:var(--bio-dark);
            margin-bottom:16px; padding-bottom:12px;
            border-bottom: 2px solid var(--bio-green-light);
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
        .result-title { font-weight:700; font-size:0.9em; color:var(--bio-green);
                        margin:20px 0 4px; text-transform:uppercase; letter-spacing:0.06em; }
        .download-btn {
            display:inline-block; background:linear-gradient(135deg, var(--bio-dark), var(--bio-green));
            color:white; padding:11px 22px; text-decoration:none; border-radius:10px;
            font-weight:700; font-size:0.92em; margin-top:14px; margin-right:8px; transition:opacity 0.2s;
        }
        .download-btn:hover { opacity:0.88; }
        .seq-block {
            background:#fff; border:1.5px solid var(--border); border-radius:14px;
            padding:20px 24px; margin-bottom:16px; font-family:"JetBrains Mono",monospace;
            font-size:0.88em; line-height:1.9; overflow-x:auto;
            box-shadow:0 2px 12px rgba(13,107,92,0.06);
        }
        .seq-block .seq-header { font-weight:700; color:var(--bio-green); font-size:0.95em; display:block; margin-bottom:6px; }
        .seq-block .seq-line   { display:block; white-space:pre; }
        .seq-block .cartoon-line { display:block; white-space:pre; letter-spacing:0; color:#2a9d8f; font-size:1em; margin-top:0; }
        .score-table { width:100%; border-collapse:collapse; margin-top:14px; font-family:"JetBrains Mono",monospace; font-size:0.9em; }
        .score-table th { background:var(--bio-green); color:white; padding:10px 16px; text-align:left; font-size:0.82em; letter-spacing:0.05em; }
        .score-table td { padding:9px 16px; border-bottom:1px solid var(--bio-green-light); }
        .score-table tr:last-child td { border-bottom:none; }
        .score-table tr:hover td { background:var(--bio-green-pale); }
        .score-good { color:#0d6b5c; font-weight:700; }
        .score-mid  { color:#e07b00; font-weight:600; }
        .score-low  { color:#c0392b; font-weight:600; }
        .pagination { display:flex; gap:10px; margin:18px 0 10px; align-items:center; }
        .page-btn {
            padding:9px 20px; border:2px solid var(--border); border-radius:9px;
            background:white; cursor:pointer; font-weight:600; font-size:0.88em;
            color:var(--bio-green); font-family:"DM Sans",sans-serif; transition:all 0.2s;
        }
        .page-btn:hover:not(:disabled) { border-color:var(--bio-teal); background:var(--bio-green-pale); }
        .page-btn:disabled { opacity:0.35; cursor:default; }
        .page-info { font-size:0.88em; color:var(--text-muted); flex:1; text-align:center; font-weight:500; }
        .jmol-wrap { background:#fff; border:1.5px solid var(--border); border-radius:14px;
                     padding:20px; margin-top:0; box-shadow:0 2px 12px rgba(13,107,92,0.08); }
        .jmol-wrap h4 { font-size:0.82em; font-weight:700; text-transform:uppercase;
                        letter-spacing:0.07em; color:var(--bio-green); margin-bottom:12px; }
        .jmol-legend { display:flex; gap:14px; flex-wrap:wrap; margin-top:12px; font-size:0.82em; color:var(--text-muted); }
        .jmol-legend span { display:flex; align-items:center; gap:6px; }
        .leg-dot { width:11px; height:11px; border-radius:50%; display:inline-block; }
        .jmol-two-col { display:grid; grid-template-columns:1fr 1fr; gap:16px; margin-top:8px; }
        @media(max-width:640px) { .jmol-two-col { grid-template-columns:1fr; } }
        footer { position:relative; z-index:1;
            background:linear-gradient(180deg, var(--bio-dark) 0%, var(--bio-deep) 100%);
            color:rgba(255,255,255,0.88); padding:22px 24px; font-size:0.88rem;
            text-align:center; border-top:1px solid rgba(13,107,92,0.4); margin-top:40px;
        }
        footer a { color:#9ae6d9; text-decoration:none; }
        footer a:hover { text-decoration:underline; }
        .info-icon {
            position:relative; display:inline-flex; align-items:center; justify-content:center;
            width:20px; height:20px; border-radius:50%; background:var(--bio-green);
            color:white; font-size:0.85em; font-weight:700; cursor:default;
            font-style:normal; flex-shrink:0; user-select:none;
        }
        .info-icon .tooltip {
            display:none; position:absolute; left:28px; top:50%; transform:translateY(-50%);
            background:#1a2e2a; color:#e8f6f3; padding:18px 20px; border-radius:12px;
            font-size:0.93em; font-weight:400; white-space:normal; z-index:100;
            box-shadow:0 6px 28px rgba(0,0,0,0.45); line-height:1.9;
            pointer-events:none; width:280px;
        }
        .info-icon .tooltip::before {
            content:""; position:absolute; right:100%; top:50%; transform:translateY(-50%);
            border:7px solid transparent; border-right-color:#1a2e2a;
        }
        .info-icon:hover .tooltip,
        .info-icon:focus .tooltip { display:block; }
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
        var _currentPage = 1;
        function showPage(p) {
            var perPage = 5;
            var total   = document.querySelectorAll('.seq-block').length;
            var pages   = Math.ceil(total / perPage);
            if (p < 1 || p > pages) return;
            _currentPage = p;
            document.querySelectorAll('.seq-block').forEach(function(b, i) {
                b.style.display = (i >= (p-1)*perPage && i < p*perPage) ? 'block' : 'none';
            });
            var info = document.getElementById('page-info');
            if (info) info.textContent = 'Seite ' + p + ' von ' + pages + ' (' + total + ' Sequenzen)';
            document.getElementById('btn-prev').disabled = (p === 1);
            document.getElementById('btn-next').disabled = (p === pages);
        }
        function prevPage() { showPage(_currentPage - 1); }
        function nextPage() { showPage(_currentPage + 1); }
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

    <div class="card">
        <div class="card-title">1. Sequenz eingeben</div>
        <div class="tab-buttons">
            <button type="button" class="seq-tab-btn tab-btn active" id="tab-text"
                    onclick="switchTab('text')">Text / FASTA</button>
            <button type="button" class="seq-tab-btn tab-btn" id="tab-pdbid"
                    onclick="switchTab('pdbid')">PDB ID</button>
            <button type="button" class="seq-tab-btn tab-btn" id="tab-upload"
                    onclick="switchTab('upload')">FASTA hochladen</button>
        </div>
        <div class="seq-tab-content tab-content active" id="content-text">
            <label style="display:flex;align-items:center;gap:7px;">
                Aminosäuresequenz
                <span class="info-icon" tabindex="0">&#9432;
                    <span class="tooltip">
                        <strong>Format 1 – mit Header:</strong><br>
                        &gt;ProteinName<br>ACDEFGHIKLM...<br><br>
                        <strong>Format 2 – ohne Header:</strong><br>
                        ACDEFGHIKLM...<br><br>
                        <strong>Mehrere Sequenzen:</strong><br>
                        Mit Headern oder durch<br>Zeilenumbrüche trennen.
                    </span>
                </span>
            </label>
            <textarea name="seq_text" rows="5" placeholder=">MeinProtein&#10;ACDEFGHIKLMNPQRSTVWY...&#10;&#10;oder einfach Sequenz ohne Header"></textarea>
        </div>
        <div class="seq-tab-content tab-content" id="content-pdbid">
            <label>PDB ID (z.B. 1TIM)</label>
            <input type="text" name="pdb_id" placeholder="z.B. 1TIM, 2HHB...">
            <small>Sequenz wird automatisch von rcsb.org geladen.</small>
        </div>
        <div class="seq-tab-content tab-content" id="content-upload">
            <label>FASTA-Datei hochladen</label>
            <input type="file" name="seq_file" accept=".fasta,.fa,.txt">
            <small>Erlaubte Formate: .fasta, .fa, .txt</small>
        </div>
    </div>

    <div class="card">
        <div class="card-title">2. Modell hochladen</div>
        <label>Trainiertes GOR-Modell</label>
        <input type="file" name="model_file" accept=".txt,.model">
        <small>Modell-Datei aus dem Training-Tool (GOR I, III oder IV)</small>
    </div>

    <div class="card">
        <div class="card-title">3. Referenzstruktur (optional) &ndash; für Q3 &amp; SOV Scores</div>
        <div class="tab-buttons">
            <button type="button" class="ref-tab-btn tab-btn active" id="ref-tab-none"
                    onclick="switchRefTab('none')">Keine</button>
            <button type="button" class="ref-tab-btn tab-btn" id="ref-tab-upload"
                    onclick="switchRefTab('upload')">Datei hochladen</button>
            <button type="button" class="ref-tab-btn tab-btn" id="ref-tab-text"
                    onclick="switchRefTab('text')">Text einfügen</button>
            <button type="button" class="ref-tab-btn tab-btn" id="ref-tab-pdbid"
                    onclick="switchRefTab('pdbid')">PDB ID</button>
        </div>
        <div class="ref-tab-content tab-content active" id="ref-content-none">
            <small style="margin-top:0">Kein Vergleich – nur Vorhersage anzeigen.</small>
        </div>
        <div class="ref-tab-content tab-content" id="ref-content-upload">
            <label>Referenz-Datei</label>
            <input type="file" name="ref_file" accept=".seclib,.txt,.db">
            <small>Erlaubte Formate: .seclib, .txt, .db</small>
        </div>
        <div class="ref-tab-content tab-content" id="ref-content-text">
            <label>Referenz einfügen</label>
            <textarea name="ref_text" rows="5" placeholder="&gt; ProteinID&#10;AS ACDEFG...&#10;PS HHHCCC..."></textarea>
        </div>
        <div class="ref-tab-content tab-content" id="ref-content-pdbid">
            <label>PDB ID als Referenz</label>
            <input type="text" name="ref_pdb_id" placeholder="z.B. 1TIM, 2HHB...">
            <small>Bekannte Sekundärstruktur wird automatisch via RCSB mmCIF geholt (Helix/Sheet/Coil).</small>
        </div>
    </div>

    <div class="card">
        <div class="card-title">4. Optionen</div>
        <div class="checkbox-row" style="margin-bottom:12px">
    <input type="checkbox" name="probabilities" id="probabilities" value="1">
    <label for="probabilities">
        Probability für die Berechnung benutzen
    </label>
        </div>

        <div class="checkbox-row" style="margin-bottom:12px">
            <input type="checkbox" name="show_scores" id="show_scores" value="1">
            <label for="show_scores">
                Scores anzeigen (-> Kein Post-Processing)
        </label>
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
    show_scores = fs.getvalue("show_scores") == "1"
    show_jmol     = fs.getvalue("show_jmol")     == "1"

    tmp_seq = tmp_model = tmp_ref = None
    seq_path = model_path = ref_path = None
    pdb_id_used = None
    error_msg = None

    if "model_file" in fs and fs["model_file"].filename:
        tmp_m = tempfile.NamedTemporaryFile(suffix=".txt", delete=False)
        tmp_m.write(fs["model_file"].file.read())
        tmp_m.close()
        model_path = tmp_m.name
        tmp_model  = model_path
    else:
        error_msg = "Bitte ein Modell hochladen."

    if not error_msg:
        if seq_mode == "text":
            seq_text = (fs.getvalue("seq_text") or "").strip()
            if not seq_text:
                error_msg = "Bitte eine Sequenz eingeben."
            else:
                seq_text = normalize_seq_input(seq_text)
                tmp_s = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False, mode='w')
                tmp_s.write(seq_text); tmp_s.close()
                seq_path = tmp_s.name; tmp_seq = seq_path

        elif seq_mode == "pdbid":
            pdb_id_used = (fs.getvalue("pdb_id") or "").strip().upper()
            if not pdb_id_used:
                error_msg = "Bitte eine PDB ID eingeben."
            else:
                print(f"<div class='status-msg info'>Lade Sequenz für PDB <strong>{html.escape(pdb_id_used)}</strong> von rcsb.org…</div>")
                seq_path, error_msg = fetch_fasta_from_pdb(pdb_id_used)
                tmp_seq = seq_path

        elif seq_mode == "upload":
            if "seq_file" in fs and fs["seq_file"].filename:
                tmp_s = tempfile.NamedTemporaryFile(suffix=".fasta", delete=False)
                tmp_s.write(fs["seq_file"].file.read()); tmp_s.close()
                seq_path = tmp_s.name; tmp_seq = seq_path
            else:
                error_msg = "Bitte eine FASTA-Datei hochladen."

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
        elif ref_mode == "pdbid":
            ref_pdb_id = (fs.getvalue("ref_pdb_id") or "").strip().upper()
            if not ref_pdb_id:
                error_msg = "Bitte eine PDB ID als Referenz eingeben."
            else:
                print(f"<div class='status-msg info'>Lade Referenzstruktur für PDB <strong>{html.escape(ref_pdb_id)}</strong> von rcsb.org…</div>")
                ref_path, _, error_msg = fetch_seclib_from_pdb(ref_pdb_id)
                tmp_ref = ref_path

    print("<div class='result-panel'>")
    print("<div class='result-panel-title'>Ergebnis</div>")

    if error_msg:
        print(f"<div class='status-msg error'><strong>Fehler:</strong> {html.escape(error_msg)}</div>")
        for t in [tmp_seq, tmp_model, tmp_ref]:
            if t and os.path.exists(t): os.unlink(t)
        print("</div>")
        print(HTML_FOOTER)
        return

    res_txt = run_gor_predict(model_path, seq_path, "html", probabilities, show_scores)

    if not isinstance(res_txt, subprocess.CompletedProcess):
        print(f"<div class='status-msg error'><strong>Systemfehler:</strong> {html.escape(str(res_txt))}</div>")
        print("</div>")
        print(HTML_FOOTER)
        return

    if res_txt.returncode != 0 or not res_txt.stdout.strip():
        stderr = res_txt.stderr
        if "keinen gültigen FASTA-Header" in stderr or "startet nicht mit" in stderr:
            msg = "Ungültige Sequenz: FASTA-Format erwartet (&gt;ID, dann Sequenz)."
        elif "loadModel" in stderr or "NumberFormatException" in stderr or "NoSuchElementException" in stderr:
            msg = "Ungültiges Modell."
        elif "IOException" in stderr:
            msg = "Fehler beim Lesen der Eingabedatei."
        elif not res_txt.stdout.strip():
            msg = "Keine Ausgabe – Sequenz muss mindestens 17 Aminosäuren lang sein."
        else:
            msg = f"Fehler bei der Vorhersage.<br><small style='font-family:monospace'>{html.escape(stderr[:300])}</small>"
        print(f"<div class='status-msg error'><strong>Fehler:</strong> {msg}</div>")
        print("</div>")
        print(HTML_FOOTER)
        return

    print("<div class='result-title'>Vorhersage:</div>")

    blocks = []

    if "<div class='gor-results'" in res_txt.stdout:

        # HTML direkt von Java anzeigen
        print(res_txt.stdout)

    else:

        current = []

        for line in res_txt.stdout.splitlines():

            if line.startswith(">"):

                if current:
                    blocks.append(current)

                current = [line]

            else:

                if line.strip():
                    current.append(line)

        if current:
            blocks.append(current)

    total = len(blocks)
    per_page = 5
    num_pages = max(1, (total + per_page - 1) // per_page)

    if blocks:

        for block in blocks:
            print("<div class='seq-block'>")
            for line in block:
                print(html.escape(line) + "<br>")
            print("</div>")

    if num_pages > 1:
        print("<div class='pagination'>")

        print(
            "<button class='page-btn' id='btn-prev' "
            "onclick='prevPage()' disabled>"
            "&#9664; Zurück</button>"
        )

        print(
            f"<span class='page-info' id='page-info'>"
            f"Seite 1 von {num_pages} ({total} Sequenzen)"
            "</span>"
        )

        next_disabled = " disabled" if num_pages == 1 else ""

        print(
            f"<button class='page-btn' id='btn-next' "
            f"onclick='nextPage()'{next_disabled}>"
            "Weiter &#9654;</button>"
        )

        print("</div>")

    for i, block in enumerate(blocks):
        display = "block" if i < per_page else "none"
        header  = block[0]
        print(f"<div class='seq-block' style='display:{display}'>")
        print(f"<span class='seq-header'>{html.escape(header)}</span>")
        for line in block[1:]:
            if line.startswith("PS "):
                print(f"<span class='seq-line'>{html.escape(line)}</span>")
                cartoon_line = get_cartoon_aligned(line)
                print(f"<span class='cartoon-line'>{html.escape(cartoon_line)}</span>")
            else:
                print(f"<span class='seq-line'>{html.escape(line)}</span>")
        print("</div>")

    b64 = base64.b64encode(res_txt.stdout.encode()).decode()
    print(f'<a href="data:text/plain;base64,{b64}" download="prediction.txt" class="download-btn">&#11015; prediction.txt herunterladen</a>')

    if ref_path:
        def normalize_prd(text):
            out = []
            for line in text.splitlines():
                if line.startswith(">"):
                    out.append("> " + line[1:].strip().split()[0])
                elif line.startswith("PS "):
                    out.append("PS " + line[3:].replace("-", "C"))
                else:
                    out.append(line)
            return "\n".join(out) + "\n"

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
            print("<div class='status-msg info'>Validierung nicht möglich – Formate oder Sequenz-IDs prüfen.</div>")

        for t in [tmp_prd.name, tmp_sum.name, tmp_det.name]:
            if os.path.exists(t): os.unlink(t)

    if show_jmol:
        pdb_ids_found = []
        for block in blocks:
            header = block[0]
            m = re.search(r'[>\s]([A-Za-z0-9]{4})(?:[_|\s]|$)', header)
            if m:
                candidate = m.group(1).upper()
                if candidate not in pdb_ids_found:
                    pdb_ids_found.append(candidate)
        if pdb_id_used and pdb_id_used not in pdb_ids_found:
            pdb_ids_found.insert(0, pdb_id_used)

        print("<div class='result-title'>3D-Struktur (JSmol):</div>")

        if not pdb_ids_found:
            print("<div class='status-msg info'><strong>Kein JSmol-Viewer:</strong> Keine PDB ID im Ergebnis gefunden. Bitte PDB ID als Sequenz-Input verwenden oder PDB ID im FASTA-Header angeben (z.B. <code>&gt;1TIM_A|...</code>).</div>")
        else:

            pred_by_chain = {}
            for block in blocks:
                header = block[0]
                m = re.search(r'[>\s]([A-Za-z0-9]{4})_([A-Za-z0-9]+)', header)
                if not m:
                    continue
                chain = m.group(2)
                ps_line = None
                as_line = None
                for line in block[1:]:
                    if line.startswith("PS "):
                        ps_line = line[3:]
                    elif line.startswith("AS "):
                        as_line = line[3:]
                if ps_line:
                    pred_by_chain[chain] = ps_line

            for idx, pid in enumerate(pdb_ids_found):
                pid_esc = html.escape(pid)
                uid_a   = f"jmol_color_{idx}"
                script_color = (
                    f"load ={pid}; background white; cartoons on; wireframe off; rockets off; "
                    "select helix; color [230,57,70]; "
                    "select sheet; color [244,162,97]; "
                    "select coil; color [173,181,189]; "
                    "spin off; zoom 80;"
                )

                show_two_col = (ref_mode == "pdbid")

                if show_two_col:
                    print("<div class='jmol-two-col'>")

                print(f"""<div class="jmol-wrap">
      <h4>Echte Struktur (H / E / C)</h4>
      <div id="{uid_a}"></div>
      <div class="jmol-legend">
        <span><span class="leg-dot" style="background:#e63946"></span> Helix</span>
        <span><span class="leg-dot" style="background:#f4a261"></span> Strang</span>
        <span><span class="leg-dot" style="background:#adb5bd"></span> Coil</span>
      </div>
    </div>
    <script>
      (function() {{
        var info = {{ width:380, height:380, use:"HTML5", color:"#f8fcfb",
          j2sPath:"https://chemapps.stolaf.edu/jmol/jsmol/j2s", script:"{script_color}" }};
        var app = Jmol.getApplet("{uid_a}_app", info);
        document.getElementById("{uid_a}").innerHTML = Jmol.getAppletHtml(app);
      }})();
    </script>""")

                if show_two_col:
                    uid_b = f"jmol_gray_{idx}"

                    ref_pdb_for_jmol = (fs.getvalue("ref_pdb_id") or "").strip().upper()
                    wrong_residues_by_chain = {}

                    if ref_pdb_for_jmol and ref_pdb_for_jmol == pid:
                        tmp_sec2, resnum_order2, err2 = fetch_seclib_from_pdb(ref_pdb_for_jmol)
                        if not err2 and tmp_sec2:
                            ref_ss_by_chain = {}
                            ref_resnums_by_chain = {}
                            with open(tmp_sec2) as f2:
                                cur_chain = None
                                for line2 in f2:
                                    line2 = line2.strip()
                                    if line2.startswith(">"):
                                        parts2 = line2[1:].strip().split("_")
                                        cur_chain = parts2[-1] if len(parts2) > 1 else "A"
                                    elif line2.startswith("PS ") and cur_chain:
                                        ref_ss_by_chain[cur_chain] = line2[3:]
                            os.unlink(tmp_sec2)


                            wrong_select_parts = []
                            for chain, pred_ss in pred_by_chain.items():
                                ref_ss = ref_ss_by_chain.get(chain, "")
                                res_nums = (resnum_order2 or {}).get(chain, [])
                                min_len = min(len(pred_ss), len(ref_ss), len(res_nums)) if res_nums else min(len(pred_ss), len(ref_ss))
                                for pos in range(min_len):
                                    if pred_ss[pos] != ref_ss[pos]:
                                        pdb_resnum = res_nums[pos] if pos < len(res_nums) else pos + 1
                                        wrong_select_parts.append(f"{pdb_resnum}:{chain}")

                            if wrong_select_parts:
                                wrong_sel = " or ".join(wrong_select_parts)
                                script_compare = (
                                    f"load ={pid}; background white; cartoons on; wireframe off; rockets off; "
                                    "select all; color [180,180,180]; "
                                    f"select {wrong_sel}; color [155,93,229]; "
                                    "spin off; zoom 80;"
                                )
                            else:
                                script_compare = (
                                    f"load ={pid}; background white; cartoons on; wireframe off; rockets off; "
                                    "select all; color [180,180,180]; "
                                    "spin off; zoom 80;"
                                )
                        else:
                            script_compare = (
                                f"load ={pid}; background white; cartoons on; wireframe off; rockets off; "
                                "select all; color [180,180,180]; "
                                "spin off; zoom 80;"
                            )
                    else:
                        script_compare = (
                            f"load ={pid}; background white; cartoons on; wireframe off; rockets off; "
                            "select all; color [180,180,180]; "
                            "spin off; zoom 80;"
                        )

                    print(f"""<div class="jmol-wrap">
      <h4>Vorhersage (grau = korrekt, lila = falsch)</h4>
      <div id="{uid_b}"></div>
      <div class="jmol-legend">
        <span><span class="leg-dot" style="background:#b0b0b0"></span> Korrekte Prediction</span>
        <span><span class="leg-dot" style="background:#9b5de5"></span> Falsche Prediction</span>
      </div>
    </div>
    <script>
      (function() {{
        var info = {{ width:380, height:380, use:"HTML5", color:"#f8fcfb",
          j2sPath:"https://chemapps.stolaf.edu/jmol/jsmol/j2s", script:"{script_compare}" }};
        var app = Jmol.getApplet("{uid_b}_app", info);
        document.getElementById("{uid_b}").innerHTML = Jmol.getAppletHtml(app);
      }})();
    </script>""")
                    print("</div>")

    for t in [tmp_seq, tmp_model, tmp_ref]:
        if t and os.path.exists(t): os.unlink(t)

    print("</div>")
    print(HTML_FOOTER)


if __name__ == "__main__":
    main()