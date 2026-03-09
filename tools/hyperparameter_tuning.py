#!/usr/bin/python3

import subprocess
import re
import csv
import itertools
import random
import shutil
import tempfile
from pathlib import Path

seclib_file     = r"C:\Users\annan\OneDrive - TUM\Desktop\propaBlock\CB513DSSP.db"
train_jar       = r"C:\Users\annan\OneDrive - TUM\Desktop\propra\src\progprakt-17\out\artifacts\train\train.jar"
predict_jar     = r"C:\Users\annan\OneDrive - TUM\Desktop\propra\src\progprakt-17\out\artifacts\predict\predict.jar"
validateGor_jar = r"C:\Users\annan\OneDrive - TUM\Desktop\propra\src\progprakt-17\out\artifacts\validateGor\validateGor.jar"
train_ratio     = 0.8
seed            = 42

pseudos     = [0.01, 0.1, 0.5, 0.8, 1.0, 2.0]
minH_values = [3, 4,5, 6]
minE_values = [2, 3,4]

out_file = "tuning_results.tsv"

text = Path(seclib_file).read_text(encoding="utf-8", errors="replace").strip()
text = text.replace('\r\n', '\n').replace('\r', '\n')
entries = [e.strip() for e in text.split("\n\n") if e.strip()]

random.seed(seed)
random.shuffle(entries)

split      = int(len(entries) * train_ratio)
train_data = entries[:split]
val_data   = entries[split:]

print(f"Train: {len(train_data)} Proteine, Val: {len(val_data)} Proteine\n")

tmpdir = Path(tempfile.mkdtemp())
train_db  = tmpdir / "tune_train.db"
val_db    = tmpdir / "tune_val.db"
val_fasta = tmpdir / "tune_val.fasta"

train_db.write_text("\n\n".join(train_data), encoding="utf-8")
val_db.write_text  ("\n\n".join(val_data),   encoding="utf-8")

def get_id(entry):
    for l in entry.split("\n"):
        if l.startswith(">"): return l[1:].strip()

def get_seq(entry):
    for l in entry.split("\n"):
        if l.startswith("AS "): return l[3:].strip()

fasta_lines = []
for e in val_data:
    fasta_lines.append(">" + get_id(e))
    fasta_lines.append(get_seq(e))
val_fasta.write_text("\n".join(fasta_lines), encoding="utf-8")

def parse_summary(path):
    scores = {}
    if not Path(path).exists():
        return scores
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            if re.match(r'^q3\s*:', s, re.IGNORECASE):
                m = re.search(r'Mean:\s*([\d.]+)', s)
                if m: scores["Q3"] = float(m.group(1))
            elif re.match(r'^sov\s*:', s, re.IGNORECASE):
                m = re.search(r'Mean:\s*([\d.]+)', s)
                if m: scores["SOV"] = float(m.group(1))
    return scores

def run_combo(method, pseudo, minH, minE):
    model = str(tmpdir / "tune_model.txt")
    pred  = str(tmpdir / "tune_pred.prd")
    summ  = str(tmpdir / "tune_summary.txt")
    det   = str(tmpdir / "tune_detailed.txt")


    r = subprocess.run([
        "java", "-jar", train_jar,
        "--db",     str(train_db),
        "--method", method,
        "--model",  model,
        "--pseudo", str(pseudo),
        "--minH",   str(minH),
        "--minE",   str(minE),
    ], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"train failed:\n{r.stderr.strip()}")


    with open(pred, "w") as pred_out:
        r = subprocess.run([
            "java", "-jar", predict_jar,
            "--model",  model,
            "--seq",    str(val_fasta),
            "--format", "txt",
            "--minH",   str(minH),
            "--minE",   str(minE),
            "--pseudo", str(pseudo),
        ], stdout=pred_out, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"predict failed:\n{r.stderr.strip()}")


    r = subprocess.run([
        "java", "-jar", validateGor_jar,
        "-p", pred, "-r", str(val_db),
        "-s", summ, "-d", det, "-f", "txt",
    ], capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"validate failed:\n{r.stderr.strip()}")

    scores = parse_summary(summ)
    return scores.get("Q3", 0.0), scores.get("SOV", 0.0)

# gor1 pseudo irrelevant ->  nur minH minE tunen
gor1_combos = [("gor1", 0.8, minH, minE)
               for minH, minE in itertools.product(minH_values, minE_values)]

# gor3 gor4 alle drei Parameter
gor34_combos = [
    (method, pseudo, minH, minE)
    for method in ["gor3", "gor4"]
    for pseudo, minH, minE in itertools.product(pseudos, minH_values, minE_values)
]

all_combos = gor1_combos + gor34_combos
total = len(all_combos)
print(f"GOR1: {len(gor1_combos)} Kombos (nur minH/minE)")
print(f"GOR3+GOR4: {len(gor34_combos)} Kombos")
print(f"Gesamt: {total} Kombinationen\n")

# beste Parameter pro Methode
best_per_method = {"gor1": None, "gor3": None, "gor4": None}
best_overall    = None

with open(out_file, 'w', newline='') as tsvfile:
    writer = csv.writer(tsvfile, delimiter='\t')
    writer.writerow(["method", "pseudo", "minH", "minE", "Q3", "SOV"])

    for idx, (method, pseudo, minH, minE) in enumerate(all_combos):
        pseudo_str = "n/a" if method == "gor1" else str(pseudo)
        print(f"[{idx+1}/{total}] {method} pseudo={pseudo_str} minH={minH} minE={minE}", end=" ... ", flush=True)

        try:
            q3, sov = run_combo(method, pseudo, minH, minE)
            print(f"Q3={q3:.1f}%  SOV={sov:.1f}%")
            writer.writerow([method, pseudo_str, minH, minE, f"{q3:.2f}", f"{sov:.2f}"])
            tsvfile.flush()

            params = (method, pseudo_str, minH, minE, q3, sov)
            if best_per_method[method] is None or q3 > best_per_method[method][4]:
                best_per_method[method] = params
            if best_overall is None or q3 > best_overall[4]:
                best_overall = params

        except (subprocess.CalledProcessError, RuntimeError) as e:
            print(f"ERROR: {e}")
            writer.writerow([method, pseudo_str, minH, minE, "ERROR", "ERROR"])

shutil.rmtree(tmpdir, ignore_errors=True)

print(f"\n{'='*55}")
for method in ["gor1", "gor3", "gor4"]:
    p = best_per_method[method]
    print(f"\nBeste Parameter [{method.upper()}]:")
    if p:
        print(f"  Pseudo  : {p[1]}")
        print(f"  minH    : {p[2]}")
        print(f"  minE    : {p[3]}")
        print(f"  Q3      : {p[4]:.2f}%")
        print(f"  SOV     : {p[5]:.2f}%")
    else:
        print("  Keine erfolgreichen Kombinationen!")

print(f"\n{'='*55}")
print(f"\nBeste Parameter [OVERALL]:")
if best_overall:
    print(f"  Methode : {best_overall[0]}")
    print(f"  Pseudo  : {best_overall[1]}")
    print(f"  minH    : {best_overall[2]}")
    print(f"  minE    : {best_overall[3]}")
    print(f"  Q3      : {best_overall[4]:.2f}%")
    print(f"  SOV     : {best_overall[5]:.2f}%")
print(f"\nAlle Ergebnisse in: {out_file}")