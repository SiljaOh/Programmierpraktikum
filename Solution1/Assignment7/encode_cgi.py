#!/usr/bin/python3
import mycgi
import sys
import re
import html
import urllib.request

sys.stderr = sys.stdout

print("Content-type: text/html\n\n")


#
#
# aus psscan.py
def parse_fasta(text):
    sequences = {}
    current_header = None
    current_seq = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current_header:
                sequences[current_header] = "".join(current_seq)
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line)
    if current_header:
        sequences[current_header] = "".join(current_seq)
    return sequences

# aus psscan.py
def prosite_to_regex(pattern):
    regex = pattern.strip().rstrip(".")
    regex = regex.replace("<", "^").replace(">", "$")
    regex = re.sub(r'\{([A-Za-z]+)\}', r'[^\1]', regex)
    regex = regex.replace("x", ".")
    regex = regex.replace("-", "")
    regex = re.sub(r'\((\d+(?:,\d+)?)\)', r'{\1}', regex)
    return regex

# aus psscan.py
def pattern_from_web(prosite_id):
    url = "https://ftp.expasy.org/databases/prosite/prosite.dat"
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            found = False
            for line in response:
                line = line.decode("utf-8").rstrip("\n")
                if line.startswith("AC") and prosite_id in line:
                    found = True
                if found and line.startswith("PA"):
                    return line[2:].strip().rstrip(".")
    except Exception as e:
        print(f"Fehler beim Abrufen: {e}")
    return None

#aus psscan.py
def scan_seq(sequences, regex):
    results = []
    search_regex = f"(?=({regex}))"
    for header, sequence in sequences.items():
        for match in re.finditer(search_regex, sequence):
            results.append((header, match.start() + 1, match.group(1)))
    return results


#
#
#


form = mycgi.Form()

fasta_input   = form.getvalue("fasta")   or ""
pattern_input = form.getvalue("pattern") or ""
web_input     = form.getvalue("web")     or ""


def h(s):
    return html.escape(str(s))



#HTML
print(f"""
<h1>PROSITE Pattern Scanner</h1>

<form method="POST" action="">

  <label>FASTA-Sequenz:</label><br>
  <input type="text" name="fasta" value="{h(fasta_input)}"
     placeholder="z.B. C-x(2,4)-C-x(3)-[LIVMFY]" size="40">
  <br><br>

  <label>Pattern (PROSITE-Format):</label><br>
  <input type="text" name="pattern" value="{h(pattern_input)}"
         placeholder="z.B. C-x(2,4)-C-x(3)-[LIVMFY]" size="40">
  <br><br>

  <b>– ODER –</b><br><br>

  <label>PROSITE-ID (wird online abgerufen):</label><br>
  <input type="text" name="web" value="{h(web_input)}"
         placeholder="z.B. PS00001" size="20">
  <br><br>

  <input type="submit" value="Scannen">
</form>
""")



#Verabreitung mit CGI
if fasta_input:
    current_pattern = None
    pattern_label   = ""

    if pattern_input:
        current_pattern = pattern_input
        pattern_label   = "Manuelles Pattern"
    elif web_input:
        current_pattern = pattern_from_web(web_input)
        if current_pattern:
            pattern_label = f"Pattern für {h(web_input)}"
        else:
            print(f'<p style="color:red">Pattern für ID „{h(web_input)}" nicht gefunden.</p>')

    if current_pattern:
        print(f'<p>Verwendetes Pattern ({pattern_label}): <code>{h(current_pattern)}</code></p>')
        try:
            regex_str = prosite_to_regex(current_pattern)
            sequences = parse_fasta(fasta_input)
            results   = scan_seq(sequences, regex_str)

            if results:
                print("<pre>")
                for header, pos, matched in results:
                    print(f"{h(header)}\t{pos}\t{h(matched)}")
                print("</pre>")
            else:
                print("<p>Keine Treffer gefunden.</p>")

        except re.error as e:
            print(f'<p style="color:red">Ungültiger Regex: {h(e)}</p>')
        except Exception as e:
            print(f'<p style="color:red">Fehler: {h(e)}</p>')

    elif not pattern_input and not web_input:
        print('<p style="color:red">Bitte Pattern oder PROSITE-ID angeben.</p>')

