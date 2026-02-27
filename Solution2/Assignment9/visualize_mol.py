#!/usr/bin/python3
import argparse
import subprocess
import webbrowser

#JMOL räumliche darstellung Moleküle, drehen,zoomen, messen, ...

parser = argparse.ArgumentParser()
parser.add_argument('--id', required=True)
parser.add_argument('--output', type=argparse.FileType('w'))
parser.add_argument('--html', action="store_true")
parser.add_argument('--colourized', action='store_true')
args = parser.parse_args()


if len(args.id) != 4:
    print("Ungültige ID")
else:
    pdb_id = args.id


if args.output:
    output_name = args.output
else:
    output_name= f"{pdb_id}.png"


jmol_commands= f"load={pdb_id}; cartoon on;"
if args.colourized:
    jmol_commands+= "color structure; "
jmol_commands_web = jmol_commands
jmol_commands+= f"write image png '{output_name}'; "

#-n damit im hintergrund läuft, -j führt die commands aus
subprocess.run(["java", "-jar", "Jmol.jar", "-n", "-j", jmol_commands])

#html für die website machen
if args.html:
    #parallel open! muss ja nebenbei öfnnen, -m nmmt modul http.server und startet, 8000 srandard
    subprocess.Popen(["python", "-m", "http.server", "8000"])

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>JSmol: {pdb_id}</title>
        <script type="text/javascript" src="jsmol/JSmol.min.js"></script>
        
        <script type="text/javascript">
        var Info = {{
            width: 500,
            height: 500,
            use: "HTML5",
            j2sPath: "jsmol/j2s",
            script: "{jmol_commands_web}"
        }};
        </script>
    </head>
    <body>
        <script type="text/javascript">
            jmolApplet0 = Jmol.getApplet("jmolApplet0", Info);
        </script>
    </body>
    </html>
    """

    # HTML-Datei kurzzeitig speichern
    args.output.write(html_content)

    # Die Datei automatisch im Standardbrowser öffnen
    webbrowser.open('http://localhost:8000/jsmol_website.html')

#schliessen?
#was wenn - als output
