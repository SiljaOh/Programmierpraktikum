#!/usr/bin/python3
import argparse
import os
import subprocess
import urllib.request

#für Liste von PDBs
parser = argparse.ArgumentParser()
parser.add_argument('--id', nargs='+', required=True)
parser.add_argument('--output')
args = parser.parse_args()


#base_ids = ["1MBN", "1BAF", "1TIM", "256b"]
#ids = base_ids + args.id
ids = args.id

if args.output:
    output_dir = args.output

for id in ids:
    if len(id) != 4:
        print("Ungültige ID")
        continue

    atoms = []
    sekund_s = set()
    #pdb für größe der box etc beuntzen
    url = f"https://files.rcsb.org/download/{id}.pdb"
    response = urllib.request.urlopen(url)
    data = response.read().decode('utf-8')
    #splitlines wg auch leerzeichen
    lines = data.splitlines()
    #atome finden
    for l in lines:
        if l.startswith("ATOM"):
            atoms.append(l)

        if l.startswith(("HELIX")):
            try:
                start = int(l[21:25].strip())
                end = int(l[33:37].strip())
                # Alle IDs in diesem Bereich hinzufügen
                for r in range(start, end + 1):
                    sekund_s.add(str(r))
            except (ValueError, IndexError):
                continue

        if l.startswith(("SHEET")):
            try:
                start = int(l[22:26].strip())
                end = int(l[33:37].strip())
                # Alle IDs in diesem Bereich hinzufügen
                for r in range(start, end + 1):
                    sekund_s.add(str(r))
            except (ValueError, IndexError):
                continue

    coordinates = []
    ca = []
    cb = []
    all_as = set()

    for a in atoms:
        #x,y,z...
        new_coor = [float(a[30:38]), float(a[38:46]), float(a[46:54])]
        coordinates.append(new_coor)

        if "CA" in a[12:16]:
            ca.append(new_coor)
        elif "CB" in a[12:16]:
            cb.append(new_coor)
        #ID, wie viele gibt es insgesamt von AS, unique
        all_as.add(a[22:26].strip())


    mins = [float('inf'), float('inf'), float('inf')]
    maxs = [float('-inf'), float('-inf'), float('-inf')]

    for c in coordinates:
        for i in range(3):
            #i = 0 (X), 1 (Y), 2 (Z)
            if c[i] < mins[i]:
                mins[i] = c[i]
            if c[i] > maxs[i]:
                maxs[i] = c[i]

    #Kantenlängen berechnen
    dx = maxs[0] - mins[0]
    dy =  maxs[1] - mins[1]
    dz = maxs[2] - mins[2]

    volume =  dx * dy * dz

    count_ss = 0
    for res in all_as:
        if res in sekund_s:
            count_ss += 1

    count_total = len(all_as)
    as_anteil = (count_ss / count_total)

    if len(ca) > 0:
        distance_ca = ((ca[0][0] - ca[-1][0]) ** 2 + (ca[0][1] - ca[-1][1]) ** 2 + (ca[0][2] - ca[-1][2]) ** 2) ** 0.5
    else:
        distance_ca = 0.0

    if len(cb) > 0:
        distance_cb = ((cb[0][0] - cb[-1][0]) ** 2 + (cb[0][1] - cb[-1][1]) ** 2 + (cb[0][2] - cb[-1][2]) ** 2) ** 0.5
    else:
        distance_cb = 0.0

    print(f"{id}\tAnteil AS in Sekundaerstruktur\t{as_anteil:.4f}")
    print(f"{id}\tAbstand C_alpha\t{distance_ca:.4f}")
    print(f"{id}\tAbstand C_beta\t{distance_cb:.4f}")
    print(f"{id}\tX-Groesse\t{dx:.4f}")
    print(f"{id}\tY-Groesse\t{dy:.4f}")
    print(f"{id}\tZ-Groesse\t{dz:.4f}")
    print(f"{id}\tVolumen\t{volume:.4f}")

    jmol_cmd = f" load ={id}; cartoon on; "

    if args.output:
        jmol_cmd += f" write image png '{output_dir}/{id}.png'; "
    else:
        jmol_cmd += f" write image png '{id}.png';"

    if os.path.exists("Jmol.jar"):
        # -n damit im hintergrund läuft, -j führt die commands aus
        subprocess.run(["java", "-jar", "Jmol.jar", "-n", "-j", jmol_cmd])
    else:
        pass
