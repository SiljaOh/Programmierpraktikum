#!/bin/bash

# search_pid.sh sucht auf mehreren Rechnern nach einem Prozess
# benutzt isrunning.sh aus Aufgabe 14
# wenn Prozess laeuft -> Hinweise in /tmp sammeln

# 1) Parameter einlesen
while getopts "h:u:p:o:" opt; do
  case "$opt" in
    h) hostfile="$OPTARG" ;;   # Datei mit Hostnamen
    u) user="$OPTARG" ;;       # Prozess-Owner
    p) process="$OPTARG" ;;    # Prozessname
    o) out="$OPTARG" ;;        # Outputdatei
    *) exit 2 ;;
  esac
done

# 2) Parameter pruefen
if [ -z "$hostfile" ] || [ -z "$user" ] || [ -z "$process" ] || [ -z "$out" ]; then
  echo "Usage: $0 -h hostnames -u user -p process -o out" >&2
  exit 2
fi

# 3) Outputdatei neu anlegen
: > "$out"

# 4) Hostliste zeilenweise lesen
while IFS= read -r hostname; do

  # leere Zeilen ueberspringen
  [ -z "$hostname" ] && continue

  # Kommentarzeilen ueberspringen
  case "$hostname" in
    \#*) continue ;;
  esac

  # 5) isrunning.sh remote starten
  # Aufgabe 14: Exit 1 = Prozess laeuft
  ssh -p 12 "fiocco@${hostname}" 'bash -s' < ./isrunning.sh "$process" "$user"
  rc=$?

  # 6) wenn Prozess laeuft -> Hinweise sammeln
  if [ "$rc" -eq 1 ]; then
    echo "Process found on ${hostname}" >> "$out"
    echo "Hints in /tmp:" >> "$out"

    # zuerst nach Prozessnamen suchen
    ssh -p 12 "fiocco@${hostname}" \
      "ls -la /tmp | grep -F \"$process\"" >> "$out"

    echo >> "$out"
  fi

done < "$hostfile"