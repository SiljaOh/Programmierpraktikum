#!/bin/bash

# 1) Parameter einlesen
while getopts "h:u:p:o:" opt; do
  case "$opt" in
    h) hostfile="$OPTARG" ;;
    u) proc_user="$OPTARG" ;;   
    p) process="$OPTARG" ;;     
    o) out="$OPTARG" ;;         # Outputdatei
    *) exit 2 ;;
  esac
done

# 2) Pflichtparameter prüfen
if [ -z "$hostfile" ] || [ -z "$proc_user" ] || [ -z "$process" ] || [ -z "$out" ]; then
  echo "Usage: $0 -h hostfile -u process_owner -p process -o out" >&2
  exit 2
fi

# 3) Outputdatei leeren (neu)
: > "$out"

# 4) SSH Login User festlegen
# -u ist der Prozess-Owner (nicht user)  
ssh_user="fiocco"    

# 5) Hostliste (pro Zeile)
while IFS= read -r hostname; do

  [ -z "$hostname" ] && continue

  case "$hostname" in
    \#*) continue ;;
  esac

# 6) isrunning.sh remote ausführen
# Exit-Code 0 bedeutet "Prozess läuft"
ssh -p 12 "${ssh_user}@${hostname}" 'bash -s' < ./isrunning.sh "$process" "$proc_user"
rc=$?

# 7) wenn Prozess läuft, Host (Rechner) + /tmp Hinweise (Dateien des Programms) sammeln
if [ "$rc" -eq 0 ]; then
  {
    echo "Process found on ${hostname}"
    echo "Hints in /tmp:"
    ssh -p 12 "${ssh_user}@${hostname}" \
      "ls -la /tmp | grep -F \"$proc_user\" | grep -F \"$process\""
    echo
  } >> "$out"
fi

done < "$hostfile"