#!/bin/bash

# is_running.sh prüft, ob ein bestimmter Prozess eines Users aktiv ist
# Exit-Code (d.h. Rückgabewert): 0 = Prozess läuft nicht, 1 = Prozess läuft
# Exit 0: Prozess läuft nicht
# Exit 1: Prozess läuft
# Exit 2: falsche Parameter
# zuerst Argumente prüfen: sind Prozessname (process name) und User (user) angegeben?
# Anzahl ($#) der Argument weniger als 2

if [ $# -ne 2 ]; then
    exit 2
fi

# Prozesse suchen
proc="$1"
usr="$2"

# Exit-Code printen
if pgrep -x -u "$usr" "$proc" >/dev/null 2>&1; then
  exit 1
else
  exit 0
fi