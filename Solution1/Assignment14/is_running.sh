#!/bin/bash

# is_running.sh prüft, ob ein bestimmter Prozess eines Users aktiv ist
# Exit-Code (d.h. Rückgabewert): 0 = Prozess läuft nicht, 1 = Prozess läuft


# zuerst Argumente prüfen: sind Prozessname (process name) und User (user) angegeben?
# Anzahl ($#) der Argument weniger als 2
if [ $# -lt 2 ]; then
    echo "More arguments are needed."
    exit 2
fi

# Prozesse suchen
processFound=$(pgrep -x -c -u "$2" "$1")

# Exit-Code printen
if [ "$processFound" -ne 0 ]; then
    exit 1
else
    exit 0
fi