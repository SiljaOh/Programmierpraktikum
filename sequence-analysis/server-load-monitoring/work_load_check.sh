#!/bin/bash

# work_load_check.sh gibt die Systemlast bzw. 1 Minuten Load eines Remote-Rechn>
# mindestenst 1 Namen eines Rechners eingeben.

# Hostname?
if [ $# -ne 1 ]; then
    exit 2
fi

# host = Rechnername des Zielservers bzw. Zielrechner
host="$1"

# In Linux: <1min> <5min> <15min> <running/total> <last_pid>
# loadavg (load average) = durchschnittliche Systemlast

ssh -p 12 "$host" "cut -d ' ' -f1 /proc/loadavg"
