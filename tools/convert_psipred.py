#!/usr/bin/python3

# Konvertiert CB513.psipred-> validateGor-Format (einzeilig, PS)
import sys

input_file  = sys.argv[1] if len(sys.argv) > 1 else 'CB513.psipred'
output_file = sys.argv[2] if len(sys.argv) > 2 else 'psipred_prediction.txt'

text = open(input_file).read().strip()
entries = [e.strip() for e in text.split('\n\n') if e.strip()]

out = []
for e in entries:
    id_ = ''
    as_ = ''
    ss_ = ''

    mode = None
    for line in e.split('\n'):
        if line.startswith('> '):
            id_ = line[2:].strip()
        elif line.startswith('>'):
            id_ = line[1:].strip()
        elif line.startswith('AS '):
            mode = 'AS'; as_ += line[3:].strip()
        elif line.startswith('SS '):
            mode = 'SS'; ss_ += line[3:].strip()
        elif mode == 'AS':
            as_ += line.strip()
        elif mode == 'SS':
            ss_ += line.strip()
    if id_ and as_ and ss_:
        out.append(f'> {id_}\nAS {as_}\nPS {ss_}')

open(output_file, 'w').write('\n\n'.join(out))
