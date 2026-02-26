#!/usr/bin/python3

def parse_csv(filepath):
    data_list = []

    with open(filepath, 'r') as f:
        header = f.readline().strip().split(',')

        for line in f:
            line = line.strip()
            if not line:
                continue
            values = line.split(',')
            row_dict = dict(zip(header, values))
            data_list.append(row_dict)

    return data_list