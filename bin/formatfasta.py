#!/usr/bin/env python

import sys

filename = sys.argv[1]

with open(filename, 'r') as infile:
    text = [x.strip() for x in infile.readlines()]

header = text[0]

dna = ''.join(text[1:]).replace("-","")

print(header + '\n' + dna)