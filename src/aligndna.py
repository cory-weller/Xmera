#!/usr/bin/env python

# walk along seq and print aligned nucleotide

import sys

filestem = sys.argv[1]

dna_filename = filestem + ".cds.fasta"
pep_filename = filestem + ".align.pep"

with open(dna_filename, 'r') as infile:
    dna = ''.join([x.strip() for x in infile.readlines()[1:]])

with open(pep_filename, 'r') as infile:
    pep = ''.join([x.strip() for x in infile.readlines()[1:]])

dna
pep

output = ''

for position in pep:
    if position == "-":
        output += "---"
    else:
        output += dna[:3]
        dna = dna[3:]

print(output)