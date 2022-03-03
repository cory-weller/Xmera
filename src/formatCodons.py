#!/usr/bin/env python

# reformats codon table, e.g. from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

import sys

codonInput = sys.argv[1]

with open(codonInput, "r") as infile:
    inputText = infile.read().replace("(","").replace(")","").split()

entries = [inputText[x:(x+5)] for x in range(0, len(inputText), 5)]

print('\t'.join(['triplet', 'amino_acid', 'fraction', 'frequency_per_thousand', 'number']))
for row in entries:
    print('\t'.join([str(x) for x in row]))