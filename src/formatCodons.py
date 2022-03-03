#!/usr/bin/env python

# reformats the codon table codonTable.txt, from
# the Codon Usage Database: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

with open("codonTable.txt", "r") as infile:
    inputText = infile.read().replace("(","").replace(")","").split()

entries = [inputText[x:(x+5)] for x in range(0, len(inputText), 5)]

for row in entries:
    print('\t'.join([str(x) for x in row]))