#!/usr/bin/env python3

import sys

RTfastaFilename = sys.argv[1]
N = int(sys.argv[2])

class RTfasta:
    def __init__(self, filename):
        self.RTs = []
        with open(filename, 'r') as infile:
            fileText = infile.read().split(">")[1:]
            for entry in fileText:
                splitEntry = entry.split("\n")
                header = splitEntry[0].strip()
                seq = ''.join(splitEntry[1:]).strip()
                self.RTs.append((header, seq))


def importSkpp():
    skppF = {}
    skppR = {}
    with open('Xmera/skpp15/skpp15-forward.fasta', 'r') as infile:
        primersF = infile.read().split(">")[1:]
        for primer in primersF:
            splitPrimer = primer.split("\n")
            header = splitPrimer[0].strip()
            seq = ''.join(splitPrimer[1:]).strip()
            N = int(header.split("-")[1])
            skppF[N] = seq
    with open('Xmera/skpp15/skpp15-reverse.fasta', 'r') as infile:
        primersF = infile.read().split(">")[1:]
        for primer in primersF:
            splitPrimer = primer.split("\n")
            header = splitPrimer[0].strip()
            seq = ''.join(splitPrimer[1:]).strip()
            N = int(header.split("-")[1])
            skppR[N] = seq
    return (skppF, skppR)

skppF, skppR = importSkpp()

i = RTfasta(RTfastaFilename)

for header, repairTemplate in i.RTs:
    oligo = skppF[N] + repairTemplate + skppR[N]
    print("%s_skpp-%s\t%s" % (header, N, oligo))