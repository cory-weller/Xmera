#!/usr/bin/env python

import sys
import random

fn = sys.argv[1]

def import_fasta(filename):
    with open(filename, 'r') as infile:
        text = infile.readlines()
    text = [x.rstrip() for x in text]
    header = text[0]
    sequence = ''.join(text[1:])
    return([header, sequence])


def permute_fasta(fasta):
    header = fasta[0]
    sequence = list(fasta[1].replace("*",""))
    random.shuffle(sequence)
    sequence.append("*")
    sequence = ''.join(sequence)
    return([header, sequence])


def fold_fasta(fasta):
    header = fasta[0]
    sequence = fasta[1]
    sequence = '\n'.join([sequence[i:i+60] for i in range(0, len(sequence), 60)])
    return(header + '\n' + sequence)

print(fold_fasta(permute_fasta(import_fasta(fn))))
