#!/usr/bin/env python

import random
import sys


seq_to_randomize = sys.argv[1]
try: random_seed = sys.argv[2]
except IndexError as warning:
    print("No random seed provided, using 1 as seed")
    random_seed = 1

random.seed(random_seed)


def randomize(dna_seq): 
    dna_seq = dna_seq.upper()
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
        '---':'', 
    } 
    inv_map = {}
    for k, v in table.items():
        inv_map[v] = inv_map.get(v, []) + [k]
    protein = "" 
    try: assert len(dna_seq) %3 == 0, "WARNING: dna sequence length is not a multiple of 3. Truncating extra nucleotides."
    except AssertionError as warning:
        print(warning)
        excess = len(dna_seq) % 3
        dna_seq = dna_seq[:(-1*excess)]
    for i in range(0, len(dna_seq), 3): 
        codon = dna_seq[i:i + 3] 
        protein+= table[codon]
    random_seq = ''.join([random.choice(inv_map[i]) for i in list(protein)])
    return random_seq

print("Randomized nucleotide sequence with same AA sequence:")
print(randomize(seq_to_randomize))