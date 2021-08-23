#!/usr/bin/env python

import itertools


firstseq = 'abcdefghijklmnopqrstuvwx'
secondseq = 'abcdefghijklmnopqrstuvwx'.upper()

class fastq_read:
    def __init__(self, lines):
        self.id, self.seq, self.line3, self.quality = lines

def read_by_lines(filename):
    with open(filename) as f:
        for lines in itertools.zip_longest(*[f]*4):
            yield(fastq_read([x.strip() for x in lines]))

def iterate_chimera_database(seq1, seq2):
    for i in range(0, len(seq1), 3):
        yield seq1[0:i] + seq2[i:]
    for i in range(0, len(seq1), 3):
        yield seq2[0:i] + seq1[i:]

def get_distance(string1, string2):
    len1 = len(string1)
    len2 = len(string2)
    assert len1 == len2, "string lengths do not match!"    
    score = 0
    for i in range(len1):
        if string1[i] == string2[i]:
            score += 1
    return(score)

def find_closest_chimera(read):
    for i in iterate_chimera_database(firstseq, secondseq):
        get_distance(read.seq, i)

for read in read_by_lines('example.fastq'):
    find_closest_chimera(read)






