#!/usr/bin/env python3

import sys
filename = sys.argv[1]

# the second Fasta file is the unmodified sequence
# being compared to this one to determine % identity
comparativeFastaFilename = sys.argv[2]

def wrap_fasta(seq, line_length=60):
    return '\n'.join([seq[i:i+line_length] for i in range(0, len(seq), line_length)])

def levenshtein(string1, string2):
    assert len(string1) == len(string2), "Strings unequal length!"
    seqLen = float(len(string1))
    seqDist = 0
    for i in zip(string1, string2):
        if i[0] != i[1]:
            seqDist += 1
    seqIdentity = "{:.2f}".format(100 * (seqLen - seqDist) / seqLen)
    return([seqLen, seqDist, seqIdentity])

def printFasta(header, seq):
    print(">%s\n" % (header) + wrap_fasta(seq))

def translate(dna_seq): 
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
        '---':'---', 
    } 
    protein ="" 
    try: assert len(dna_seq) %3 == 0, "WARNING: dna sequence length is not a multiple of 3. Truncating extra nucleotides."
    except AssertionError as warning:
        print(warning)
        excess = len(dna_seq) % 3
        dna_seq = dna_seq[:(-1*excess)]
    for i in range(0, len(dna_seq), 3): 
        codon = dna_seq[i:i + 3] 
        protein+= table[codon]     
    return protein 


class fasta:
    help = 'stores type, header, and sequence information for FASTA files'
    def __init__(self, filename):
        with open(filename, 'r') as infile:
            text = infile.readlines()
        self.header = text[0].replace(">","").strip()
        self.seq = ''.join(x.rstrip() for x in text[1:]).upper()
        codons = [self.seq[x:(x+3)] for x in range(0,len(self.seq),3)]
        homopolymerDict = {
        'A|AAA|A' : 'A|AAG|A',
        '|AAA|AAA|' : '|AAG|AAA|',
        '|AAA|AA' : '|AAG|AA',
        'AA|AAA|' : 'AA|AAG|',
        'T|TTT|T' : 'T|TTC|T',
        '|TTT|TTT|' : '|TTC|TTT|',
        '|TTT|TT' : '|TTC|TT',
        'TT|TTT|' : 'TT|TTC|',
        'G|GGG|G' : 'G|GGC|G',
        '|GGG|GGG|' : '|GGC|GGT|',
        '|GGG|GG' : '|GGC|GG',
        'GG|GGG|' : 'GG|GGT|',
        'C|CCC|C' : 'C|CCT|C',
        '|CCC|CCC|' : '|CCA|CCG|',
        '|CCC|CC' : '|CCA|CC',
        'CC|CCC|' : 'CC|CCG|'
        }
        splitCodons = '|'.join(codons)
        for _ in range(3):
            for x in homopolymerDict:
                splitCodons = splitCodons.replace(x, homopolymerDict[x])
        self.optimizedSeq = splitCodons.replace("|", "")
        self.seqProtein = translate(self.seq)
        self.optimizedSeqProtein = translate(self.optimizedSeq)



infileFasta = fasta(filename)

comparativeFasta = fasta(comparativeFastaFilename)

seqName = filename.replace(".fasta","")
seqLen, seqDist, seqIdentity = levenshtein(infileFasta.optimizedSeq, comparativeFasta.seq)

fastaHeader = "%s|length=%s|identity=%s|dist=%s" % (seqName, int(seqLen), seqIdentity, seqDist)

printFasta(fastaHeader, infileFasta.optimizedSeq)