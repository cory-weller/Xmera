#!/usr/bin/env python

## FUNCTIONS

def read_text(filename):
    with open(filename, 'r') as infile:
        rawText = infile.read()
    return rawText

def sampleCodon(df, aminoAcid):
    tmp = df[df['AA'] == aminoAcid]
    return(tmp.sample(1, weights=tmp['freq']).iat[0,0])

def read_fasta(filename):
    with open(filename, 'r') as infile:
        seq = ''.join([x.strip() for x in infile.readlines()[1:]])
        return(seq)

def wrap_fasta(seq, line_length=60):
    return '\n'.join([seq[i:i+line_length] for i in range(0, len(seq), line_length)])


def translate(dna_seq): 
    codonTable = { 
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
    } 
    dna_seq = dna_seq.upper()
    protein ="" 
    try: assert len(dna_seq) %3 == 0, "WARNING: dna sequence length is not a multiple of 3. Truncating extra nucleotides."
    except AssertionError as warning:
        print(warning)
        excess = len(dna_seq) % 3
        dna_seq = dna_seq[:(-1*excess)]
    for i in range(0, len(dna_seq), 3): 
        codon = dna_seq[i:i + 3] 
        protein+= codonTable[codon]     
    return protein

def getDistance(string1, string2):
    dist = 0
    for i in zip(string1, string2):
        if i[0] != i[1]:
            dist +=1
    return(dist)

def split_codons(nucleotide_seq):
    assert (len(nucleotide_seq)%3 == 0), "%s length is not a multiple of three!" % (nucleotide_seq)
    return([nucleotide_seq[x:(x+3)] for x in range(0, len(nucleotide_seq), 3)])


def getSNPs(table):
    snps = {}
    for codon1 in table:
        snps[codon1] = []
        for codon2 in table:
            if (getDistance(codon1, codon2) == 1):
                snps[codon1].append(codon2)
    return(snps)




## CLASSES

# unused
# class fasta:
#     help = 'stores type, header, and sequence information for FASTA files'
#     def __init__(self, filename):
#         with open(filename, 'r') as infile:
#             text = infile.readlines()
#         self.header = text[0].replace(">","").strip()
#         self.seq = ''.join(x.rstrip() for x in text[1:])





# PARSE ARGUMENTS

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('fasta_1_filename', type=str,
    help="File name for first (i.e. upstream) sequence.")
parser.add_argument('fasta_2_filename', type=str,
    help="File name for second (i.e. downstream) sequence.")
parser.add_argument('codon_table', type=str,
    help='''File name for tab-delimited 4-column codon table.
        Should include no header, with columns representing
        triplet, amino-acid code, relative frequency,
        fraction per 1000, and count. See README for help
        finding and generating this file.''')
parser.add_argument("-f", "--filter-freq",
    type=float,
    help='''Minimum required fraction for including a codon.
    Triplets with (amino-acid specific) frequency
    below this provided value will be excluded.
    Default: 0.10 (i.e. include all triplets with 10%%
    or greater represention for its amino acid).''')

args = parser.parse_args()

print(args.filter_freq)

#python3.8 mut.py fasta1 fasta2 condontable 

# IMPORT MODULES
import pandas as pd

# RUN CODE

## 
codons = pd.read_csv(   "codonTableCleaned.tsv", 
             sep = "\t",
             header = None,
            )

codons.columns = ["triplet", "AA", "freq", "fkp", "count"]


putativeCodonsTable = codons[codons.freq >= args.filter_freq]
putativeCodonsTable['triplet'] = putativeCodonsTable['triplet'].str.replace("U","T")


assert putativeCodonsTable['AA'].nunique() == 21, (
    "Codon frequency filter of %s eliminates some amino acids. Try a lower filter." % (args.filter_freq)
)

allAminoAcids = list(putativeCodonsTable['AA'].unique())

# retrieve rows with maximum frequency, grouped by amino acid
abundantCodons = putativeCodonsTable.loc[putativeCodonsTable.groupby('AA')['count'].idxmax()]


# sampleCodon(putativeCodonsTable, "F")

fasta_1 = read_fasta(args.fasta_1_filename)
fasta_2 = read_fasta(args.fasta_2_filename)

# Assert sequences are identical length
assert (len(fasta_1) == len(fasta_2)), "Sequences not same length!"

# Assert amino acid sequences are identical
aa_1 = translate(fasta_1)
aa_2 = translate(fasta_2)
assert (aa_1 == aa_2), "Amino acid sequences not identical!"

# Load upstream and downstream seq

with open('upstream.dna','r') as infile:
    upstream = infile.read().strip()

with open('downstream.dna','r') as infile:
    downstream = infile.read().strip()

codons_1 = split_codons(fasta_1)
codons_2 = split_codons(fasta_2)

# create mutated codon dictionary
# inserted_codons = abundantCodons[["AA","triplet"]].set_index('AA').to_dict()

inserted_codons = list(abundantCodons['triplet'])

## calculate "single-SNP" mutations 
codonTable = { 
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
} 
singleSNPcodons = {}
for codon1 in codonTable: # initialize codon dictionary
    if codon1 not in singleSNPcodons: # add codon if not yet in dictionary
        singleSNPcodons[codon1] = []
    for codon2 in codonTable:
        codon2_AA = codonTable[codon2]
        if getDistance(codon1, codon2) == 1:  # if distance between codons is one SNP
            if codon2_AA not in singleSNPcodons[codon1]:   
                singleSNPcodons[codon1].append(codon2_AA)   # add the amino acid that is 1 SNP away


def insert_codons(seq1, seq2, pad_left, pad_right, homology_length, triplets_to_insert):
    seqlen = len(seq1)
    for codon_position in range(0, len(seq1), 3):
        # note that seq2 is the wild type
        replaced_codon = seq2[codon_position : codon_position + 3]
        replaced_AA = translate(replaced_codon)
        left_padding_length = max([homology_length - codon_position, 0])
        if left_padding_length > 0:
            left_padding_seq = pad_left[-left_padding_length:].lower()
        else:
            left_padding_seq = ''
        
        right_padding_length = max([homology_length - (len(seq1) - codon_position), 0])
        if right_padding_length > 0:
            right_padding_seq = pad_right[:right_padding_length].lower()
        else:
            right_padding_seq = ''

        for i in triplets_to_insert:
            inserted_AA = translate(i)
            if inserted_AA in singleSNPcodons[replaced_codon]:
                singleSNPmutation = True
            else:
                singleSNPmutation = False
            left_seq  = seq1[max(0, codon_position - homology_length) : codon_position]
            right_seq = seq2[codon_position + 3 : min(seqlen, codon_position + homology_length + 3) ]
            yield(' '.join([replaced_codon, "("+replaced_AA+")", "to", inserted_AA, str(singleSNPmutation), left_padding_seq + left_seq, i, right_seq + right_padding_seq]))

mutagenic_repair_templates = insert_codons(fasta_1, fasta_2, upstream, downstream, 50, inserted_codons)

with open('test.out', 'w') as outfile:
    for i in mutagenic_repair_templates:
        outfile.write(i + '\n')


exit


# Take seq 1
# Take seq 2

# confirm same length

# replace i'th codon
# for i in range(0,(seq_length / 3)):
#seq_1[0:(3*i)] + <all codons> + seq_2[(3*i + 3):]

#30 0:30