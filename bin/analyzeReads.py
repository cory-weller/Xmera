#!/usr/bin/env python

import itertools
import regex as re
import gzip
import bisect

infile = "example.fastq"

known_upstream =  'ACTACCGGCTGATATCATCG'
#post_seq = 'GATCCCTGAGTAACCGGTTC'
known_downstream = 'CCTGAGTAACCGGTTC'
fullTolerance = 2
knownUpstreamTolerance = 1
knownDownstreamTolerance = 1

expectedLength = 150

# find pre_seq, error <= 1, get indices (pre_start, pre_end)
# find post_seq, error <= 1, get indices (post_start, post_end)
# if either match multiply, ignore and give warning
# Extract variable region as pre_end : post_start
# take first N reads as barcodeUpstream
# take last N reads as barcodeDownstream
#


class fastq_read:
    def __init__(self, lines):
        self.id, self.seq, self.line3, self.quality = lines

class read:
    def __init__(self, name, rawSeq):
        self.name = name
        self.rawSeq = rawSeq
        self.upBarcode = ''
        self.downBarcode = ''
        self.repairTemplate = ''
        self.count = 1
    #@property
    #def increment(self):
    #    self.count += 1
    def extractSeq(self, fullPattern, upstreamPattern, downstreamPattern):
        self.matches = re.search(fullPattern, self.rawSeq)
        if self.matches == None:
            print("no matches")
            pass
        elif len(self.matches.captures()) > 1:
            print("more than one pattern match, do what now?")
            pass
        else:
            self.extractedSeq = self.matches.captures()[0]
            #self.extractedStart, self.extractedStop = self.extractedSeq.span()
            self.preMatches = re.search(upstreamPattern, self.extractedSeq)
            self.preMatchString = self.preMatches.captures()[0]
            self.preMatchEnd = self.preMatches.ends()[0]
            self.postMatches = re.search(downstreamPattern, self.extractedSeq)
            self.postMatchString = self.postMatches.captures()[-1]
            self.postMatchStart = self.postMatches.starts()[-1]
            self.RTseq = self.extractedSeq[self.preMatchEnd : self.postMatchStart]

# a = read("testRead", 'TAANAGCNAAGCACCTTTCGAGAGGACGATGCCCGTGTCTAAATGATTCGACCAGCCTAAGAATGTTCAACGGCCCACTACCGGCTGATATCATCACTTTCTATTTGGGTTCGGAGGCAGGAGGGTCTAGTTTGGAGGTTGAGGACCTCAGCCTGGAGAACGCAATAGGTAAGCCATTCGATAAAACATTAGACACTACTTAGACCATTTGTGGCACCAGCTTGGAGCCAAGACTAAATCCTTAGTTCAGGATTTGAAGAATTACGAACTTTGCTGCAGTATCTCTCGGACTAGCGGCCGTCAACCTGTCTCCAAAGCCTGAGTAACCGGTTCGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAGCGCGCGTAATACGACTCACTATAGGGCGAATTGGGTACCGGCCGCAAATTAAAGCCTTCGAGCGTCCCAAAACCTTCTCAAGCAAGGTTTTCAGTATAATGTTACATGCGTACACGCGTCTGTACAGAAAAAAAAGAAAAATTTGAAATATAAATAACGTTCTTAATACTAACATAACTATAAAAAAATAAATAGGGACCTAGACTTCAGGTTGTCTAACTCCTTCCTTTTCGGTTAGAGCGGATGTGGGGGGAGGGCGTGAACGTAAGCGTGACATAACTAATTACATGACTCGAAAACATAAAAAACAAAAAAGCACCACCGACTCGGTGCCACTTTTCAAGTTGATAACGGACTACCCAAANNA')

# re.findall(pattern, seqs[0][1])

# a.extractSeq(fullPattern, upstreamPattern, downstreamPattern)


class reads:
    def __init__(self, f):
        self.all = []
        if f.endswith('.gz'):
            self.gzip = True
            self.filename = f.rstrip(".gz")
        else:
            self.gzip = False
            self.filename = f
        self.extension = self.filename.split('.')[-1].lower()
        if self.extension in ["fastq", "fq"]:
            self.readType = "fastq"
        elif self.extension in ["fasta", "fa"]:
            self.readType = "fasta"
        self.entries = {}
    def importGzipFastq(self):
        print("Importing gzipped fastq sequences from %s" % self.filename)
        with gzip.open(self.filename, 'rb') as f:
            for lines in itertools.zip_longest(*[f]*4):      
                rawSeq = lines[1].strip()
                if rawSeq not in self.all:
                    self.all.append(read(readName, rawSeq))
    def importFastq(self):
        print("Importing fastq sequences from %s" % self.filename)
        with open(self.filename, 'r') as f:
            for lines in itertools.zip_longest(*[f]*4):      
                rawSeq = lines[1].strip()
                if rawSeq not in self.all:
                    readName = lines[0].strip()
                    self.all.append(read(readName, rawSeq))
    def importGzipFasta(self):
        print("Importing gzipped fasta sequences from %s" % self.filename)
        with gzip.open(self.filename, 'rb') as f:
            samples = f.read().split(">")[1:]
            samples = [x.strip().split("\n") for x in samples]
            samples = [(x[0], ''.join(x[1:])) for x in samples]
    def importFasta(self):
        print("Importing fasta sequences from %s" % self.filename)
        with open(self.filename, 'r') as f:
            samples = f.read().split(">")[1:]
            samples = [x.strip().split("\n") for x in samples]
            samples = [(x[0], ''.join(x[1:])) for x in samples]
            for i in samples:
                readName = i[0]
                rawSeq = i[1]
                self.all.append(read(readName, rawSeq))

a = reads("example.fastq")
if a.gzip == True:
    if a.extension == "fasta":
        a.importGzipFasta()
    elif a.extension == "fastq":
        a.importGzipFastq()
else:
    if a.extension == "fasta":
        a.importFasta()
    elif a.extension == "fastq":
        a.importFastq()

fullPattern = re.compile("""(%s[ACTGN]{%s,}%s){e<=%s}""" % (known_upstream, expectedLength, known_downstream, fullTolerance))
upstreamPattern = re.compile("""(%s){e<=%s}""" % (known_upstream, knownUpstreamTolerance))
downstreamPattern = re.compile("""(%s){e<=%s}""" % (known_downstream, knownDownstreamTolerance))

b = a.all[0]

b.extractSeq(fullPattern, upstreamPattern, downstreamPattern)

# iterate through reads and extract known flanking regions
for i in a.all:
    regex_pattern


print("wait")

exit()

# input = readsFile(inFileName)



# def importReads(filename):
#     if getFiletype(filename) in ["fastq", "fq"
#     with open(filename) as f:
#         for lines in itertools.zip_longest(*[f]*4):
#             yield read(lines[1])
#             #yield(fastq_read([x.strip() for x in lines]))

# def importFasta(filename):
#     # opens a .fasta file and returns a list of 2-tuples
#     # [('headerOne', 'ACTGACTGAC'), ('headerTwo', 'GTACCATGA')]
#     with open(filename) as f:
#         # split into separate entries
#         samples = f.read().split(">")[1:]
#         samples = [x.strip().split("\n") for x in samples]

# def importFastq(filename):


# def iterate_chimera_database(seq1, seq2):
#     for i in range(0, len(seq1), 3):
#         yield seq1[0:i] + seq2[i:]
#     for i in range(0, len(seq1), 3):
#         yield seq2[0:i] + seq1[i:]

# def get_distance(string1, string2):
#     len1 = len(string1)
#     len2 = len(string2)
#     assert len1 == len2, "string lengths do not match!"    
#     score = 0
#     for i in range(len1):
#         if string1[i] == string2[i]:
#             score += 1
#     return(score)

# def levenshtein(s1, s2):
#         if s1 == s2: return 0
#         elif len(s1) == 0: return len(s2)
#         elif len(s2) == 0: return len(s1)
#         v0 = [None] * (len(s2) + 1)
#         v1 = [None] * (len(s2) + 1)
#         for i in range(len(v0)):
#             v0[i] = i
#         for i in range(len(s1)):
#             v1[0] = i + 1
#             for j in range(len(s2)):
#                 cost = 0 if s1[i] == s2[j] else 1
#                 v1[j + 1] = min(v1[j] + 1, v0[j + 1] + 1, v0[j] + cost)
#             for j in range(len(v0)):
#                 v0[j] = v1[j]
#         return v1[len(s2)]

# def find_closest_chimera(read):
#     for i in iterate_chimera_database(firstseq, secondseq):
#         get_distance(read.seq, i)


# seqs = importFasta('example.fasta')

# pattern = re.compile("""%s{e<=2}.*%s{e<=%s}""" % (pre_seq, post_seq, tolerance))
# pattern = re.compile("""%s{e<=2}[ACTGN]{200}""" % (pre_seq))

# pattern = re.compile("""(%s[ACTGN]{190,}%s){e<5}""" % (pre_seq, post_seq))

# re.findall(pattern, seqs[0][1])

# GATCCCTGAGTAACCGGTTC

# post_seq = 'CCTGAGTAACCGGTTC'
# post_seq = 'CCCTGAGTAACCGGTTC'
# post_seq = 'TCCCTGAGTAACCGGTTC'
# post_seq = 'ATCCCTGAGTAACCGGTTC'
# post_seq = 'GATCCCTGAGTAACCGGTTC'

# # iterate through file and add reads to dictionary of all reads




