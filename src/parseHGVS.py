#!/usr/bin/env python

import regex as re
import sys

variant_filename = sys.argv[1]

class hgvs:
    help = "hgvs variant object"
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength):
        global codons
        self.upstreamSeq = upstreamSeq
        self.downstreamSeq = downstreamSeq
        self.armLength = armLength
        self.refSeq = list(seq1)
        if ":" in hgvs_text:
            self.reference = hgvs_text.split(":")[0]
            self.description = hgvs_text.split(":")[1]
        else:
            self.acc = "NA"
            self.description = hgvs_text
        self.seqType, self.posEdit = self.description.split(".")
    def roundUp(self, n):
        return(n + 3 - n%3)
    def roundDown(self, n):
        return(n - n%3)
    def translate(self, dna_seq): 
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
    def buildRepairTemplate(self, leftArm, middle, rightArm, upstreamSeq, downstreamSeq, armLength):
        self.leftHA = leftArm[-armLength : ]
        self.nLeft = armLength - len(self.leftHA)
        self.leftHA = upstreamSeq[(len(upstreamSeq) - self.nLeft):] + self.leftHA
        self.middle = middle
        self.rightHA = rightArm[:armLength]
        self.nRight = armLength - len(self.rightHA)
        self.rightHA = self.rightHA + downstreamSeq[:self.nRight]
        self.oligo = self.leftHA + self.middle.upper() + self.rightHA
        return self.oligo


class insertion(hgvs):
    # e.g.  g.123_124insACTG
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength):
        super().__init__(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        i, self.insert = (re.split('ins', self.posEdit)[:2])
        i = [int(x) for x in i.split("_")]
        self.L = min(i) - 1
        self.R = max(i) - 1
        self.afterIns = seq2[self.R : ]
        assert self.L + 1 == self.R, "insertion indices %s and %s are not adjacent" % (self.L + 1, self.R + 1)
        self.insertionLength = len(self.insert)
        self.leftArmLength = self.L - (self.L % 3)
        self.leftArm = seq1[:self.leftArmLength]
        self.afterLeftArm = seq2[self.leftArmLength : ]
        self.nNeededL = (1 + self.L - self.leftArmLength) % 3
        self.nNeededR = -(self.nNeededL + len(self.insert)) % 3
        self.middlePre = self.afterLeftArm[:self.nNeededL] + self.insert + self.afterIns[:self.nNeededR]
        assert len(self.middlePre) % 3 == 0, "middle fragment %s not multiple of 3 %s" % (self.middlePre)
        self.middle = ''.join([codons[x] for x in self.translate(self.middlePre)])
        assert len(self.leftArm) % 3 == 0, "left fragment not a multiple of 3"
        self.rightArm = self.afterIns[self.nNeededR:]
    def getRT(self):
        a = self.buildRepairTemplate(   self.leftArm, 
                                        self.middle, 
                                        self.rightArm, 
                                        self.upstreamSeq, 
                                        self.downstreamSeq, 
                                        self.armLength)
        return a


class duplication(hgvs):
    # treat like an insertion, where insertion is the duplicated region
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength):
        super().__init__(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        i = [int(x) for x in re.split('dup', self.posEdit)[0].split("_")]
        self.dupStart = min(i) - 1
        self.dupEnd = max(i) - 1
        self.insert = seq2[self.dupStart : self.dupEnd + 1]
        self.L = self.dupEnd 
        self.R = self.dupEnd + 1
        self.afterIns = seq2[self.R : ]
        assert self.L + 1 == self.R, "insertion indices %s and %s are not adjacent" % (self.L + 1, self.R + 1)
        self.insertionLength = len(self.insert)
        self.leftArmLength = self.L - (self.L % 3)
        self.leftArm = seq1[:self.leftArmLength]
        self.afterLeftArm = seq2[self.leftArmLength : ]
        self.nNeededL = (1 + self.L - self.leftArmLength) % 3
        self.nNeededR = -(self.nNeededL + len(self.insert)) % 3
        self.middlePre = self.afterLeftArm[:self.nNeededL] + self.insert + self.afterIns[:self.nNeededR]
        assert len(self.middlePre) % 3 == 0, "middle fragment %s not multiple of 3 %s" % (self.middlePre)
        self.middle = ''.join([codons[x] for x in self.translate(self.middlePre)])
        assert len(self.leftArm) % 3 == 0, "left fragment not a multiple of 3"
        self.rightArm = self.afterIns[self.nNeededR:]
    def getRT(self):
        a = self.buildRepairTemplate(   self.leftArm, 
                                        self.middle, 
                                        self.rightArm, 
                                        self.upstreamSeq, 
                                        self.downstreamSeq, 
                                        self.armLength)
        return a

 
class deletion(hgvs):
    # e.g.  g.123_128del
    # e.g.  g.123del
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength):
        super().__init__(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        i = [int(x) for x in re.split('del', self.posEdit)[0].split("_")[:2]]
        self.delStart = min(i) - 1
        self.delEnd = max(i) - 1
        self.delStartAdj = self.delStart - (self.delStart % 3)
        self.nNeeded = (3 - (self.delStart % 3)) % 3
        self.leftArm = seq1[: self.delStartAdj]
        self.rightArm = seq2[self.delEnd + 1 + self.nNeeded :]
        self.middlePre = seq2[self.delStartAdj : self.delStart] + seq2[self.delEnd + 1 : (self.delEnd + self.nNeeded + 1)]
        assert len(self.middlePre) % 3 == 0, "middle fragment not multiple of 3"
        self.middle = ''.join([codons[x] for x in self.translate(self.middlePre)])
    def getRT(self):
        a = self.buildRepairTemplate(   self.leftArm, 
                                        self.middle, 
                                        self.rightArm, 
                                        self.upstreamSeq, 
                                        self.downstreamSeq, 
                                        self.armLength)
        return a

class substitution(hgvs):
    pass


def generateHGVS(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength):
        if ":" in hgvs_text:
            description = hgvs_text.split(":")[1]
        else:
            description = hgvs_text
        posEdit = description.split(".")[1]
        if ">" in posEdit:
            return substitution(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        elif "del" in posEdit:
            return deletion(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        elif "ins" in posEdit:
            return insertion(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        elif "dup" in posEdit:
            return duplication(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, armLength)
        else:
            if "+" in hgvs_text or "-" in hgvs_text:
                return()
            else:
                raise RuntimeError("HGVS variant %s is not properly formatted or is not a simple substitution, deletion, insertion, or duplication")


# get most abundant codon per AA:
def getAbundantCodons():
    with open("seqs/codons.txt", "r") as infile:
        codons = {}
        for line in infile:
            triplet, AA, relativeFreq, _, _ = line.strip().split()
            triplet = triplet.replace("U", "T")
            if AA not in codons or relativeFreq > codons[AA][1]:
                codons[AA] = [triplet, relativeFreq]
        for i in codons:
            codons[i] = codons[i][0]
    return(codons)

codons = getAbundantCodons()


def importSingleSeqFasta(filename):
    with open(filename, 'r') as infile:
        seq = infile.readlines()[1:]
        seq = ''.join([x.strip() for x in seq])
        return(seq)


upstream = importSingleSeqFasta('seqs/ERCC4.upstream.fasta').lower()
ercc4prime  = importSingleSeqFasta('seqs/ERCC4.min_homology.fasta').lower()
ercc4 = importSingleSeqFasta('seqs/ERCC4.cds.fasta').upper()
downstream = importSingleSeqFasta('seqs/ERCC4.downstream.fasta').lower()
HAlength = 80

with open(variant_filename, 'r') as infile:
    for line in infile:
        variant_text = line.strip()
        if variant_text == '':
            continue
        i =  generateHGVS(variant_text, ercc4prime, ercc4, upstream, downstream, HAlength)
        RT = i.getRT()
        print(">ercc4:" + i.description + "\n" +  RT)

sys.exit("Exiting successfully")
