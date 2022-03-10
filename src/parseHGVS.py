#!/usr/bin/env python

import regex as re
import sys
import math

geneName = sys.argv[1]
firstHomolog = sys.argv[2]
secondHomolog = sys.argv[3]
hgvsIdentifier = sys.argv[4]
codonsFilename = sys.argv[5]
totalLength = int(sys.argv[6])
desiredArmLength = int(totalLength / 2)

firstHomologFilename = "%s/%s-%s.fasta" % (geneName, geneName, firstHomolog)
secondHomologFilename = "%s/%s-%s.fasta" % (geneName, geneName, secondHomolog)
downstreamFilename = "%s/%s-upstream.fasta" % (geneName, geneName)
upstreamFilename = "%s/%s-downstream.fasta" % (geneName, geneName)
hgvsFilename = "%s/%s-%s.txt" % (geneName, geneName, hgvsIdentifier)

class hgvs:
    help = "hgvs variant object"
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength):
        global codons
        self.upstreamSeq = upstreamSeq
        self.downstreamSeq = downstreamSeq
        self.desiredArmLength = desiredArmLength
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
    def buildRepairTemplate(self, leftArm, middle, rightArm, upstreamSeq, downstreamSeq, basesInserted, desiredArmLength):
        self.armLength = desiredArmLength
        self.leftHA = leftArm[-self.armLength : ]
        self.nLeft = self.armLength - len(self.leftHA)
        self.leftHA = upstreamSeq[(len(upstreamSeq) - self.nLeft):] + self.leftHA
        self.middle = middle
        self.rightHA = rightArm[:self.armLength]
        self.nRight = self.armLength - len(self.rightHA)
        self.rightHA = self.rightHA + downstreamSeq[:self.nRight]
        self.oligo = self.leftHA + self.middle.upper() + self.rightHA
        # Ensure oligo isn't too long
        self.extraBases = len(self.oligo) - (2*self.armLength)
        if self.extraBases > 0:
            self.trimEach = math.trunc(self.extraBases/2)
            # remove N/2 from start and end:
            self.oligo = self.oligo[self.trimEach : (-1*self.trimEach)]
            if self.extraBases % 2 != 0:
                self.oligo = self.oligo[1:]
        return self.oligo


class insertion(hgvs):
    # e.g.  g.123_124insACTG
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength):
        super().__init__(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        i, self.insert = (re.split('ins', self.posEdit)[:2])
        self.basesInserted = len(self.insert)
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
                                        self.basesInserted,
                                        self.desiredArmLength)
        return a


class duplication(hgvs):
    # treat like an insertion, where insertion is the duplicated region
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength):
        super().__init__(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        i = [int(x) for x in re.split('dup', self.posEdit)[0].split("_")]
        self.dupStart = min(i) - 1
        self.dupEnd = max(i) - 1
        self.basesInserted = self.dupEnd - self.dupStart + 1
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
                                        self.basesInserted, 
                                        self.desiredArmLength)
        return a

 
class deletion(hgvs):
    # e.g.  g.123_128del
    # e.g.  g.123del
    def __init__(self, hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength):
        super().__init__(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        i = [int(x) for x in re.split('del', self.posEdit)[0].split("_")[:2]]
        self.delStart = min(i) - 1
        self.delEnd = max(i) - 1
        self.basesInserted = -1 * (self.delEnd - self.delStart + 1)
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
                                        self.basesInserted, 
                                        self.desiredArmLength)
        return a

class substitution(hgvs):
    pass


def generateHGVS(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength):
        if ":" in hgvs_text:
            description = hgvs_text.split(":")[1]
        else:
            description = hgvs_text
        posEdit = description.split(".")[1]
        if ">" in posEdit:
            return substitution(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        elif "del" in posEdit:
            return deletion(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        elif "ins" in posEdit:
            return insertion(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        elif "dup" in posEdit:
            return duplication(hgvs_text, seq1, seq2, upstreamSeq, downstreamSeq, desiredArmLength)
        else:
            if "+" in hgvs_text or "-" in hgvs_text:
                return()
            else:
                raise RuntimeError("HGVS variant %s is not properly formatted or is not a simple substitution, deletion, insertion, or duplication")


# get most abundant codon per AA:
def getAbundantCodons():
    with open(codonsFilename, "r") as infile:
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


upstream = importSingleSeqFasta(upstreamFilename).lower()
firstHomolog = importSingleSeqFasta(firstHomologFilename).lower()
secondHomolog = importSingleSeqFasta(secondHomologFilename).upper()
downstream = importSingleSeqFasta(downstreamFilename).lower()

with open(hgvsFilename, 'r') as infile:
    for line in infile:
        variant_text = line.strip()
        if variant_text == '':
            continue
        i =  generateHGVS(variant_text, firstHomolog, secondHomolog, upstream, downstream, desiredArmLength)
        RT = i.getRT()
        print(">%s:%s\n%s" % (geneName, i.description, RT))

# sys.exit("Exiting successfully")