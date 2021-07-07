#!/usr/bin/env python

import regex as re
import copy
import sys

variant_filename = sys.argv[1]

class hgvs:
    help = "hgvs variant object"
    def __init__(self, hgvs_text, ref_seq):
        self.refSeq = list(ref_seq)
        if ":" in hgvs_text:
            self.reference = hgvs_text.split(":")[0]
            self.description = hgvs_text.split(":")[1]
        else:
            self.acc = "NA"
            self.description = hgvs_text
        self.seqType, self.posEdit = self.description.split(".")
        if ">" in self.posEdit:
            # NC_0001.1:g.123A>G
            self.varType = "substitution"
            # split on boundry of letters and characters
            self.pos, self.edit = re.split('(\d+)', self.posEdit)[1:]
            self.pos = [int(self.pos)]
            self.before, self.after = self.edit.split(">")
        elif "del" in self.posEdit:
            # NC_0001.1:g.123_127del
            self.varType = "deletion"
            i = [int(x) for x in re.split('del', self.posEdit)[0].split("_")[:2]]
            self.pos = [min(i), max(i)]
        elif "ins" in self.posEdit:
            self.varType = "insertion"
            i, self.insert = (re.split('ins', self.posEdit)[:2])
            i = [int(x) for x in i.split("_")]
            self.pos = [min(i), max(i)]
        elif "dup" in self.posEdit:
            self.varType = "duplication"
            i = re.split('dup', self.posEdit)[0]
            i = [int(x) for x in i.split("_")]
            self.pos = [min(i), max(i)]
        else:
            if "+" in hgvs_text or "-" in hgvs_text:
                return()
            else:
                raise RuntimeError("HGVS variant %s is not properly formatted or is not a simple substitution, deletion, insertion, or duplication")
    def generateVar(self):
        assert hasattr(self, 'varType'), "%s is a malformed description, has no variant type" % self.description
        self.varSeq = copy.deepcopy(self.refSeq)
        if self.varType == "substitution":
            # position 10 in hgvs is index 9 in python
            assert self.refSeq[self.pos[0] - 1] == self.before, "'before' nucleotide %s should be %s" % (self.before, self.refSeq[self.pos[0] - 1])
            self.varSeq[self.pos[0] - 1] = self.after
        elif self.varType == "deletion":
            del self.varSeq[self.pos[0] - 1 : self.pos[1]]
        elif self.varType == "insertion":
            i = self.varSeq[:self.pos[0]] + list(self.insert) + self.varSeq[self.pos[1]-1:]
            self.varSeq = i
        elif self.varType == "duplication":
            if len(self.pos) == 1:
                i = ''.join(self.refSeq[self.pos[0] - 1])
                self.varSeq[self.pos[0] - 1] = i*2
            elif len(self.pos) == 2:
                self.insert = self.refSeq[self.pos[0] - 1 : self.pos[1]]
                i = self.varSeq[:self.pos[0]] + list(self.insert) + self.varSeq[self.pos[1]-1:]
                self.varSeq = i
        return(''.join(self.varSeq))
    def generateRT(self, shuffledSeq, armLength):
        self.shuffledSeq = list(shuffledSeq)
        self.armLength = armLength
        self.offset = (self.pos[0] - 1)%3 + 1
        self.boundary = self.pos[0] - self.offset
        leftArmEnd = self.pos[0] - self.offset
        self.leftRTarm = self.shuffledSeq[leftArmEnd - self.armLength: leftArmEnd]
        if self.varType == "substitution" or self.varType == "deletion":
            self.rightRTarm = self.varSeq[leftArmEnd : leftArmEnd + self.armLength + 3]
        elif self.varType == "insertion" or self.varType == "duplication":
            self.rightRTarm = self.varSeq[leftArmEnd : leftArmEnd + self.armLength +  3 + len(self.insert)]
        self.RT = ''.join(self.leftRTarm + self.rightRTarm)
        return(self.RT)

def importSingleSeqFasta(filename):
    with open(filename, 'r') as infile:
        seq = infile.readlines()[1:]
        seq = ''.join([x.strip() for x in seq])
        return(seq)

ercc4 = importSingleSeqFasta('ercc4.wt.fasta')

ercc4prime  = importSingleSeqFasta('ercc4.shuffle.fasta').lower()

#a = hgvs("NC_0001.1:g.123T>G", 100*'AG')
#b = hgvs("NC_0001.1:g.27_52del", 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')
#c = hgvs("NC_0001.1:g.26_27insTEEHEE", 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')
#d = hgvs("NC_0001.1:g.123_124dupTT", 'aaaaaa')
#e = hgvs("NC_0001.1:g.26_27dupTEEHEE", 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789')
#f = hgvs("c.1A>G", ercc4)

with open(variant_filename, 'r') as infile:
    for line in infile:
        variant_text = line.strip()
        if variant_text == '':
            continue
        i = hgvs(variant_text, ercc4)
        i.generateVar()
        i.generateRT(ercc4prime, 50)
        print(">ercc4:" + i.description + "\n" +  i.RT)
