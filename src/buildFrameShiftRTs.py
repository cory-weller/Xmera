#!/usr/bin/env python

## FUNCTIONS

def read_text(filename):
    with open(filename, 'r') as infile:
        rawText = infile.read()
    return rawText

def read_fasta(filename):
    with open(filename, 'r') as infile:
        seq = ''.join([x.strip() for x in infile.readlines()[1:]])
        return(seq)

def wrap_fasta(seq, line_length=80):
    return '\n'.join([seq[i:i+line_length] for i in range(0, len(seq), line_length)])




def is_file(arg, provided):
    if provided == None:
        sys.stderr.write("Error: %s file is not provided.\n" % (arg))
        return False
    elif not os.path.isfile(provided):
        sys.stderr.write("Error: %s '%s' does not exist or is misspelled.\n" % (arg, provided))
        return False
    else:
        return True

## CLASSES
class frameshift_repair_template:
    help = 'repair template for generating an intended mutant chimera withh a frame shift of -1, -2, +1 or +2'
    def __init__(self, seq1, seq2, upstream, downstream, totalLength, frameshift):
        self.seq1 = seq1
        self.seq2 = seq2
        self.upstream = upstream
        self.downstream = downstream
        self.totalLength = totalLength
        self.frameshift = frameshift
        self.linkerLength = self.frameshift + 3
        self.armLength = math.floor((self.totalLength - self.linkerLength) / 2)
    def buildRepairTemplate(self, codonPosition):
        self.codonPosition = codonPosition
        self.originalCodon = self.seq2[codonPosition : codonPosition + 3]
        self.leftArm = self.seq1[:codonPosition][-self.armLength:]
        self.rightArm = self.seq2[(codonPosition+3):][:self.armLength]
        self.leftShortBases = self.armLength - len(self.leftArm)
        self.rightShortBases = self.armLength - len(self.rightArm)
        if self.leftShortBases > 0:
            self.leftPad = self.upstream[-(self.leftShortBases):]
        else:
            self.leftPad = ''
        if self.rightShortBases > 0:
            self.rightPad = self.downstream[:self.rightShortBases]
        else:
            self.rightPad = ''
        if self.frameshift == -3:
            self.linker = ''
        if self.frameshift == -2:
            self.linker = self.originalCodon[2:]
        if self.frameshift == -1:
            self.linker = self.originalCodon[1:]
        if self.frameshift > 0:
            self.linker = self.frameshift*'C' + self.originalCodon
        self.repairTemplate = self.leftPad + self.leftArm + self.linker + self.rightArm + self.rightPad
        if self.frameshift > 0:
            self.sign = '+' 
        elif self.frameshift < 0:
            self.sign = ''
        self.header = ">%s:%s_shift_%s%s" % (args.prepend, self.codonPosition, self.sign, self.frameshift)
        if args.wrap == True:
            self.formatted = """%s\n%s""" % (self.header, wrap_fasta(self.repairTemplate))
        else:
            self.formatted = """%s\n%s\n%s\n%s\n%s\n%s""" % (self.header, self.leftPad, self.leftArm, self.linker, self.rightArm, self.rightPad)
        return(self.formatted)
    def __iter__(self):
        for codonPosition in range(0, len(self.seq1), 3):
            yield self.buildRepairTemplate(codonPosition)



# PARSE ARGUMENTS
if __name__ == "__main__":

    import argparse
    import sys
    import os


    parser = argparse.ArgumentParser()
    parser.add_argument('--first', 
        type=str,
        nargs='?',
        const=1,
        help="File name for first (i.e. upstream) sequence.")
    parser.add_argument('--second',
        type=str,
        nargs='?',
        const=1,
        help="File name for second (i.e. downstream) sequence.")
    parser.add_argument("--upstream",
        type=str,
        nargs='?',
        const=1,
        help='''fasta file containing upstream DNA sequence''')
    parser.add_argument("--length",
        type=int,
        nargs='?',
        const=1,
        default=189,
        help='''length of combined homology arms and inserted codon (default:189)''')
    parser.add_argument("--downstream",
        type=str,
        nargs='?',
        const=1,
        help='''fasta file containing downstream DNA sequence''')
    parser.add_argument("--prepend",
        type=str,
        nargs='?',
        const=1,
        default='',
        help='''string used to prepend fasta header names, if desired''')
    parser.add_argument("--wrap",
        action='store_true',
        help='''if this argument is used, repair templates will be concatenated and wrapped at 80 chars.
                Otherwise, the various elements of repair templates will be printed on separate lines:
                leftPadding, leftArm, frame-shift linker, rightArm, rightPadding.''')
    parser.add_argument("--frameshift",
        type=int,
        nargs='?',
        const=1,
        default=0,
        help='''[Integer]: size of frameshift to introduce. Allows -3, -2, -1, 1, 2''')
    args = parser.parse_args()

    missing_files = 0
    arg_list = ["first", "second", "upstream", "downstream"]
    provided_list = [args.first, args.second, args.upstream, args.downstream]
    for arg, provided in zip(arg_list, provided_list):
        if not is_file(arg, provided):
            missing_files += 1
    if missing_files > 0:
        sys.exit("Aborting due to missing files!\n")

        
    #python3.8 mut.py fasta1 fasta2 condontable 

    # IMPORT MODULES
    import math

    # RUN CODE
    fasta_1 = read_fasta(args.first).upper()
    fasta_2 = read_fasta(args.second).upper()

    # Assert sequences are identical length
    assert (len(fasta_1) == len(fasta_2)), "Sequences not same length!"


    # Load upstream and downstream seq
    upstream = read_fasta(args.upstream)
    downstream = read_fasta(args.downstream)


    mutagenic_repair_templates = frameshift_repair_template(fasta_1, fasta_2, upstream, downstream, args.length, args.frameshift)
    
    for RT in mutagenic_repair_templates:
        print(RT)


    exit