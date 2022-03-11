#!/usr/bin/env python

####################
# DEFINE FUNCTIONS #
####################

def read_text(filename):
    with open(filename, 'r') as infile:
        rawText = infile.read()
    if debug:
        print("Read in text:\n" + rawText)
    return rawText

def unfold_alignment(rawText):
    text = rawText.splitlines()[3:]
    
    lineCount = 0
    seq1 = ''
    seq2 = ''
    alnScore = ''
    
    for line in text:
        if lineCount %4 == 0:
            seq1add = line.split()[1]
            line_len = len(seq1add)
            seq1 += seq1add
        elif lineCount %4 == 1:
            seq2add = line.split()[1]
            line_len = len(seq2add)
            seq2 += seq2add
        elif lineCount %4 == 2:
            alnAdd = line.lstrip()
            alnPadding = line_len - len(alnAdd)
            alnScore += ' ' * alnPadding + alnAdd
        elif lineCount %4 == 3:
            pass
        lineCount += 1
    # modify score string
    new_score = ''
    for aa1, aa2, score in zip(seq1, seq2,alnScore):
        if aa1 == "-" and aa2 != "-":
            new_score += 's'
        elif aa1 != "-" and aa2 == "-":
            new_score += 'f'
        else:
            new_score += score
    return [seq1, seq2, new_score]

def get_homology_regions(alignment_scores, length, specificity):
    if specificity == "high":
        pattern = r"[*]"
    elif specificity == "medium":
        pattern = r"[*;]"
    elif specificity == "low":
        pattern = r"[*;.]"
    pattern += "{%s,}" % length
    boundaries = [(m.start(0), m.end(0)) for m in re.finditer(pattern, alignment_scores)]
    return boundaries

def get_non_homology_regions(homology_regions, alignment_length):
    if len(homology_regions) == 0:
        return [(0, alignment_length)]
    else:
        return ( [(0,homology_regions[0][0])] +
                [ (homology_regions[x][1], homology_regions[x+1][0]) for x in range(0,len(homology_regions)-1)] +
                [(homology_regions[-1][1], alignment_length)] )

def get_overall_score(score):
    gaps = re.findall('[fs]+', score)
    nonperfect_align = re.findall(r'[^\*fs]', score)
    perfect_align = re.findall(r'\*', score)
    total_gap_length = sum([len(x) for x in gaps])
    total_align_length = len(nonperfect_align) + len(perfect_align)
    overall_score = float(len(perfect_align)) / (len(perfect_align) + len(nonperfect_align) + len(gaps))
    return [overall_score, total_align_length, total_gap_length, len(perfect_align), len(nonperfect_align), len(gaps)]

def get_index_offset(seq):
    seq_offset = []
    index_value = -1
    for i in seq:
        if i != "-":
            index_value += 1
            seq_offset.append(index_value)
        else:
            seq_offset.append("NA")
    return seq_offset

def return_combos(dict1, dict2, ranges, max_indel):
    for i in ranges:
        for start in range(i[0], i[1]):
            for end in range(start, i[1]):
                if abs(end - start) <= max_indel or max_indel == -1:
                    value1 = dict1[start]
                    value2 = dict2[end]
                    if not (value1 == 'NA' or value2 == 'NA'):
                        yield(value1, value2)
                else:
                    continue

def permute_alignment(alignment, n):
    alignment_list = list(alignment)
    for _ in range(n):
        random.shuffle(alignment_list)
        yield(''.join(alignment_list))

def wrap_fasta(seq, line_length=60):
    return '\n'.join([seq[i:i+line_length] for i in range(0, len(seq), line_length)])

def permute_peptide_fasta(fasta_text):
    split_text = fasta_text.rstrip().split('\n')
    header = split_text[0]
    seq = list(''.join(split_text[1:]).rstrip('*'))
    random.shuffle(seq)
    permuted_seq = header + '\n' + wrap_fasta(''.join(seq) + '*') + '\n'
    if debug == True:
        print(permuted_seq)
    return permuted_seq

def run_clustal(b_alignment):
    output = subprocess.run(["clustalo", "-i", "-", "--outfmt", "clu"], input=b_alignment, stdout=PIPE).stdout.decode()
    if debug:
        print(output)
    return output

def run_alignment(n_permutations, pep1, pep2, length, specificity):
    if n_permutations == 0:
        yield clustal_alignment(0, pep1, pep2, length, specificity)
    elif n_permutations > 0:
        for _ in range(n_permutations):
            yield clustal_alignment(random.randint(1,2), pep1, pep2, length, specificity)

def generate_repair_templates(  all_recombination_points, 
                                gene1, 
                                gene2,
                                pep1, 
                                pep2, 
                                dna1, 
                                dna2, 
                                upstream, 
                                downstream, 
                                homology_length, 
                                five_prime_padding, 
                                three_prime_padding,
                                primer_length,
                                oligo_length):
    for recombination_point in all_recombination_points:
        if not args.deletion or (args.deletion == True and recombination_point[1] >= recombination_point[0]):
            yield(repair_template(recombination_point, gene1, gene2, pep1, pep2, dna1, dna2, upstream, downstream, homology_length, five_prime_padding, three_prime_padding, primer_length, oligo_length))

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

def printRTs(RTs, uniqueType):
    global unique_chimeras
    for template in RTs:
        if uniqueType == "all":
            print(template.rt_formatted)
        else:
            if uniqueType == "protein":
                chimera = template.pepChimera
            elif uniqueType == "dna":
                chimera = template.dnaChimera
            if chimera not in unique_chimeras:
                unique_chimeras.append(chimera)
                #print(template.pepChimera)
                print(template.rt_formatted)

def warning_on_one_line(message, category, filename, lineno, file=None, line=None):
        return '%s:%s: %s:%s\n' % (filename, lineno, category.__name__, message)

##################
# DEFINE CLASSES #
##################

class clustal_alignment:
    help = 'clustal alignment of permuted protein sequence'
    # zero for unpermuted, 1 for permute seq 1, 2 for permute seq 2
    def __init__(self, permute, pep1, pep2, length, specificity): # 0 for unpermuted, 1 for permute seq 1
        if permute == 0:
            self.pep1 = pep1
            self.pep2 = pep2
        elif permute == 1:
            self.pep1 = pep1
            self.pep2 = permute_peptide_fasta(pep2)
        elif permute == 2:
            self.pep1 = permute_peptide_fasta(pep1)
            self.pep2 = pep2
        self.clustal = run_clustal(
                                    (self.pep1 + self.pep2).encode() 
                                    )
        self.aln1, self.aln2, self.scores = unfold_alignment(self.clustal)  # unfolded clustal alignment
        self.homology_regions = get_homology_regions(self.scores, length, specificity)
        self.non_homology_regions = get_non_homology_regions(self.homology_regions, len(self.aln1))
        self.offset1 = get_index_offset(self.aln1)
        self.offset2 = get_index_offset(self.aln2)
        self.non_homology_combos = list(return_combos(self.offset1, self.offset2, self.non_homology_regions, args.max_indel))
        self.homology_combos = list(return_combos(self.offset1, self.offset2, self.homology_regions, args.max_indel))
        self.aligned_combos = [(self.offset1[i], self.offset2[i]) for i in range(len(self.scores)) if self.scores[i] not in ["f","s"]]

class repair_template:
    help = 'repair template for generating an intended chimera'
    def __init__(self, recombination_point, gene1, gene2, pep1, pep2, dna1, dna2, upstream, downstream, homology_length, five_prime_padding, three_prime_padding, primer_length, oligo_length):
        self.idx1 = recombination_point[0] * 3
        self.idx2 = recombination_point[1] * 3
        self.dna1 = dna1
        self.dna2 = dna2
        self.pep1 = pep1
        self.pep2 = pep2
        self.gene1 = gene1.split("/")[-1]
        self.gene2 = gene2.split("/")[-1]
        self.homology_length = int(homology_length)
        self.dnaSeq1 = dna1[0:self.idx1]
        self.dnaSeq2 = dna2[self.idx2:]
        self.dnaChimera = self.dnaSeq1 + self.dnaSeq2
        self.pepChimera = translate(self.dnaChimera)
        self.homologyLeft = self.dnaSeq1[-(self.homology_length):]
        self.homologyRight =  self.dnaSeq2[:(self.homology_length)]
        # # .padLeft and .padRight indicate # of nucleotides missing from left or right repair template
        self.padLeft = self.homology_length - len(self.homologyLeft)
        self.padRight = self.homology_length - len(self.homologyRight)
        if self.padLeft > 0:
            self.homologyLeft = upstream[-(self.padLeft):] + self.homologyLeft
        if self.padRight > 0:
            self.homologyRight = self.homologyRight + downstream[:(self.padRight)]
        if five_prime_padding == '' and three_prime_padding == '':
            self.rt = self.homologyLeft + self.homologyRight
        else:
            self.paddingTotal = int(oligo_length) - ( 2*int(primer_length) + 2*int(homology_length))
            self.paddingLength = int(self.paddingTotal /2)
            if self.paddingLength == 0:
                self.paddingLeft = ''
                self.paddingRight = ''
            else:
                self.paddingLeft = five_prime_padding[:self.paddingLength]
                self.paddingRight = three_prime_padding[-self.paddingLength:]
            self.rt = self.paddingLeft + self.homologyLeft + self.homologyRight + self.paddingRight
        self.rt_formatted = ">%s:Start-%s|%s:%s-End|RT:%s_nt_each\n%s" % (self.gene1, self.idx1, self.gene2, self.idx2+1, self.homology_length, wrap_fasta(self.rt))
# oligo length = 2 * primer length + 2 * homology length + 2 * padding
class fasta:
    help = 'stores type, header, and sequence information for FASTA files'
    def __init__(self, filename):
        with open(filename, 'r') as infile:
            text = infile.readlines()
        self.header = text[0].replace(">","").strip()
        self.seq = ''.join(x.rstrip() for x in text[1:])


#######################################################################################################################
#                                                        MAIN
#######################################################################################################################

if __name__ == "__main__":

    ##################
    # IMPORT MODULES #
    ##################

    import os
    import sys
    import re
    import random
    import argparse
    import subprocess
    import warnings
    warnings.formatwarning = warning_on_one_line

    from subprocess import PIPE

    parser = argparse.ArgumentParser()
    parser.add_argument('gene1_filestem', type=str,
                        help='File stem for gene 1. Ensure proper naming of <GENE1_FILESTEM>.fasta')
    parser.add_argument('gene2_filestem', type=str,
                        help='File stem for gene 2. Ensure proper naming of <GENE2_FILESTEM>.fasta')
    parser.add_argument("--debug", 
                        help="Assist with debugging by increasing output verbosity",
                        action="store_true")
    parser.add_argument("--deletion", 
                        help='''Only generate RTs that are a deletion (second sequence start point is >= the
                        end point of the first''',
                        default=False,
                        action="store_true")
    parser.add_argument("--max-indel",
                        type=int,
                        nargs='?',
                        const=1,
                        default=0,
                        help='''Integer. When generating chimeric protein sequences, INDEL is 
                                the maximum number of codons added or removed by
                                transitioning between two points offset from their
                                codon partner in the alignment. See <figure> for 
                                explanation. To include all possible transitions, set
                                to -1. Default: 0 (only transition between codons at
                                their partner in the alignment, with no offset).''')
    parser.add_argument("--threshold-length",
                        type=int,
                        nargs='?',
                        const=1,
                        default=4,
                        help='''Integer. Minimum number of consecutive amino acids above the
                                specificity threshold to determine a region of
                                confident homology (defined by -S). Only used when
                                manually testing alignments--routine use does not
                                require this argument. Default: 4''')
    parser.add_argument("--primer-length",
                        type=int,
                        nargs='?',
                        const=1,
                        default=15,
                        help='''Integer. Defines length of primer seq that will be
                                added at later stage (on boths sides of repair template).
                                Total oligo array oligonucleotide length will be equal to
                               <Repair Template length> + 2*<primer length> + <padding> to
                                equal specified <Total Oligo Length>''')
    parser.add_argument("--oligo-length",
                        type=int,
                        nargs='?',
                        const=1,
                        default=210,
                        help='''Integer. Defines total length (max possible) from oligo library
                                synthesis.''')
    parser.add_argument("--five-prime-padding",
                        type=str,
                        nargs='?',
                        const=1,
                        default='',
                        help='''Defines nucleotide sequence to pad 5' (upstream) space if RT length
                                is insufficient to reach total oligo length. Nucleotides will be 
                                inserted 5' (upstream) the left RT homology arm. Default: No padding.''')
    parser.add_argument("--three-prime-padding",
                        type=str,
                        nargs='?',
                        const=1,
                        default='',
                        help='''Defines nucleotide sequence to pad 3' (downstream) space if RT length
                                is insufficient to reach total oligo length. Nucleotides will be 
                                inserted 3' (downstream) the right RT homology arm. Default: No padding.''')
    parser.add_argument("--repair-template-length",
                        type=int,
                        nargs='?',
                        const=1,
                        default=120,
                        help='''Integer. Length of repair template. Half the length will
                                be allocated to the first sequence, and half 
                                to the second sequence. Will be rounded down
                                by one if an odd integer is supplied.  Default: 50''')
    parser.add_argument("--specificity",
                        type=str,
                        nargs='?',
                        const=1,
                        default="high",
                        help='''String, accepts 'low' 'medium' or 'high'.
                                Defines degree of amino acid similarity from clustal alignment
                                for defining regions of confident homology. Used with -L.
                                'low' includes [.:*], 'medium includes' [:*], 'high' includes [*].
                                Default: high.''')
    parser.add_argument("--mode",
                        type=str,
                        nargs='?',
                        const=1,
                        default="aligned",
                        help='''String, accepts 'aligned' 'strict' 'all' and 'extensive'.
                                aligned (default): only transition at aligned codons
                                strict: only transition at aligned codons with confident homology (see: --threshold-length and --specificity)
                                all: every possible transition, bounded by regions of homology if --max-indel is nonzero
                                extensive: every possible transition (N choose 2)''')
    parser.add_argument("--permutations",
                        type=int,
                        nargs='?',
                        const=1,
                        help = '''Interger. Number of alignment permutations to run. Default: 100.''')
    parser.add_argument("--unique", 
                        help='''String, accepts 'protein' or 'dna' or 'all'. Default: all sequences (no deduplication).
                                Specifying 'protein' will deduplicate oligo sequences, retaining one representative
                                for a given protein sequence. 'dna' will retain identical dna sequences,
                                with the possibility of duplicate resulting protein sequences (which may be of
                                interest for codon bias).''',
                        type=str,
                        nargs='?',
                        const=1,
                        default="all")
    parser.add_argument("--upstream", 
                        help='''String. Defines the file name for fasta sequence upstream of the homologous genes.''',
                        type=str,
                        nargs='?',
                        const=1,
                        default=None)
    parser.add_argument("--downstream", 
                        help='''String. Defines the file name for fasta sequence downstream of the homologous genes.''',
                        type=str,
                        nargs='?',
                        const=1,
                        default=None)
    args = parser.parse_args()
    if args.debug:
        debug = True
        print("debug mode on, increasing verbosity.")
    else:
        debug = False

    if debug == True:
        print("""running with length = %s and specificity = %s""" % (args.threshold_length, args.specificity))

#######################################################################################################################

    ######################
    # VALIDATE ARGUMENTS #
    ######################

    try: assert args.threshold_length > 0, "ERROR: length must be > 0"
    except AssertionError as error: sys.exit(error)

    try: assert args.specificity in ["low", "medium", "high"], "ERROR: -s, --specificity must be 'low', 'medium', or 'high'"
    except AssertionError as error: sys.exit(error)

    try: assert args.unique in ["protein", "dna", "all"], "ERROR: --unique must be 'prot' or 'dna'"
    except AssertionError as error: sys.exit(error)

    try: assert args.upstream != None, "ERROR: --upstream required.'"
    except AssertionError as error: sys.exit(error)

    try: assert args.downstream != None, "ERROR: --downstream required.'"
    except AssertionError as error: sys.exit(error)

    try: assert args.mode in ["aligned", "strict", "all", "extensive"], "ERROR: --mode must be 'aligned' 'strict' 'all' or 'extensive'"
    except AssertionError as error: sys.exit(error)

    try: assert not (args.mode in ['aligned', 'strict'] and args.max_indel != 0), "ERROR: nonzero --max-indel incompatible with strict or aligned mode. Run with --mode all"
    except AssertionError as error: sys.exit(error)



#######################################################################################################################

    ##################
    # VALIDATE INPUT #
    ##################

    missing_files = 0
    for filename in [   args.gene1_filestem + '.fasta',
                        args.gene2_filestem + '.fasta',
                        args.upstream,
                        args.downstream]:
        if not os.path.isfile(filename):
            sys.stderr.write("Error: File %s does not exist\n" % filename)
            missing_files += 1
    if missing_files > 0:
        sys.exit("Aborting due to missing files\n")
    if args.repair_template_length % 2 != 0:
        args.repair_template_length -= 1
    homology_length = args.repair_template_length / 2

#######################################################################################################################

    ############
    # RUN CODE #
    ############

    dna1 = fasta(args.gene1_filestem + '.fasta')
    dna2 = fasta(args.gene2_filestem + '.fasta')

    pep1 = lambda: None
    pep2 = lambda: None

    pep1.header = args.gene1_filestem
    pep1.seq = translate(dna1.seq)
    pep2.header = args.gene2_filestem
    pep2.seq = translate(dna2.seq)


    pep_txt1 = '>%s\n%s\n' % (args.gene1_filestem, pep1.seq)
    pep_txt2 = '>%s\n%s\n' % (args.gene2_filestem, pep2.seq)

    
    upstream = fasta(args.upstream)
    downstream = fasta(args.downstream)

    alignment = list(run_alignment(0, pep_txt1, pep_txt2, args.threshold_length, args.specificity))[0]

    combos = []
    unique_chimeras = []

    if args.mode == "strict":
        combos += alignment.homology_combos
    elif args.mode == "aligned":
        combos += alignment.aligned_combos
    elif args.mode == "extensive":
        combos += [ (i,j) for i in range(len(pep1.seq)) for j in range(len(pep2.seq)) ]
    else:
        combos += alignment.aligned_combos
        combos += alignment.non_homology_combos
        combos += alignment.homology_combos

    # convert to unique list
    combos = list(set(combos))


    RTs = generate_repair_templates(combos, args.gene1_filestem, args.gene2_filestem, pep1.seq, pep2.seq, dna1.seq, dna2.seq, upstream.seq, downstream.seq, homology_length, args.five_prime_padding, args.three_prime_padding, args.primer_length, args.oligo_length)

    printRTs(RTs, args.unique)


    sys.exit()