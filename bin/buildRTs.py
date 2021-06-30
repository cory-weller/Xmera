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
    output = subprocess.run(["./clustalo", "-i", "-", "--outfmt", "clu"], input=b_alignment, stdout=PIPE).stdout.decode()
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

# Unused
# def write_permutations(permutations, n_permutations, outfile_name="permutations.out"):
#     print("Writing %s permutations to %s" % (n_permutations, outfile_name))
#     with open(outfile_name, "w") as outfile:
#         current_i = 0
#         for i in permutations:
#             current_i += 1
#             if current_i % 10 == 0:
#                 print(current_i)
#             outfile.write('\n'.join(unfold_alignment(i)) + '\n\n')

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

class repair_template:
    help = 'repair template for generating an intended chimera'
    def __init__(self, recombination_point, gene1, gene2, pep1, pep2, dna1, dna2, upstream, downstream, homology_length, five_prime_padding, three_prime_padding, primer_length, oligo_length):
        self.idx1 = recombination_point[0] * 3
        self.idx2 = recombination_point[1] * 3
        self.dna1 = dna1
        self.dna2 = dna2
        self.pep1 = pep1
        self.pep2 = pep2
        self.gene1 = gene1
        self.gene2 = gene2
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
    from subprocess import PIPE

    parser = argparse.ArgumentParser()
    parser.add_argument('gene1_filestem', type=str,
                        help='File stem for gene 1. Ensure proper naming of <GENE1_FILESTEM>.pep.fasta and <GENE1_FILESTEM>.cdna.fasta')
    parser.add_argument('gene2_filestem', type=str,
                        help='File stem for gene 2. Ensure proper naming of <GENE2_FILESTEM>.pep.fasta and <GENE2_FILESTEM>.cdna.fasta')
    parser.add_argument("--debug", 
                        help="Assist with debugging by increasing output verbosity",
                        action="store_true")
    parser.add_argument("--all", 
                        help='''By default, this program will generate a list for only
                                potential chimeras where transitions between homologs
                                occurs at a point of confident homology. Setting the
                                --all flag will also generate chimeras with transitions
                                outside of regions of confident homology. Used with
                                -n, --indel. See <figure> for explanation.''',
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
#    parser.add_argument("--padding",
#                        type=int,
#                        nargs='?',
#                        const=1,
#                        default=1000,
#                        help='''Integer. Defines the number of nucleotides upstream and downstream
#                                of the genomic DNA for the gene of interest in the
#                                user-provided dna fasta file. Default value: 1000 bp,
#                                i.e. genomic dna +/- 1 kilobase.''')
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
                                 Only used when manually testing alignments--routine use 
                                does not require this argument. Default: high.''')
    parser.add_argument("--strict", 
                        help='''Specify to only include repair templates in regions
                               of high homology confidence''',
                        action="store_true")
    parser.add_argument("--permutations",
                        type=int,
                        nargs='?',
                        const=1,
                        help = '''Interger. Number of alignment permutations to run. Default: 100.''')
    parser.add_argument("--unique", 
                        help='''String, accepts 'protein' or 'dna'. 
                                Default 'prot' will deduplicate oligo sequences, retaining one representative
                                for a given protein sequence. 'dna' will retain identical dna sequences,
                                with the possibility of duplicate resulting protein sequences (which may be of
                                interest for codon bias).''',
                        type=str,
                        nargs='?',
                        const=1,
                        default="prot")
    parser.add_argument("--flanking", 
                        help='''String. Defines the file stem signifying flanking regions upstream and downstream the
                                homologous genes, i.e. <FILE_STEM>.upstream and <FILE_STEM>.downstream.''',
                        type=str,
                        nargs='?',
                        const=1,
                        default="flanking")
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

    try: assert args.unique in ["protein", "dna"], "ERROR: --unique must be 'prot' or 'dna'"
    except AssertionError as error: sys.exit(error)
    
#######################################################################################################################

    ##################
    # VALIDATE INPUT #
    ##################

    missing_files = 0
    for filename in [   args.gene1_filestem + '.fasta',
                        args.gene2_filestem + '.fasta',
                        args.flanking +  ".upstream.fasta", 
                        args.flanking + ".downstream.fasta"]:
        if not os.path.isfile(filename):
            print("Error: File %s does not exist" % filename)
            missing_files += 1
    if not os.path.isfile("clustalo"):
        print("Error: clustalo binary not found")
        missing_files += 1
    if missing_files > 0:
        sys.exit("Aborting due to missing files")
    if args.repair_template_length % 2 != 0:
        args.repair_template_length -= 1
    homology_length = args.repair_template_length / 2

#######################################################################################################################

    ############
    # RUN CODE #
    ############

    #pep1 = read_text(args.gene1_filestem + '.pep.fasta').rstrip() + "\n"    # Ensures file ends with newline. Required in this format for clustal
    #pep2 = read_text(args.gene2_filestem + '.pep.fasta').rstrip() + "\n"    # Ensures file ends with newline Required in this format for clustal



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

    



    upstream = fasta(args.flanking + ".upstream.fasta")
    downstream = fasta(args.flanking + ".downstream.fasta")

    if args.all:
        transition_points = [ [i,j] for i in range(len(pep1.seq)) for j in range(len(pep2.seq)) ]
        all_RTs = generate_repair_templates(transition_points, 
                                            args.gene1_filestem, 
                                            args.gene2_filestem, 
                                            pep1.seq, pep2.seq, 
                                            dna1.seq, dna2.seq, 
                                            upstream.seq, 
                                            downstream.seq, 
                                            homology_length, 
                                            args.five_prime_padding, 
                                            args.three_prime_padding, 
                                            args.primer_length, 
                                            args.oligo_length
                                            )
        for template in all_RTs:
            print(template.rt_formatted)
        sys.exit("exiting")

    alignment = list(run_alignment(0, pep_txt1, pep_txt2, args.threshold_length, args.specificity))[0]

    
    non_homology_RTs = generate_repair_templates(alignment.non_homology_combos, args.gene1_filestem, args.gene2_filestem, pep1.seq, pep2.seq, dna1.seq, dna2.seq, upstream.seq, downstream.seq, homology_length, args.five_prime_padding, args.three_prime_padding, args.primer_length, args.oligo_length)
    homology_RTs = generate_repair_templates(alignment.homology_combos, args.gene1_filestem, args.gene2_filestem, pep1.seq, pep2.seq, dna1.seq, dna2.seq, upstream.seq, downstream.seq, homology_length, args.five_prime_padding, args.three_prime_padding, args.primer_length, args.oligo_length)

    unique_chimeras = []

    if args.unique == "protein":
        for template in homology_RTs:
            if(template.pepChimera) not in unique_chimeras:
                unique_chimeras.append(template.pepChimera)
                print(template.rt_formatted)
        if not args.strict:
            for template in non_homology_RTs:
                if(template.pepChimera) not in unique_chimeras:
                    unique_chimeras.append(template.pepChimera)
                    print(template.rt_formatted)

    elif args.unique == "dna":
        for template in homology_RTs:
                print(template.rt_formatted)
        if not args.strict:
            for template in non_homology_RTs:
                print(template.rt_formatted)




  #  print(alignment.non_homology_combos)


    #print(alignment.non_homology_combos)

    sys.exit()

    if not args.permutations == None:
        print("Building permutation generator for %s permutations..." % (args.permutations))
        permutations = run_alignment(args.permutations, pep1, pep2, args.threshold_length, args.specificity)
        print("Built permutation generator!")
        for i in permutations:
            print(i.homology_regions)
            print(i.homology_combos)

    elif args.permutations == None:
        alignment = run_alignment(0, pep1, pep2, args.threshold_length, args.specificity)
        i = list(alignment)[0]
        if not os.path.isfile("true_alignment.out"):
            with open("true_alignment.out", 'w') as outfile:
                outfile.write(                           
                 '\t'.join(["specificity", "match_length", "homology_regions", "alignment_length", "n_chimeras"]) + '\n'
                        )
        with open("true_alignment.out", 'a') as outfile:
            outfile.write('\t'.join([str(x) for x in [args.specificity, args.threshold_length, len(i.homology_regions), len(i.aln1), len(i.homology_combos)]]) + '\n')

    sys.exit("Exited successfully!")

#######################################################################################################################
# NOTES / TODO

# Round 1 permutations to determine area to focus in on
# Round 2 permutations to get best FDR
# Option to just get indel = 0 chimeras
# Option to only make chimeras within high confidence homology regions
