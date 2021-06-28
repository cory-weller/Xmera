# Generating codon-shuffled versions of the CAD2 allele

## Retrieving PCA1 and CAD2 coding sequences

PCA1 retrieved from https://www.yeastgenome.org/locus/S000000499

Renamed from `S288C_YBR295W_PCA1_coding.fsa` to `PCA1.cds.fasta`

CAD2 allele retrieved from https://www.ebi.ac.uk/ena/browser/view/AB027571

Renamed from `AB027571.1.fasta` to `CAD2.cds.fasta`:

```
mv S288C_YBR295W_PCA1_coding.fsa PCA1.cds.fasta
mv AB027571.1.fasta CAD2.cds.fasta
```

Codon table retrieved from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932&aa=1&style=N

Body of the above table's data saved as text file `codonTable.txt` 

The contents of `codonTable.txt` is parsed and processed by [`formatCodons.py`](formatCodons.py) into `codons.tab`:

```
( if [ ! -f "codons.tab" ]; then python3 formatCodons.py > "codons.tab"; fi )
```

## Generate aligned sequences
```
echo ">PW5_PCA1_pep" > PW5_PCA1.pep
python3 ./translate.py PW5_PCA1.cds.fasta >> PW5_PCA1.pep

echo ">BY_PCA1_pep" > BY_PCA1.pep
python3 ./translate.py BY_PCA1.cds.fasta >> BY_PCA1.pep

cat BY_PCA1.pep | fold -w 60 > align.in.pep
cat PW5_PCA1.pep | fold -w 60 >> align.in.pep

clustalo -i align.in.pep > align.pep

lineNo2=$(grep -n ">" align.pep | tail -n 1 | cut -d ":" -f 1)
let lineNo1="${lineNo2}-1"

# save separate files PW5_PCA1.align.pep and BY_PCA1.align.pep
head -n ${lineNo1} align.pep > BY_PCA1.align.pep
tail -n +${lineNo2} align.pep > PW5_PCA1.align.pep

echo ">PW5_PCA1" > PW5_PCA1.align.dna
python3 ./aligndna.py PW5_PCA1 | fold -w 60 >> PW5_PCA1.align.dna

echo ">BY_PCA1" > BY_PCA1.align.dna
python3 ./aligndna.py BY_PCA1 | fold -w 60 >> BY_PCA1.align.dna

Rscript shuffle.R BY_PCA1.align.dna PW5_PCA1.align.dna
```


## Generate shuffled sequences
```
Rscript shuffle.R BY_PCA1.align.dna PW5_PCA1.align.dna
```
The script [`shuffle.R`](shuffle.R) generates four new `fasta` files with various % identity shared with `PCA1.cds.fasta` (in addition to the original `CAD2.cds.fasta` which is 3651 nt long, has a Levenshtein distance of 52 for 98.58% identity with PCA1):
 * `CAD2.cds.min.fasta`
 * `CAD2.cds.low.fasta`
 * `CAD2.cds.med.fasta`
 * `CAD2.cds.high.fasta`

 ## Removing homopolymers

 Then, I used [`homopolymers.py`](homopolymers.py) to remove homopolymers. See `class fasta` within the script for exact replacements.

```
./homopolymers.py PW5_PCA1.align.min.fasta BY_PCA1.align.dna > PW5_PCA1.min_homology.fasta
./homopolymers.py PW5_PCA1.align.low.fasta BY_PCA1.align.dna > PW5_PCA1.low_homology.fasta
./homopolymers.py PW5_PCA1.align.medium.fasta BY_PCA1.align.dna > PW5_PCA1.medium_homology.fasta
./homopolymers.py PW5_PCA1.align.high.fasta BY_PCA1.align.dna > PW5_PCA1.high_homology.fasta
./homopolymers.py PW5_PCA1.align.max.fasta BY_PCA1.align.dna > PW5_PCA1.max_homology.fasta

for file in PW5_PCA1.*homology.fasta; do
    python3 ./formatfasta.py $file | fold -w 60 >> PW5_PCA1_gene_blocks.fasta
done

```

The following code block will show presence of 5+ repeat homopolymers in the original `CAD2.cds.fasta`:
```
for file in CAD2.cds.fasta; do
    for nucleotide in A C T G; do
        cat ${file} | tr -d "\n" | grep -E "[${nucleotide}]{5,}"
    done
done
```

Whereas when checking the newly generated `CAD2` files, no homopolymers are found:
```
for file in CAD2.high.fasta CAD2.medium.fasta CAD2.low.fasta CAD2.min.fasta; do
    for nucleotide in A C T G; do
        cat ${file} | tr -d "\n" | grep -E "[${nucleotide}]{5,}"
    done
done
```

## Final levels of identity for variable CAD2
| CAD2 version                             | Percent identity| Levenshtein distance |
|------------------------------------------|-----------------|----------------------|
| [`CAD2.cds.fasta`](CAD2.cds.fasta)       | 98.58%          |  52                  |
| [`CAD2.high.fasta`](CAD2.high.fasta)     | 86.88%          |  479                 |
| [`CAD2.medium.fasta`](CAD2.medium.fasta) | 74.23%          |  941                 |
| [`CAD2.low.fasta`](CAD2.low.fasta)       | 65.38%          |  1264                |
| [`CAD2.min.fasta`](CAD2.min.fasta)       | 55.49%          |  1625                |

Alignment of all protein sequences [`CAD2.pep.aln`](CAD2.pep.aln) shows expected 100% identity.
# Preparing to generate repair templates

Add nucleotide sequences upstream and downstream of the original genome location used for chimeragenesis.  i.e. when chimerizing PCA1 with CAD2 allele, it will be done around the original PCA1 site. This will require two files, `PCA1.upstream` and `PCA1.downstream`. 

I took the `Genomic DNA +/- 1kb` fasta from [SGD](https://www.yeastgenome.org/locus/S000000499#sequence). The first 1000 nucleotides made `PCA1.upstream` and last 1000 nucleotides made `PCA1.downstream`.

# Copy required files to this directory
```
cp ../codon_shuffle/CAD2.min.fasta ./CAD2.min.cds.fasta
cp ../codon_shuffle/CAD2.low.fasta ./CAD2.low.cds.fasta
cp ../codon_shuffle/CAD2.medium.fasta ./CAD2.medium.cds.fasta
cp ../codon_shuffle/CAD2.high.fasta ./CAD2.high.cds.fasta
cp ../codon_shuffle/CAD2.cds.fasta .
```

## Generate Repair Templates that show the method works
Using longest repair templates (80 bp for each, or 160 bp total) and most diverged CAD2 sequence
```
mkdir -p 01_test_method
python3 ./chimera.py BY_PCA1 PW5_PCA1.min --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/BY-PW5.min.RT-160.lax.fasta
python3 ./chimera.py PW5_PCA1.min BY_PCA1 --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/PW5-BY.min.RT-160.lax.fasta

python3 ./chimera.py BY_PCA1 PW5_PCA1.min --strict --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/BY-PW5.min.RT-160.strict.fasta
python3 ./chimera.py PW5_PCA1.min BY_PCA1 --strict --flanking BY_PCA1 --repair-template-length 160 --unique protein > 01_test_method/PW5-BY.min.RT-160.strict.fasta
# 262 rts each * 2 = 524
```

## Look at variety of repair template lengths and sequence homology
This will be a total of 10 transformations (PCA1-CAD2 and CAD2-PCA1 orientations, with 5 levels of sequence homology)
```
mkdir -p 02_RT_length
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.min.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.min :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.low.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.low :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.medium.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.medium :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.high.RT-all.lax.fasta  ::: BY_PCA1 ::: PW5_PCA1.high :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.orig.RT-all.lax.fasta ::: BY_PCA1 ::: PW5_PCA1.max ::: 40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.min_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.min ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.low_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.low ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.medium_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.medium ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.high_BY_PCA1.RT-all.lax.fasta  ::: PW5_PCA1.high ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160

parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.min.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.min :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.low.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.low :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.medium.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.medium :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.high.RT-all.strict.fasta  ::: BY_PCA1 ::: PW5_PCA1.high :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/BY_PCA1_PW5_PCA1.orig.RT-all.strict.fasta ::: BY_PCA1 ::: PW5_PCA1.max ::: 40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.min_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.min ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.low_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.low ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.medium_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.medium ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.high_BY_PCA1.RT-all.strict.fasta  ::: PW5_PCA1.high ::: BY_PCA1 :::  40 44 50 58 68 80 94 110 128 148 160
parallel -j 1 python3 ./chimera.py {1} {2} --strict --flanking BY_PCA1 --repair-template-length {3} --unique protein --primer-length 15 --oligo-length 190 --five-prime-padding gcgacaacggtttaggtgggtacgggtccccattccttatatagaaatggcatgttagatcggagcttccaaatcacgat --three-prime-padding accgttctttgttggaagaatagctaagcgcagggacttcccgaatctcggtattatcccggtaagtgtggactatattt  > 02_RT_length/PW5_PCA1.orig_BY_PCA1.RT-all.strict.fasta ::: PW5_PCA1.max ::: BY_PCA1 ::: 40 44 50 58 68 80 94 110 128 148 160

```

## Look at the affects of chimerizing at synonymous sites
```
mkdir -p 03_synonymous_RT
python3 ./chimera.py  BY_PCA1 PW5_PCA1.min --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/BY-PW5.min.RT-160-syn.lax.fasta
python3 ./chimera.py  PW5_PCA1.min BY_PCA1 --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/PW5-BY.min.RT-160-syn.lax.fasta

python3 ./chimera.py  BY_PCA1 PW5_PCA1.min --strict --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/BY-PW5.min.RT-160-syn.strict.fasta
python3 ./chimera.py  PW5_PCA1.min BY_PCA1 --strict --flanking BY_PCA1 --repair-template-length 160 --unique dna --oligo-length 190 --repair-template-length 160 > 03_synonymous_RT/PW5-BY.min.RT-160-syn.strict.fasta

# 1198 RTs each * 2 = 2396
```

## Generate chimeras between AcrIIa2 and AcrIIa2b (no codon shuffling needed)
```

```

## Generate all possible chimeras (no alignment needed) for AcrIIa2b and AcrIIa4
```
python3 chimera.py AcrIIa2b AcrIIa4 --all --unique dna --flanking acr --repair-template-length 160 --primer-length 15 --oligo-length 190 > 04_Acrs/AcrIIa2b_AcrIIa4.RT-160.fasta
python3 chimera.py AcrIIa4 AcrIIa2b --all --unique dna --flanking acr --repair-template-length 160 --primer-length 15 --oligo-length 190 > 04_Acrs/AcrIIa4_AcrIIa2b.RT-160.fasta
# 23,056
```

## Add skpp15 primers

Retrieved `skpp15` primers. PCR primer pairs (15-mers) obtained directly from Sri Kosuri (@skosuri). Modified from Elledge barcodes, https://elledge.hms.harvard.edu/?page_id=638

```
wget https://github.com/lasersonlab/pribar/raw/master/skpp15-forward.faa
wget https://github.com/lasersonlab/pribar/raw/master/skpp15-reverse.faa
```

Built reverse complement of of the reverse primers at https://www.bioinformatics.org/sms2/rev_comp.html and removed blank lines with `sed -i '/^$/d' skpp15-reverse-complemented.faa` 

Convert FASTA format to single-line csv

```
for file in $(ls 01_test_method/*.fasta 02_RT_length/*.fasta 03_synonymous_RT/*.fasta); do
    cat ${file} | tr "\n" "@" | sed 's/@>/\n/g' | sed 's/each@/each,/g' | tr -d "@" > ${file%.fasta}.csv
done
```

Concatenate skpp15 primers to RT sequences

```
N=100
for file in $(ls 01_test_method/*.csv 02_RT_length/*.csv 03_synonymous_RT/*.csv); do
    let N=N+1
    let lineNo=N*2
    forward=$(sed -n ${lineNo}p skpp15-forward.faa)
    reverse=$(sed -n ${lineNo}p skpp15-reverse-complemented.faa)
    awk -F "," -v f=${forward} -v r=${reverse} '{print $1,f$2r}' ${file} > ${file%.csv}.skpp${N}.RT.csv
done
```

Zip skpp15 RT oligos

```
find . | grep skpp | grep csv | zip skpp15.RTs.zip -@

```