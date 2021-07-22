# Generating codon-shuffled versions of the CAD2 allele
# README

## Adding primer sequences
the script `addPrimers.py` will add sublibrary primer sequences and format as a tab-delimited
two-column file. The first column will contain the repair template header plus the numeric ID of
the skpp primer used. The second column will contain the full oligonucleotide sequence.

```bash
python3 Xmera/bin/addPrimers.py <RT FASTA> <skpp primer #>
```