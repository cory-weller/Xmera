#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)

fiveFOAfilename <- '5FOApool.splitreads.tab'
glucoseFilename <- 'glucosepool.splitreads.tab'

dat <- fread(RTfilename, header=FALSE)
splitreadsHeader <- c('plasmidL', 'barcodeL', 'sublibraryL', 'RT', 'sublibraryR', 'barcodeR', 'plasmidR')
setnames(dat, splitreadsHeader)

# 2.48 M (2484722) Reads with perfect sublibrary seqs

# Exclude any seqs with N in repair template
dat <- dat[! RT %like% 'N']
# 2.48 M (2483901) Reads remaining

# Exclude any seqs where barcodeL doesn't end with 'GGTT' or barcodeR doesn't begin with 'AACC'
dat <- dat[barcodeL %like% 'GGTT$' & barcodeR %like% '^AACC']
# 2.42 M (2429735) Reads remaining

# Calculate sequenced repair template length
dat[, RTlength := nchar(RT)]
dat[, RTlengthBin := cut(RTlength, breaks=seq(0,250,5))]

# Log10 transformation
ggplot(dat[, .N, by=RTlengthBin], aes(x=RTlengthBin, y=N)) +
geom_bar(stat='identity', position='dodge') +
scale_y_log10() + 
labs(   x='Repair Template Length',
        y='Log10(Frequency)',
        title='Repair Template Length Distribution'
    ) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Normal
ggplot(dat[, .N, by=RTlengthBin], aes(x=RTlengthBin, y=N)) +
geom_bar(stat='identity', position='dodge') +
labs(   x='Repair Template Length',
        y='Frequency',
        title='Repair Template Length Distribution'
    ) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Calculate inteded repair template length
printedOligos <- fread('ERCC4.RT.txt', col.names=c('ID','oligo'))

# Calculate inteded RT length
printedOligos[, RTlength := nchar(oligo) - 30]

nrow(dat)
nrow(dat[RTlength == 163])
nrow(dat[RTlength == 163]) / nrow(dat)

# Make new length bins for better visualization
dat[, RTlengthBin := cut(RTlength, breaks=c(0,160, 161, 162, 163, 164, Inf))]
dat[, RTlengthBin := as.numeric(RTlengthBin)]
dat[RTlengthBin == 1, RTlengthLabel := '<=160']
dat[RTlengthBin == 2, RTlengthLabel := '161']
dat[RTlengthBin == 3, RTlengthLabel := '162']
dat[RTlengthBin == 4, RTlengthLabel := '163']
dat[RTlengthBin == 5, RTlengthLabel := '164']
dat[RTlengthBin == 6, RTlengthLabel := '>=165']
dat[, RTlengthLabel := factor(RTlengthLabel, levels=c('<=160','161','162','163','164','>=165'))]


# Log10 transformation of modified bins
ggplot(dat[, .N, by=RTlengthLabel], aes(x=RTlengthLabel, y=N)) +
geom_bar(stat='identity', position='dodge') +
labs(   x='Repair Template Length',
        y='Frequency',
        title='Sequenced Repair Template Length Distribution'
    ) +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Only concern ourselves with 163 bp oligos for now:
printedOligos <- printedOligos[RTlength == 163]
printedOligos[, RT := substring(oligo, 16, 178)]
printedOligos[, RT := toupper(RT)]

# only includ perfect matches
perfMatch <- dat[RT %in% printedOligos[,RT]]
perfMatch[, homologyL := substr(RT, 1, 80)]
perfMatch[, homologyR := substr(RT, 84, 163)]
perfMatch[, arms := paste0(homologyL, homologyR)]

perfMatch[, .N, by=list(RT)]

# abundance of 6158 unique RTs observed in 5FOA pool:
RTdistribution <- perfMatch[, .N, by=RT][order(N)]
ggplot(RTdistribution, aes(x=1:6158, y=N)) + geom_bar(stat='identity', position='dodge') + scale_y_log10()

# abundance of 931 unique homology arms (collapsing together all codons for a given position)
armsDistribution <- perfMatch[, .N, by=list(arms)][order(N)]
ggplot(armsDistribution, aes(x=1:931, y=N)) + geom_bar(stat='identity', position='dodge') + scale_y_log10()