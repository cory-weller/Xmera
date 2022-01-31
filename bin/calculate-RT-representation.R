#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)

fiveFOAfilename <- 'data/processed/5FOApool-splitreads.tab'
glucoseFilename <- 'data/processed/glucosepool-splitreads.tab'
galactoseFilename <- 'data/processed/galactosepool-splitreads.tab'


splitreadsHeader <- c('plasmidL', 'barcodeL', 'sublibraryL', 'RT', 'sublibraryR', 'barcodeR', 'plasmidR')


glu <- fread(glucoseFilename)
setnames(glu, splitreadsHeader)
glu[, type := "GLU"]

gal <- fread(galactoseFilename)
setnames(gal, splitreadsHeader)
gal[, type := "GAL"]

fiveFOA <- fread(fiveFOAfilename)
setnames(fiveFOA, splitreadsHeader)
fiveFOA[, type := "fiveFOA"]

dat <- rbindlist(list(glu, gal, fiveFOA))

rm(glu)
rm(gal)
rm(fiveFOA)
gc()

# Number of reads per sample:
dat[, .N, by=type]
#       type       N
# 1:     GLU 3358450
# 2:     GAL 3612224
# 3: fiveFOA 2545790

# Calculate number of reads per sample
GLU_reads <- dat[, .N, by=type][type=='GLU', N]
GAL_reads <- dat[, .N, by=type][type=='GAL', N]
fiveFOA_reads <- dat[, .N, by=type][type=='fiveFOA', N]


# Exclude any seqs with N in repair template
dat <- dat[! RT %like% 'N']

# Remaining reads per sample:
dat[, .N, by=type]
#       type       N
# 1:     GLU 3357347
# 2:     GAL 3611067
# 3: fiveFOA 2544950



# Exclude any seqs where barcodeL doesn't end with 'GGTT' or barcodeR doesn't begin with 'AACC'
dat <- dat[barcodeL %like% 'GGTT$' & barcodeR %like% '^AACC']


dat[, .N, by=type]
#       type       N
# 1:     GLU 3271184
# 2:     GAL 3527316
# 3: fiveFOA 2489438



# Calculate sequenced repair template length
dat[, RTlength := nchar(RT)]

dat[RTlength==163, .N, by=type]
#       type       N
# 1:     GLU 1598595
# 2:     GAL 1844764
# 3: fiveFOA 1297521

dat[, type := factor(type, levels=c("GLU", "GAL", "fiveFOA"))]

# Filter to only include 163 bp RTs
dat <- dat[RTlength == 163]

# Remove unnecessary cols to cut down on mem use
dat[, c('plasmidL', 'plasmidR', 'sublibraryL') := NULL]
dat[, sublibraryR := as.numeric(as.factor(sublibraryR))]



# Load printed oligos
printedOligos <- fread('../ERCC4.RT.txt', col.names=c('ID','oligo'))
printedOligos[, RTlength := nchar(oligo) - 30]


# Only concern ourselves with 163 bp oligos for now:
printedOligos <- printedOligos[RTlength == 163]
printedOligos[, RT := substring(oligo, 16, 178)]
printedOligos[, RT := toupper(RT)]



# only includ perfect matches
dat <- dat[RT %in% printedOligos[,RT]]

dat[, .N, by=type][type==GLU
#       type       N
# 1:     GLU 1145298
# 2:     GAL 1336999
# 3: fiveFOA  933125
# ~72% of 163 bp oligos are perfect match to 'complete list' of RTs

dat[, .N, by=list(RT, type)]

dat[, .N, by=list(type, RT)]

tmp <- dat[, .N, by=list(type,RT)]
o <- foreach(threshold=c(1, 10, 100, 1000), .combine='rbind') %do% {
    tmp2 <- tmp[N >= threshold][, .N, by=type]
    tmp2[, 'threshold' := threshold]
    return(tmp2)
}

o[, threshold := factor(threshold)]

o[type=="GLU", N := N*(1e6/GLU_reads)]
o[type=="GAL", N := N*(1e6/GAL_reads)]
o[type=="fiveFOA", N := N*(1e6/fiveFOA_reads)]

ggplot(o, aes(x=type, y=N, color=threshold, group=threshold)) + geom_point() + geom_line() +
labs(x="sample", y="Unique RTs", color="minimum_threshhold") +
theme_few(12)


perfMatch[, homologyL := substr(RT, 1, 80)]
perfMatch[, homologyR := substr(RT, 84, 163)]
perfMatch[, arms := paste0(homologyL, homologyR)]









#dat[, RTlengthBin := cut(RTlength, breaks=seq(0,250,5))]

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


# abundance of 6158 unique RTs observed in 5FOA pool:
RTdistribution <- perfMatch[, .N, by=RT][order(N)]
ggplot(RTdistribution, aes(x=1:6158, y=N)) + geom_bar(stat='identity', position='dodge') + scale_y_log10()

# abundance of 931 unique homology arms (collapsing together all codons for a given position)
armsDistribution <- perfMatch[, .N, by=list(arms)][order(N)]
ggplot(armsDistribution, aes(x=1:931, y=N)) + geom_bar(stat='identity', position='dodge') + scale_y_log10()