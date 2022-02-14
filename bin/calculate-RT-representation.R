#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(foreach)
library(ggthemes)

fiveFOAfilename <- 'data/processed/5FOApool-splitreads.tab.gz'
glucoseFilename <- 'data/processed/glucosepool-splitreads.tab.gz'
galactoseFilename <- 'data/processed/galactosepool-splitreads.tab.gz'
SNPPoolFilename <- 'data/processed/SNPs-pool-splitreads.tab.gz'
notSNPPoolFilename <- 'data/processed/not-SNPs-pool-splitreads.tab.gz'
indelPoolFilename <- 'data/processed/indel-pool-splitreads.tab.gz'
aggregateFilename <- 'data/processed/combined-splitreads.tab.gz'
librariesFilename <- 'src/libraries.csv'

splitreadsHeader.g <- c('plasmidL', 'barcodeL', 'sublibraryL', 'RT', 'sublibraryR', 'barcodeR', 'plasmidR')
colOrder.g <- c('barcodeL', 'RT', 'barcodeR', 'timepoint', 'library', 'treatment')

libs.g <- fread(librariesFilename)
setkey(libs.g, library, treatment)



processReads <- function(readsFilename, libraries, timepoint, treatments) {
    dat.all <- fread(readsFilename)
    setnames(dat.all, splitreadsHeader.g)

    output <- foreach(library.i=libraries, .combine='rbind') %do% {
        foreach(treatment.i=treatments, .combine='rbind') %do% {
            libs.i <- libs.g[.(library.i, treatment.i)]
            barcodeLlength <- libs.i[,LbarcodeLength]
            barcodeRlength <- libs.i[,RbarcodeLength] 
            dat <- merge(dat.all, libs.i, by=c('sublibraryL', 'sublibraryR'))
            dat[, c('plasmidL', 'sublibraryL', 'sublibraryR', 'plasmidR', 'primer') := NULL]
            dat[, barcodeLsplit := tstrsplit(barcodeL, split=libs.i[,'LbarcodeEnd'])[1]]
            dat[, barcodeRsplit := tstrsplit(barcodeR, split=libs.i[,'RbarcodeStart'])[2]]
            dat <- dat[nchar(barcodeLsplit) == barcodeLlength]
            dat <- dat[nchar(barcodeRsplit) == barcodeRlength]

            dat[, c('barcodeL', 'barcodeR') := NULL]
            setnames(dat, 'barcodeLsplit', 'barcodeL')
            setnames(dat, 'barcodeRsplit', 'barcodeR')
            dat[, c('LbarcodeLength','RbarcodeLength','LbarcodeEnd','RbarcodeStart') := NULL]
            dat[, 'timepoint' := timepoint]
            dat[, 'library' := library.i]
            dat[, 'treatment' := treatment.i]
            setcolorder(dat, colOrder.g)
            return(dat[])
        }
    }
    return(output)
}

if(! file.exists(aggregateFilename)) {
    # get reads at each stage of experiment
    glu <- processReads(glucoseFilename, c('SNP','not-SNP','indel'), c('GLU'), 'all')
    gal <- processReads(galactoseFilename, c('SNP','not-SNP','indel'), c('GAL'), 'all')
    fiveFOA <- processReads(fiveFOAfilename, c('SNP','not-SNP','indel'), c('fiveFOA'), 'all')

    # get plasmid pool reads
    plasmidPoolSNP <- processReads(SNPPoolFilename, 'SNP', 'plasmid-pool', 'all')
    plasmidPoolNotSNP <- processReads(notSNPPoolFilename, 'not-SNP', 'plasmid-pool', 'all')
    plasmidPoolIndel <- processReads(indelPoolFilename, 'indel', 'plasmid-pool', 'all')

    dat <- rbindlist(list(glu, gal, fiveFOA, plasmidPoolSNP, plasmidPoolNotSNP, plasmidPoolIndel)) 
    fwrite(dat, file=aggregateFilename, quote=F, row.names=F, col.names=T)
} else {
    dat <- fread(aggregateFilename)
}

dat.tmp <- dat[, .N, by=list(barcodeL, barcodeR, timepoint)][, .N, by=list(timepoint)]
dat.tmp[, timepoint := factor(timepoint, levels=c('plasmid-pool','GLU','GAL','fiveFOA'))]
g <- ggplot(dat.tmp, aes(x=timepoint, y=N)) + geom_point() +
labs(x='timepoint', y='N unique barcodes')

ggsave(g, file='uniqueBCs.png')

# Number of reads per sample:
dat[, .N, by=list(timepoint, library)]
#        timepoint library       N
#  1:          GLU     SNP 1869391
#  2:          GLU not-SNP  453316
#  3:          GLU   indel   45842
#  4:          GAL     SNP 1176513
#  5:          GAL not-SNP 1118400
#  6:          GAL   indel   45991
#  7:      fiveFOA     SNP  675307
#  8:      fiveFOA not-SNP  989537
#  9:      fiveFOA   indel  132251
# 10: plasmid-pool     SNP  836022
# 11: plasmid-pool not-SNP  546227
# 12: plasmid-pool   indel  427393


# Exclude any seqs with N in repair template
dat <- dat[! RT %like% 'N']


# Calculate sequenced repair template length
dat[, RTlength := nchar(RT)]




# Load printed oligos
printedOligos <- fread('../ERCC4.RT.txt', col.names=c('ID','oligo'))

# exclude first 15 and last 15 nucleotides
printedOligos[, endpoint := nchar(oligo) - 15]
printedOligos[, RT := substring(oligo, 16, endpoint)]
printedOligos[, RT := toupper(RT)]



# only includ perfect matches
dat <- dat[RT %in% printedOligos[,RT]]

tmp <- dat[, .N, by=list(timepoint,RT)]
o <- foreach(threshold=c(1, 10, 100, 1000), .combine='rbind') %do% {
    tmp2 <- tmp[N >= threshold][, .N, by=timepoint]
    tmp2[, 'threshold' := threshold]
    return(tmp2)
}

o[, threshold := factor(threshold)]
o[, timepoint := factor(timepoint, levels=c('plasmid-pool','GLU','GAL','fiveFOA'))]


o.notunique <- o
o.notunique[, type := 'unfiltered']

g <- ggplot(o, aes(x=timepoint, y=N, color=threshold, group=threshold)) + geom_point() + geom_line() +
labs(x="sample", y="Unique RTs", color="minimum_threshhold") +
theme_few(12) +
geom_hline(yintercept=nrow(printedOligos), linetype='dashed', alpha=0.6)

ggsave(g, file='RTrepresentation.png')




dat.cts <- dat[, .N, by=list(RT, barcodeL, barcodeR, timepoint)][order(N)]
dat.cts[, xval := 1:.N, by=list(timepoint)]
dat.cts[, timepoint := factor(timepoint, levels=c('plasmid-pool','GLU','GAL','fiveFOA'))]
g2 <- ggplot(dat.cts, aes(x=xval, y=N)) + geom_line() + scale_y_log10() + facet_grid(.~timepoint) +
labs(x='unique combinations of RT and barcode', title='A small number of barcodes dominate the pools') +
theme_few(12)
ggsave(g2, file='RTtakeover.png')





# same as above, but filtering for unique barcodes
dat <- unique(dat)

tmp <- dat[, .N, by=list(timepoint,RT)]
o <- foreach(threshold=c(1, 10, 100, 1000), .combine='rbind') %do% {
    tmp2 <- tmp[N >= threshold][, .N, by=timepoint]
    tmp2[, 'threshold' := threshold]
    return(tmp2)
}

o[, threshold := factor(threshold)]
o[, timepoint := factor(timepoint, levels=c('plasmid-pool','GLU','GAL','fiveFOA'))]

o.unique <- o
o.unique[, type := 'unique-filtered-barcode']

o <- rbindlist(list(o.unique, o.notunique))

g <- ggplot(o, aes(x=timepoint, y=N, color=threshold, group=threshold)) + geom_point() + geom_line() +
labs(x="", y="RTs sampled at least (threshold) times", color="threshold") +
theme_few(12) +
geom_hline(yintercept=nrow(printedOligos), linetype='dashed', alpha=0.6) +
facet_grid(.~type)

ggsave(g, file='RTrepresentation.png')



















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
