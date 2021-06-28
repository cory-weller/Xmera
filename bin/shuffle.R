#!/usr/bin/env Rscript
# This script will shuffle fasta_2 to be various levels of similarity
# to the unmodified fasta_1 sequence

args <- commandArgs(trailingOnly=TRUE)

# set freq_filter to exclude codons below a certain (amino-acid specific) frequency
freq_filter <- 0.0

fasta_1_filename <- args[1]
fasta_2_filename <- args[2]
splitName <- strsplit(fasta_2_filename, split="\\.")[[1]]
fileStem <- paste0(splitName[1:(length(splitName)-1)], collapse=".")


codons_filename <- "codons.tab"

fasta_1 <- readLines(fasta_1_filename)
fasta_1 <- fasta_1[2:length(fasta_1)]
fasta_1 <- paste0(fasta_1, collapse="")
fasta_1 <- strsplit(fasta_1, "(?<=.{3})", perl = TRUE)[[1]]

fasta_2 <- readLines(fasta_2_filename)
fasta_2 <- fasta_2[2:length(fasta_2)]
fasta_2 <- paste0(fasta_2, collapse="")
fasta_2 <- strsplit(fasta_2, "(?<=.{3})", perl = TRUE)[[1]]

library(data.table)
library(foreach)

dat <- fread(codons_filename, col.names=c("codon", "AA", "fraction_per_AA", "count_per_thousand", "count"))
dat[, codon := gsub("U", "T", codon)]

# recalculate fraction_per-AA with greater number of digits
dat[, fraction_per_AA := count/sum(count), by=AA]

# Remove codons with fraction_per_AA lower than <freq_filter>%
dat <- dat[fraction_per_AA > freq_filter]
setkey(dat, codon, AA)

AAs <- unique(dat[,AA])
codonTable <- unique(dat[, c("codon","AA")])

o <- CJ("codon1"=dat[, codon], "codon2"=dat[, codon])

dat2 <- merge(o, dat, by.x="codon1", by.y="codon")
setnames(dat2, "codon1", "codon")

dat2[, "aa1" := ifelse(substr(codon, 1, 1) == substr(codon2, 1, 1), FALSE, TRUE)]
dat2[, "aa2" := ifelse(substr(codon, 2, 2) == substr(codon2, 2, 2), FALSE, TRUE)]
dat2[, "aa3" := ifelse(substr(codon, 3, 3) == substr(codon2, 3, 3), FALSE, TRUE)]
dat2[, dist := aa1 + aa2 + aa3]
dat2[, c("aa1","aa2","aa3") := NULL]

setkey(dat2, codon, AA)
dat <- dat2[fraction_per_AA >= 0.1]

setnames(dat2, "codon", "new_codon")
setnames(dat2, "codon2", "old_codon")

translate <- function(i) {
    return(codonTable[codon==i, AA])
}

# shuffle sequence 2 to be <method> divergent from 1.
# method = "max" will return maximally divergent sequence
# method = "min" will return minimally divergent sequence (i.e. as similar as possible on nucleotide level)
codonShuffle <- function(seq_1_codon_static, seq_2_codon_to_shuffle, method) {
    dt <- dat2[AA %in% translate(seq_2_codon_to_shuffle) & old_codon == seq_1_codon_static]
    if(nrow(dt) == 0) {
        return(seq_2_codon_to_shuffle)
    } else {
        if(method == "max") {
            return(dt[which(dist==max(dist))][sample(.N, size=1, prob=fraction_per_AA), new_codon])
        } else if(method == "high") {
            return(dt[sample(.N, size=1, prob=((1 + dist)^1.7 * fraction_per_AA)), new_codon])
        } else if(method == "min") {
            return(dt[which(dist==min(dist))][sample(.N, size=1, prob=fraction_per_AA), new_codon])
        } else if(method == "low") {
            return(dt[sample(.N, size=1, prob=((4 - dist)^4 * fraction_per_AA)), new_codon])
        } else if(method == "medium") {
            return(dt[sample(.N, size=1, prob=fraction_per_AA), new_codon])
        } else {
            stop("Invalid method. Use method = 'max' 'high' 'medium' 'low' or 'min'")
        }
    }
}

getStringDistance <- function(a, b) {
    a <- strsplit(a, split = "")[[1]]
    b <- strsplit(b, split = "")[[1]]
    # Returns length of differences between A and B if they are of same length
    if(length(a) != length(b)) {
        stop("strings are not of equal length!")
    } else {
        return(sum(a != b))
    }
}

fasta_1_collapsed <- paste0(fasta_1, collapse="")
fasta_2_collapsed <- paste0(fasta_2, collapse="")

set.seed(1)

shuffled_max <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="max"), collapse="")
shuffled_min <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="min"), collapse="")
shuffled_low <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="low"), collapse="")
shuffled_high <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="high"), collapse="")
shuffled_medium <- paste0(mapply(codonShuffle, seq_1_codon_static = fasta_1, seq_2_codon_to_shuffle = fasta_2, method="medium"), collapse="")

# n_min <- getStringDistance(shuffled_min, fasta_1_collapsed)
# percent_min <- format(100 * (1 - (n_min / nchar(fasta_1_collapsed))), digits = 5, format = "f")

# n_low <- getStringDistance(shuffled_low, fasta_1_collapsed)
# percent_low <- format(100 * (1 - (n_low / nchar(fasta_1_collapsed))), digits = 5, format = "f")

# n_medium <- getStringDistance(shuffled_medium, fasta_1_collapsed)
# percent_medium <- format(100 * (1 - (n_medium / nchar(fasta_1_collapsed))), digits = 5, format = "f")

# n_high <- getStringDistance(shuffled_high, fasta_1_collapsed)
# percent_high <- format(100 * (1 - (n_high / nchar(fasta_1_collapsed))), digits = 5, format = "f")

# n_max <- getStringDistance(shuffled_max, fasta_1_collapsed)
# percent_max <- format(100 * (1 - (n_max / nchar(fasta_1_collapsed))), digits = 5, format = "f")

# plot(x=0:4, y=c(n_min, n_low, n_medium, n_high, n_max))


# headers renamed here to refer to IDENTITY instead of divergence here. Max identity = min divergence.
write(paste0(">", fileStem, "_min_identity", collapse=""), file=paste0(fileStem, ".min.fasta", collapse=""))
write(paste0(">", fileStem, "_max_identity", collapse=""), file=paste0(fileStem, ".max.fasta", collapse=""))
write(paste0(">", fileStem, "_high_identity", collapse=""), file=paste0(fileStem, ".high.fasta", collapse=""))
write(paste0(">", fileStem, "_low_identity", collapse=""), file=paste0(fileStem, ".low.fasta", collapse=""))
write(paste0(">", fileStem, "_medium_identity", collapse=""), file=paste0(fileStem, ".medium.fasta", collapse=""))

write(shuffled_medium, file=paste0(fileStem, ".medium.fasta", collapse=""), append=TRUE)
write(shuffled_max, file=paste0(fileStem, ".min.fasta", collapse=""), append=TRUE)
write(shuffled_min, file=paste0(fileStem, ".max.fasta", collapse=""), append=TRUE)
write(shuffled_low, file=paste0(fileStem, ".high.fasta", collapse=""), append=TRUE)
write(shuffled_high, file=paste0(fileStem, ".low.fasta", collapse=""), append=TRUE)
