#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(BSgenome, Biostrings, GenomicRanges, GenomicFeatures, seqinr, rtracklayer)


overlap <- function(exons, circRNAs) {
  # find all exons within the circRNAs and reformat result
  o <- findOverlapPairs(exons, circRNAs, type = "within")
  ids <- names(o@second)
  o <- o@first
  names(o) <- ids
  # combine overlapping parts on circRNA level
  return(reduce(split(o, ids)))
}

# fasta, gtf, circRNA hits
args = commandArgs(trailingOnly = TRUE)

# load fasta
print("loading fasta...")
fasta <- readDNAStringSet(args[1])
# load gtf
print("loading gtf...")
gtf <- readGFFAsGRanges(args[2])
# load circRNA hits
print("loading circRNA hits...")
circ.expression <- read.table(args[3], sep = "\t", header = T)
# build Granges
ci.expression <- circ.expression[circ.expression$type == "ciRNA",]
circ.expression <- circ.expression[circ.expression$type == "circRNA",]
circ.ranges <- makeGRangesFromDataFrame(circ.expression[,1:4])
ci.ranges <- makeGRangesFromDataFrame(ci.expression[,1:4])
# set names according to row names
names(circ.ranges) <- row.names(circ.expression)
names(ci.ranges) <- row.names(ci.expression)

# extract exons
exons <- gtf[gtf$type == "exon",]
# extract introns
introns <- gaps(exons)

# splice sequences
print("splicing circRNAs...")
# find exonic circRNA parts that match with gtf
print("exonic circRNAs")
circ.exons <- overlap(exons, circ.ranges)
# get exonic sequences
circ.sequences <- DNAStringSet(sapply(getSeq(fasta, circ.exons), unlist))
# get intronic sequences
print("intronic circRNAs")
circ.introns <- overlap(introns, ci.ranges)
ci.sequences <- DNAStringSet(sapply(getSeq(fasta, circ.introns), unlist))
# combine sequences
all.sequences <- c(circ.sequences, ci.sequences)
print("writing final circRNA fasta")
# write sequences to file
writeXStringSet(all.sequences, filepath = "circRNAs.fa")
