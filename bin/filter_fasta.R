#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(BSgenome, Biostrings, GenomicRanges, GenomicFeatures, seqinr, rtracklayer)

# fasta, gtf, circRNA hits
args = commandArgs(trailingOnly = TRUE)

# load fasta
print("loading fasta...")
fasta <- readDNAStringSet(args[1])
# load circRNA hits
print("loading circRNA hits...")
circ.expression <- read.table(args[2], sep = "\t", header = T)
# save subset of fasta entries
fasta <- fasta[rownames(circ.expression)]
# save to file
writeXStringSet(fasta, filepath = "circRNAs_filtered.fa")