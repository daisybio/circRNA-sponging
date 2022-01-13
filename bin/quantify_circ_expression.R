#!/usr/bin/env Rscript

# run psirc for every sample counts

args = commandArgs(trailingOnly = TRUE)

# save psirc index location
index <- args[1]
# read samplesheet
samplesheet <- read.table(args[2], header = T, sep = "\t")
# read circRNA counts
circ.counts <- read.table(args[3], header = T, sep = "\t")