#!/usr/bin/env Rscript

# run psirc for every sample counts

args <- commandArgs(trailingOnly = TRUE)

# save psirc index location
index <- args[1]
# read samplesheet
samplesheet <- read.table(args[2], header = T, sep = "\t")
# read circRNA counts
circ.counts <- read.table(args[3], header = T, sep = "\t")
# change colnames
colnames(circ.counts) <- sapply(gsub("\\.", "-", colnames(circ.counts)), "[", 1)
# create tmp
dir.create("tmp", showWarnings = F)

# TODO: implement paired read option
# run psirc for every sample
for (i in 1:nrow(samplesheet)) {
  sample <- samplesheet[i, "sample"]
  fastq <- samplesheet[i, "totalRNA1"]
  cmd <- paste(c("psirc-quant quant -i", index, "-o", "./tmp", "--single", "-l", "15", "-s", "20", fastq), collapse = " ")
  cat(cmd)
  # system(cmd)
}