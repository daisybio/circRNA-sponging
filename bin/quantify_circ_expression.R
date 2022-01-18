#!/usr/bin/env Rscript

# run psirc for every sample counts

args <- commandArgs(trailingOnly = TRUE)

# save psirc index location
index <- args[1]
# read samplesheet
samplesheet <- read.table(args[2], header = T, sep = "\t")
# read circRNA counts
circ.counts <- read.table(args[3], header = T, sep = "\t")
# PARAMS
fragment.length <- 76
standard.dev <- 20
# change colnames
colnames(circ.counts) <- sapply(gsub("\\.", "-", colnames(circ.counts)), "[", 1)
# create tmp
dir.create("tmp", showWarnings = F)
# determine annotation
annotation <- "circBaseID" %in% colnames(circ.counts)
# set row names
if (annotation) {
  rownames(circ.counts) <- circ.counts$circBaseID
} else {
  rownames(circ.counts) <- paste0(circ.counts$chr, ":", circ.counts$start, "-", circ.counts$stop, "_", circ.counts$strand)
}

# TODO: implement paired read option
# run psirc for every sample
for (i in 1:nrow(samplesheet)) {
  sample <- samplesheet[i, "sample"]
  output <- file.path("tmp", sample)
  fastq <- samplesheet[i, "totalRNA1"]
  # write to tmp/sample
  cmd <- paste(c("psirc-quant quant -i", index, "-o", output, "--single", "-l", fragment.length, "-s", standard.dev, fastq), collapse = " ")
  std <- system(cmd, ignore.stdout = T, intern = T)           # invoke psirc command
  # read created abundances
  abundance <- read.table(file.path(output, "abundance.tsv"), sep = "\t", header = T)
  if (annotation) {
    abundance.circ <- abundance[grepl("circ", abundance$target_id),]
  } else {
    abundance.circ <- abundance[grepl("chr", abundance$target_id),]
  }
  # set counts to quantified levels
  circ.counts[abundance.circ$target_id, sample] <- abundance.circ$est_counts
  cat(eval(round(i/nrow(samplesheet), 2)*100), " %", "\r")
}
# remove tmp
# unlink("tmp", recursive = T)
# write output to disk
o <- paste0(strsplit(basename(args[3]), "\\.")[[1]][1], "_quant", ".tsv")
cat("writing output file to ", o, "...\n")
write.table(circ.counts, file = o, sep = "\t")
print("done")
