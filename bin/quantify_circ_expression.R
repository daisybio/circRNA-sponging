#!/usr/bin/env Rscript

library(Biostrings)
library(argparser)

args <- commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
parser <- add_argument(parser, "--index", help = "Psirc index file name/file location", default = "psirc.index")
parser <- add_argument(parser, "--samplesheet", help = "Samplesheet containing metadata for samples")
parser <- add_argument(parser, "--circ_counts", help = "circRNA counts file")
parser <- add_argument(parser, "--psirc-quant", help = "Path to psirc-quant executable")
# ONLY NEEDED IF INDEX IS NOT ALREADY CONSTRUCTED
parser <- add_argument(parser, "--transcriptome", help = "Transcriptome for given genome (cDNA)", default = NULL)
parser <- add_argument(parser, "--circ_fasta", help = "Fasta file for circRNAs in circRNA counts file", default = NULL)

argv <- parse_args(parser, argv = args)

# save psirc index location
index <- argv$index

# build index if it not already exists
if (!file.exists(index)) {
  # not all paramaters are given
  if (!file.exists(argv$transcriptome) | !file.exists(argv$circ_fasta)) {
    stop("Not all options required for index creation are given")
  }
  # read transcriptome
  transcriptome <- readDNAStringSet(argv$transcriptome)
  # read circ fasta
  circ.fasta <- readDNAStringSet(argv$circ_fasta)
  # annotate circRNA fasta with C marker
  names(circ.fasta) <- paste0(names(circ.fasta), "\tC")
  # concat circRNAs with transcriptome
  all.fasta <- c(transcriptome, circ.fasta)
  out <- "cDNA_circRNA.fa"
  # write combined file to disk
  writeXStringSet(all.fasta, filepath = out)
  # build psirc index with combined fasta file
  cmd <- paste(c("psirc-quant index -i", index, out), collapse = " ")
  system(cmd)
}
# read samplesheet
samplesheet <- read.table(argv$samplesheet, header = T, sep = "\t")
# single or paired end
single_end <- !"totalRNA2" %in% colnames(samplesheet)
# read circRNA counts
circ.counts <- read.table(argv$circ_counts, header = T, sep = "\t")
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

circ.quant <- circ.counts
# run psirc for every sample
for (i in 1:nrow(samplesheet)) {
  sample <- samplesheet[i, "sample"]
  output <- file.path("tmp", sample)
  fastq <- samplesheet[i, "totalRNA1"]
  if (single.end) {
    # write to tmp/sample
    cmd <- paste(c("psirc-quant quant -i", index, "-o", output, "--single", "-l", fragment.length, "-s", standard.dev, fastq), collapse = " ")
  } else {
    fastq2 <- samplesheet[i, "totalRNA2"]
    cmd <- paste(c("psirc-quant quant -i", index, "-o", output, fastq, fastq2), collapse = " ")
  }
  std <- system(cmd, ignore.stdout = T, intern = T)           # invoke psirc command
  # read created abundances
  abundance <- read.table(file.path(output, "abundance.tsv"), sep = "\t", header = T)
  if (annotation) {
    abundance.circ <- abundance[grepl("circ", abundance$target_id),]
  } else {
    abundance.circ <- abundance[grepl("chr", abundance$target_id),]
  }
  # set counts to quantified levels
  circ.quant[abundance.circ$target_id, sample] <- abundance.circ$est_counts
  cat(eval(round(i/nrow(samplesheet), 2)*100), " %", "\r")
}

# filter quantified data

# remove tmp
# unlink("tmp", recursive = T)
# write output to disk
o <- paste0(strsplit(basename(args[3]), "\\.")[[1]][1], "_quant", ".tsv")
cat("writing output file to ", o, "...\n")
write.table(circ.quant, file = o, sep = "\t")
print("done")

# extract counts
circ.counts.c <- circ.counts[order(circ.counts$circBaseID), -c(1:8)]
circ.quant.c <- circ.quant[order(circ.quant$circBaseID), -c(1:8)]
# calculate false positive rate
tmp <- circ.counts.c - circ.quant.c
fp <- length(tmp[tmp ==circ.counts.c])
fp.rate <- fp/prod(dim(circ.counts.c))
cat("False positive rate:", fp.rate, "according to psirc quantification\n")
