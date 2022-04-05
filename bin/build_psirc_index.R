#!/usr/bin/env Rscript

library(Biostrings)
library(argparser)

args <- commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
parser <- add_argument(parser, "--index", help = "Psirc index file name/file location", default = "psirc.index")
parser <- add_argument(parser, "--psirc_quant", help = "Path to psirc-quant executable", default = "psirc-quant")
# ONLY NEEDED IF INDEX IS NOT ALREADY CONSTRUCTED
parser <- add_argument(parser, "--transcriptome", help = "Transcriptome for given genome (cDNA)", default = NULL)
parser <- add_argument(parser, "--circ_fasta", help = "Fasta file for circRNAs in circRNA counts file", default = NULL)

argv <- parse_args(parser, argv = args)

# save psirc index location
index <- argv$index

# build index if it not already exists
if (!file.exists(index)) {
  print("BUILDING PSIRC INDEX")
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
  # remove duplicates
  all.fasta <- all.fasta[!duplicated(names(all.fasta))]
  out <- "cDNA_circRNA.fa"
  # write combined file to disk
  writeXStringSet(all.fasta, filepath = out)
  # build psirc index with combined fasta file
  system(paste(c(argv$psirc_quant, "index", "-i", index, "--make-unique", out), collapse = " "))
}