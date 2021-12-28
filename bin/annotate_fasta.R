#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(seqinr)

fasta <- read.fasta(file = args[1], seqtype = "DNA")
file.name <- "circRNAs.fa"
# annotate if given
if (length(args)==2 & file.exists(args[2])) {
  circ.annotation <- read.table(args[2], header = T, sep = "\t")
  circ.annotation.map <- data.table(Target=paste0(circ.annotation$chr, ":", circ.annotation$start, "-", circ.annotation$stop, "_", circ.annotation$strand),
                                    ID=circ.annotation$circRNA.ID)
  filtered_fasta <- fasta[circ.annotation.map$Target]
  names(filtered_fasta) <- circ.annotation.map$ID
  fasta <- filtered_fasta
  file.name <- "circRNAs_annotated.fa"
}
write.fasta(sequences = fasta, names = names(fasta), nbchar = 60, file.out = file.name)
