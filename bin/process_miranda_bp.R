#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

library(data.table)

miranda.bp = read.table(args[1], header = T, sep = "\t")
if (length(args)==2) {
  circ.annotation <- read.table(args[2], header = T, sep = "\t")
  circ.annotation.map <- data.table(Target=paste0(circ.annotation$chr, ":", circ.annotation$start, "-", circ.annotation$stop, "_", circ.annotation$strand),
                                    ID=circ.annotation$circRNA.ID)
  miranda.bp$ID <- circ.annotation.map[match(miranda.bp$Target, circ.annotation.map$Target), "ID"]
  f <- is.na(miranda.bp$ID)
  miranda.bp$ID[f] <- miranda.bp$Target[f]
  miranda.bp$Target <- miranda.bp$ID
}
# create contigency table
write.table(as.data.frame.matrix(table(miranda.bp$Target, miranda.bp$miRNA)), file = "miranda_counts_sponge.tsv", sep = "\t")
