#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

library(ensembldb)
library(dplyr)
library(data.table)

det_strand <- function(x) {
  if (x == "+") {
    return("1")
  }
  if (x == "-") {
    return("-1")
  }
  if (x == 1) {
    return("+")
  }
  if (x == -1) {
    return("-")
  }
}

split_encoding <- function(coded_vector) {
  split1 <- strsplit(coded_vector, "_")
  strand <- split1[[1]][2]
  split2 <- strsplit(split1[[1]][1], ":")
  chr <- substr(split2[[1]][1], 4, 5)
  split3 <- strsplit(split2[[1]][2], "-")
  start <- split3[[1]][1]
  end <- split3[[1]][2]
  return(paste(c(chr, start, end, strand), collapse = ":"))
}
# INPUT DATA
miranda_data <- data.frame(read.table(args[1], sep = "\t", header = T))
gtf_db <- ensembldb::EnsDb(args[2])

# MAKE GENOMIC RANGES
targets <- miranda_data$Target
targets <- targets[!duplicated(targets)]
targets <- lapply(targets, split_encoding)

targets.gr <- lapply(targets, function(x) {result=strsplit(x, ":")}) %>%
  unlist %>%
  matrix(ncol = 4, byrow = T) %>%
  as.data.frame %>%
  select(chrom=V1, start=V2, end=V3, strand=V4) %>%
  filter(complete.cases(.)) %>%
  makeGRangesFromDataFrame

genes <- genes(gtf_db)
# get overlaps
overlaps <- findOverlaps(targets.gr, genes, select = "first")
overlap_df <- as.data.frame(genes, row.names = NULL)
targets_df <- as.data.frame(targets.gr)
targets_df[c("gene_id", "gene_symbol")] <- overlap_df[overlaps, c("gene_id", "gene_name")]
f <- rownames(targets_df[is.na(targets_df$gene_symbol),])
targets_df[f, "gene_id"] <- targets_df[f, "gene_symbol"]
targets_df <- targets_df[complete.cases(targets_df),]
targets_df <- data.table(Target = paste0("chr", targets_df$seqnames, ":", targets_df$start, "-", targets_df$end, "_", targets_df$strand),
                         gene_id = targets_df$gene_id, gene_name = targets_df$gene_symbol)
# ANNOTATE DATA WITH ENSEMBL AND GENE NAME
annotated_pos <- match(miranda_data$Target, targets_df$Target)
miranda_data$gene_id <- targets_df[annotated_pos,]$gene_id
miranda_data$gene_symbol <- targets_df[annotated_pos,]$gene_name
# create contigency table
write.table(as.data.frame.matrix(table(miranda_data$gene_id, miranda_data$miRNA)), file = "miranda_counts_sponge.tsv", sep = "\t")
