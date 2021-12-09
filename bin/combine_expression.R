#!/usr/bin/env Rscript

library(ensembldb)
library(tximport)

args = commandArgs(trailingOnly = TRUE)

# read gene expression counts
samples_loc <- args[1]

# metadata
samplesheet <- read.table(file = args[2], sep = "\t", header = TRUE)

# get quant files
quant.files <- file.path(samples_loc, samplesheet$sample, "salmon", "quant.sf")
names(quant.files) <- samplesheet$sample

# read input gtf
txdb <- ensembldb::ensDbFromGtf(args[3])
tx <- ensembldb::EnsDb(txdb)
tx2gene <- ensembldb::transcripts(tx, return.type = "data.frame", columns = c("gene_name", "gene_id"))
# failures <- rownames(tx2gene[is.na(tx2gene$gene_name),])
# tx2gene[failures, "gene_name"] <- tx2gene[failures, "gene_id"]
tx2gene <- tx2gene[, c("tx_id", "gene_name")]
# txi object
txi <- tximport::tximport(quant.files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T)
# write total gene expression over samples to file
counts <- as.data.frame(txi$counts)
colnames(counts) <- samplesheet$sample
write.table(counts, file = "gene_expression.tsv", sep = "\t")
saveRDS(txi, file = "txi.RDS")