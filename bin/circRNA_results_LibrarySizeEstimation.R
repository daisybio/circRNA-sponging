#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}

expression_raw_path = args[1]
output_dir = args[2]

library(DESeq2)
library(data.table)

expression_raw <- read.table(expression_raw_path, sep = "\t", header=T, stringsAsFactors = F, check.names = F)
samples <- colnames(expression_raw)[-c(1:5)]
row.names(expression_raw) <- paste0(expression_raw$chr, ":", expression_raw$start, "-", expression_raw$stop,"_", expression_raw$strand, "!", expression_raw$gene_symbol)
circRNA_names <- expression_raw[,c(1,2,3,4,5)]
circRNA_names$order <- 1:nrow(circRNA_names)
expression_raw <- expression_raw[,-c(1,2,3,4,5)]

# normalize using DeSeq2, Library Size Estimation
meta <- data.frame(samples)
row.names(meta) <- meta$samples 
data <- expression_raw
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- DESeq2::counts(dds, normalized=TRUE)

# add circRNA position back to counts table
merged_data <- merge(circRNA_names, normalized_counts, by = "row.names")
merged_data <- merged_data[order(merged_data$order), ]
normalized_data <- subset(merged_data, select = -c(order, Row.names))

write.table(normalized_data, paste0(output_dir, "circRNA_counts_normalized.tsv", sep = "/"), quote = F, sep = "\t", row.names = F)


