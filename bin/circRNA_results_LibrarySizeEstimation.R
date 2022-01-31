#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}

expression_raw_path = args[1]
samples <- gsub("-", ".", read.table(file = args[2], sep = "\t", header = T)[,"sample"]) # get samples from samplesheet
output_dir = args[3]

suppressWarnings(library(DESeq2, data.table))
'%!in%' <- function(x,y)!('%in%'(x,y))

expression_raw <- read.table(expression_raw_path, sep = "\t", header=T, stringsAsFactors = F, check.names = F)
row.names(expression_raw) <- paste0(expression_raw$chr, ":", expression_raw$start, "-", expression_raw$stop,"_", expression_raw$strand, "!", expression_raw$gene_symbol)
circRNA_names <- expression_raw[,colnames(expression_raw) %!in% samples]
circRNA_names$order <- 1:nrow(circRNA_names)
expression_raw <- expression_raw[,samples]

# normalize using DeSeq2, Library Size Estimation
meta <- data.frame(samples)
row.names(meta) <- meta$samples 
data <- expression_raw
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = data + 1, colData = meta, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)

normalized_counts <- DESeq2::counts(dds, normalized=TRUE)

# add circRNA position back to counts table
merged_data <- merge(circRNA_names, normalized_counts, by = "row.names")
merged_data <- merged_data[order(merged_data$order), ]
normalized_data <- subset(merged_data, select = -c(order, Row.names))

write.table(normalized_data, paste0("circRNA_counts_normalized.tsv"), quote = F, sep = "\t", row.names = F)
