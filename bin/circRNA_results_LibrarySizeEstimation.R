#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("Three argument must be supplied", call.=FALSE)
}

expression_raw_path = args[1]
samplesheet <- read.table(file = args[2], header = T) # get samples from samplesheet
samples <- samplesheet$sample
output_dir = args[3]

suppressWarnings(library(DESeq2, data.table))
'%!in%' <- function(x,y)!('%in%'(x,y))

expression_raw <- read.table(expression_raw_path, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
expression_raw$key <- paste0(expression_raw$chr, ":", expression_raw$start, "-", expression_raw$stop,"_", expression_raw$strand)
expression_raw <- expression_raw[!duplicated(expression_raw$key),]
row.names(expression_raw) <- expression_raw$key
circRNA_names <- expression_raw[,colnames(expression_raw) %!in% samples]
circRNA_names$order <- 1:nrow(circRNA_names)
expression_raw <- expression_raw[,samples]
expression_raw[is.na(expression_raw)] <- 0

# normalize using DeSeq2, Library Size Estimation
meta <- data.frame(samples)
row.names(meta) <- meta$samples 
data <- as.matrix(expression_raw)
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(data + 1), colData = meta, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)

normalized_counts <- DESeq2::counts(dds, normalized=TRUE)
# add circRNA data to counts
normalized_data <- cbind(circRNA_names, normalized_counts)
normalized_data <- subset(normalized_data, select = -c(order, key))
write.table(normalized_data, paste0("circRNA_counts_normalized.tsv"), quote = F, sep = "\t")
