#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}

raw_counts_path = args[1]
output_dir = args[2]

library(DESeq2)
library(data.table)
 
raw_counts <- read.table(raw_counts_path, sep = "\t",
                         header = T, stringsAsFactors = F,  check.names = F)
samples <- colnames(raw_counts)[-c(1)]

row.names(raw_counts) <- raw_counts$miRNA
data <- copy(raw_counts[,-c(1)])
miRNA_names <- data.frame(miRNA = raw_counts$miRNA, order = 1:nrow(raw_counts))


# normalize using DeSeq2, Library Size Estimation
meta <- data.frame(samples)
row.names(meta) <- meta$samples 
all(colnames(data) %in% rownames(meta))
all(colnames(data) == rownames(meta))

dds <- DESeq2::DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ 1)
dds <- DESeq2::estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- DESeq2::counts(dds, normalized=TRUE)

# add miRNA IDs back to counts table
merged_data <- merge(miRNA_names, normalized_counts, by.x = "miRNA", by.y = "row.names")
merged_data <- merged_data[order(merged_data$order), ]
norm_data <- subset(merged_data, select = -c(order))

write.table(norm_data, paste0("miRNA_counts_all_samples_libSizeEstNorm.tsv"), quote = F, sep = "\t", row.names = F)


