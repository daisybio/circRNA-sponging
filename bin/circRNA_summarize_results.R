#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
output_dir = args[2]

dataset <- read.table(dataset_path, sep = "\t", header=T, stringsAsFactors = F)
samples <- dataset$sample

finaldata <- NULL
N_of_circRNAs_raw <- c()
for (i in 1:length(samples)){
  sample <- samples[i]
  path <- paste0(output_dir,"/samples/", sample , "/circRNA_detection/circExplorer2/", sample ,"_circularRNA_known.txt")
  CIRCexplorer2_output <- data.table(read.table(path, sep = "\t"))
  
  N_of_circRNAs_raw[i] <- nrow(CIRCexplorer2_output)
  names(N_of_circRNAs_raw)[i] <- sample
  
  expression_raw <- CIRCexplorer2_output[,c(1,2,3,6,13,15,16)]
  colnames(expression_raw) <- c("chr", "start", "stop", "strand", "counts", "gene_symbol", paste("isoform_", sample, sep = ""))
  expression_raw <- expression_raw[,c(1,2,3,4,5,6)]
  
  # compact and remove duplicates
  compact_raw <- data.table(circRNA=paste0(strsplit(expression_raw$chr, "_")[[1]][1], ":", expression_raw$start, "-", expression_raw$stop,"_", expression_raw$strand, "!", expression_raw$gene_symbol), counts = expression_raw$counts)
  compact_raw <- compact_raw[, max(counts), by=circRNA]
  
  compact_raw$chr <- sapply(strsplit(as.character(compact_raw$circRNA),':'), "[", 1)
  compact_raw$start <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(compact_raw$circRNA),':'), "[", 2),'-'), "[", 1))
  compact_raw$stop <- as.numeric(sapply(strsplit(sapply(strsplit(as.character(compact_raw$circRNA),'-'), "[", 2),'_'), "[", 1))
  compact_raw$strand <- sapply(strsplit(sapply(strsplit(as.character(compact_raw$circRNA),'_'), "[", 2), '!'), "[", 1)
  compact_raw$gene_symbol <- sapply(strsplit(sapply(strsplit(as.character(compact_raw$circRNA),'_'), "[", 2), '!'), "[", 2)
  
  expression <- compact_raw[,c("chr", "start", "stop", "strand", "gene_symbol", "V1")]
  colnames(expression)[6] <- sample
  
  if(is.null(finaldata)){
    finaldata <- expression
  } else {
    finaldata <- merge(finaldata, expression, by = c("chr", "start", "stop", "strand", "gene_symbol"), all = T)
  }
}
finaldata[is.na(finaldata)] <- 0
write.table(finaldata, paste0("circRNA_counts_raw.tsv"), quote = F, sep = "\t", row.names = F)
