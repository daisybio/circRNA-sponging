#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Five argument must be supplied", call.=FALSE)
}
expression_norm_path = args[1]
output_dir = args[2]
samples_percentage = as.numeric(args[3]) # default 0.2, minimum percentage of samples, a miRNA has to be expressed in is to pass filtering
read_cutoff = as.numeric(args[4]) # default 5, minimum number of reads, a  miRNA is required to have to pass filtering
count_mode = args[5] # either use tpm or regular counts

expression_norm <- read.table(expression_norm_path, sep = "\t", header=T, stringsAsFactors = F, check.names = F)

samples <- colnames(expression_norm)[-c(1)]

# filter data: counts > 5 in at least 20% of samples
if(length(samples) < 5){
  stop("Cannot perform filtering on less than 5 samples")
}
sample_nr_cutoff <- ceiling(samples_percentage *length(samples))

rows_to_keep <- c()
for (i in 1:nrow(expression_norm)){
  number_of_samples_containing_this_miRNA <- 0
  for (j in 5:ncol(expression_norm)){
    if(expression_norm[i,j] >= read_cutoff){
      number_of_samples_containing_this_miRNA <- number_of_samples_containing_this_miRNA + 1
    }
  }
  if(number_of_samples_containing_this_miRNA >= sample_nr_cutoff){
    rows_to_keep <- append(rows_to_keep, i)
  }
  
}
filtered_data <- expression_norm[rows_to_keep,]

# convert counts to tpm
if (count_mode == "tpm"){
  print("converting counts to tpm...")
  len <- nrow(filtered_data)
  filtered_data <- filtered_data/len
  filtered_data <- t(t(filtered_data)*1e6/colSums(filtered_data))
}
write.table(filtered_data, paste0("miRNA_counts_filtered.tsv"), quote = F, sep = "\t", row.names = F)

