#!/usr/bin/env Rscript
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("At least two argument must be supplied", call.=FALSE)
}
expression_norm_path = args[1]
output_dir = args[2]
samples_percentage = 0.2 # default 0.2, minimum percentage of samples, a circRNA has to be expressed in is to pass filtering
read_cutoff = 5 # default 5, minimum number of reads, a circRNA is required to have to pass filtering

expression_norm <- read.table(expression_norm_path, sep = "\t", header=T, stringsAsFactors = F, check.names = F)

samples <- colnames(expression_norm)[-c(1:4)]

# filter data: counts > 5 in at least 20% of samples
if(length(samples) < 5){
  stop("Cannot perform filtering on less than 5 samples")
}
sample_nr_cutoff <- max(floor(samples_percentage *length(samples)),1)

rows_to_keep <- c()
for (i in 1:nrow(expression_norm)){
  number_of_samples_containing_this_circRNA <- 0
  for (j in 5:ncol(expression_norm)){
    if(expression_norm[i,j] >= read_cutoff){
      number_of_samples_containing_this_circRNA <- number_of_samples_containing_this_circRNA + 1
    }
  }
  if(number_of_samples_containing_this_circRNA >= sample_nr_cutoff){
    rows_to_keep <- append(rows_to_keep, i)
  }
  
}
filtered_data <- expression_norm[rows_to_keep,]
write.table(filtered_data, paste0("circRNA_count_filtered.tsv"), quote = F, sep = "\t", row.names = F)

