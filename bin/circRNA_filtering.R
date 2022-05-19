#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Six arguments must be supplied", call.=FALSE)
}

expression_norm_path = args[1]
samples <- read.table(file = args[2], sep = "\t", header = T)[,"sample"] # get samples from samplesheet
output_dir = args[3]
samples_percentage = as.numeric(args[4]) # default 0.2, minimum percentage of samples, a circRNA has to be expressed in is to pass filtering
read_cutoff = as.numeric(args[5]) # default 5, minimum number of reads, a circRNA is required to have to pass filtering

expression_norm <- read.table(expression_norm_path, sep = "\t", header=T, stringsAsFactors = F, check.names = F)

# filter data: counts > 5 in at least 20% of samples
if(length(samples) < 5){
  stop("Cannot perform filtering on less than 5 samples")
}
sample_nr_cutoff <- ceiling(samples_percentage *length(samples))
# select expressions only
expr_only <- expression_norm[,samples]
rows_to_keep <- c()
for (i in 1:nrow(expr_only)){
  number_of_samples_containing_this_circRNA <- 0
  for (j in 1:ncol(expr_only)){
    if(expr_only[i,j] >= read_cutoff){
      number_of_samples_containing_this_circRNA <- number_of_samples_containing_this_circRNA + 1
    }
  }
  if(number_of_samples_containing_this_circRNA >= sample_nr_cutoff){
    rows_to_keep <- append(rows_to_keep, i)
  }
}
filtered_data <- expression_norm[rows_to_keep,]
# write final output
write.table(filtered_data, file = "circRNA_counts_filtered.tsv", quote = F, sep = "\t", row.names = F)
