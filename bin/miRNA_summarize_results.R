#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
output_dir = args[2]

library(ggplot2)

dataset <- read.table(dataset_path, sep = "\t", header=T, stringsAsFactors = F)
samples <- dataset$sample

dataset_counts_raw <- NULL
for(i in 1:length(samples)){
  # get results
  sample <- samples[i]
  sample_folder <- paste0(output_dir, "/samples/",  sample, "/miRNA_detection/")
  miRDeep2_output <- read.table(dir(sample_folder, full.names=T, pattern="miRNAs_expressed"), sep = "\t", header=F, stringsAsFactors = F)
  names(miRDeep2_output) <- c("miRNA", "counts", "precursor", "total", "seq", "norm")
  
  raw_counts <- miRDeep2_output[,c(1,2)]
  
  # sum up counts which come from the same miRNA because of different precursors
  raw_counts <- aggregate(raw_counts$counts, by=list(miRNA=raw_counts$miRNA), FUN=sum)
  names(raw_counts) <- c("miRNA", sample)
  if(is.null(dataset_counts_raw)){
    dataset_counts_raw <- raw_counts
  } else {
    dataset_counts_raw <- merge(dataset_counts_raw, raw_counts, all = T, by = "miRNA")    
  }
  
}
dataset_counts_raw[is.na(dataset_counts_raw)]  <- 0
write.table(dataset_counts_raw, paste0("miRNA_counts_raw.tsv"), quote = F, sep = "\t", row.names = F)
