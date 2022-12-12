#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Two argument must be supplied", call.=FALSE)
}
dataset_path = args[1]
input_dir = args[2]
output_dir = "./"

dataset <- read.table(dataset_path, sep = "\t", header = T, stringsAsFactors = F)
samples <- dataset$sample

finaldata <- NULL
N_of_circRNAs_raw <- c()
for (i in 1:length(samples)){
  sample <- samples[i]
  path <- file.path(input_dir,"samples", sample, "circRNA_detection", "circExplorer2", paste0(sample ,"_circularRNA_known.txt"))
  message("Processing CIRCexplorer2 output for sample: ", sample)
  # skip if no circRNAs found
  if(file.size(path) == 0){
    warning(paste0("Sample ", sample, " has no circRNAs that could be detected with CIRCexplorer2. Sample will be removed from results"))
    next
  }
  file <- read.table(path, sep = "\t")
  CIRCexplorer2_output <- data.table(file)
  
  N_of_circRNAs_raw[i] <- nrow(CIRCexplorer2_output)
  names(N_of_circRNAs_raw)[i] <- sample
  
  expression_raw <- CIRCexplorer2_output[,c(1,2,3,6,13,15,16,14)]
  colnames(expression_raw) <- c("chr", "start", "stop", "strand", "counts", "gene_symbol", "isoform", "type")
  expression_raw <- expression_raw[,c(1,2,3,4,5,6,8)]
  
  # REMOVE SCAFFOLDS/ CONTIG CHROMOSME ENTRIES
  expression_raw <- expression_raw[!grepl("_", expression_raw$chr, fixed = T),]
  
  # compact and remove duplicates
  compact_raw <- data.table(circRNA=paste(expression_raw$chr,
                                           expression_raw$start,
                                           expression_raw$stop,
                                           expression_raw$strand,
                                           expression_raw$gene_symbol,
                                           expression_raw$type, sep = ":"),
                            counts = expression_raw$counts)
  compact_raw <- compact_raw[, max(counts), by=circRNA]
  # rebuild original data
  data <- strsplit(as.character(compact_raw$circRNA),':')
  compact_raw$chr <- sapply(data, "[", 1)
  compact_raw$start <- sapply(data, "[", 2)
  compact_raw$stop <- sapply(data, "[", 3)
  compact_raw$strand <- sapply(data, "[", 4)
  compact_raw$gene_symbol <- sapply(data, "[", 5)
  compact_raw$type <- sapply(data, "[", 6)
  # build ID
  rownames(compact_raw) <- paste0(compact_raw$chr, ":", compact_raw$start, "-", compact_raw$stop, "_", compact_raw$strand)
  
  expression <- compact_raw[,c("chr", "start", "stop", "strand", "gene_symbol", "type", "V1")]
  colnames(expression)[7] <- sample
  # first entry
  if(is.null(finaldata)){
    finaldata <- expression
  } else {
    finaldata <- merge(finaldata, expression, by = c("chr", "start", "stop", "strand", "gene_symbol", "type"), all = T)
  }
}
# remove NAs
finaldata[is.na(finaldata)] <- 0
# write to disk
write.table(finaldata, file.path(output_dir, "circRNA_counts_raw.tsv"), quote = F)
