#!/usr/bin/env Rscript
# install required packages if they do not already exist
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "Rsubread"))

# start script
library("Rsubread", "DESeq2")

args = commandArgs(trailingOnly = TRUE)

# read sam file(s)
sam_files <- list.files(args[1], pattern = "\\.sam$", full.names = TRUE, recursive = TRUE)

# metadata
metaData <- read.table(file = args[2], sep = "\t", header = TRUE)
# genome version
genome_version <- args[3]
# gtf
gtf_file <- args[4]
# read mode
parse_boolean <- function(isT) {
  if (isT == "true") {
    return (TRUE)
  } else {
    return (FALSE)
  }
}
# check for paired reads
isPaired <- parse_boolean(args[5])
# create ouptut data and plots
create_outputs <- function(d, results, n, marker, file_name) {
  # write data to disk
  write.table(cbind(ENS_ID=rownames(results), results), file = paste(file_name, "tsv", sep = "."), quote = FALSE, sep = "\t", col.names = NA)
  # create normalized counts plots for top n genes
  if (length(results@rownames) >= n) {
    r_top = as.list(results@rownames)[1:n]
    for (gene in r_top) {
      pl = DESeq2::plotCounts(d, gene = gene, intgroup = marker)
      gene_name = paste(file_name, gene, sep = "_")
      png(filename = paste(gene_name, "png", sep = "."))
      plot(pl)
      dev.off
    }
  }
  # PCA
  vsdata <- DESeq2::vst(d, blind = FALSE)
  PCA_plot <- DESeq2::plotPCA(vsdata, intgroup = marker)
  plot_name <- paste(file_name, "pca", sep = "_")
  png(filename = paste(plot_name, "png", sep = "."))
  plot(PCA_plot)
  dev.off
}
# convert sam to countsObject
countsObject <- Rsubread::featureCounts(sam_files, annot.inbuilt = genome_version, annot.ext = gtf_file,
                            isGTFAnnotationFile = TRUE, isPairedEnd = isPaired)
# prepare deseq2 counts data
countsData <- countsObject[["counts"]]
# save general gene expression of all samples
write.table(countsData, file = paste("gene_expression", "tsv", sep = "."), quote = FALSE, sep = "\t", col.names = NA)
# reformat ids
new_names <- list()
for (name in row.names(countsData)) {
  general_name <- strsplit(name, "\\.")[[1]][1]
  new_names <- c(new_names, general_name)
}
row.names(countsData) <- new_names

# dds object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countsData,
                                      colData = metaData,
                                      design = ~condition)

# differential expression analysis
dds <- DESeq2::DESeq(dds)
# results
res <- DESeq2::results(dds)
# sort by p-value
res <- res[order(res$padj),]
# create summary
DESeq2::summary(res)

# WRITE OUTPUTS
# total_RNA
create_outputs(d = dds, results = res, n = 5, marker = "condition", file_name = "total_rna")
# only circRNA, if file loc is given
if (length(args) == 6) {
  circ_RNAs <- read.table(file = args[6], sep = "\t", header = TRUE)
  ens_ids <- c_RNAs$Ensembl_gene_ID
  # circ rnas only
  filtered_res <- res[row.names(res) %in% ens_ids,]
  create_outputs(d = dds, results = filtered_res, n = 5, marker = "condition", file_name = "circ_rna_only")
}