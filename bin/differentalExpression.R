#!/usr/bin/env Rscript
# install required packages if they do not already exist
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "Rsubread"))
# start script
library("Rsubread", "DESeq2")

args = commandArgs(trailingOnly = TRUE)
args <- c("/Users/leonschwartz/Desktop/Bioinformatik/local_data/pipeline_expample_output/output/samples/",
          "/Users/leonschwartz/Desktop/Bioinformatik/local_data/pipeline_expample_output/subset.tsv",
          "hg38",
          "/Users/leonschwartz/Desktop/Bioinformatik/local_data/references/gencode.v35.primary_assembly.annotation.gtf",
          "FALSE")
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
  if (isT == "TRUE") {
    return (TRUE)
  } else {
    return (FALSE)
  }
}
isPaired <- parse_boolean(args[5])
# convert sam to countsFile

countsObject <- Rsubread::featureCounts(sam_files, annot.inbuilt = genome_version, annot.ext = gtf_file,
                            isGTFAnnotationFile = TRUE, isPairedEnd = isPaired)
# prepare deseq2 counts data
countsData <- countsObject[["counts"]]

# write output file
# write.table(countsData, file = "countsFile.tsv", quote = FALSE, sep = "\t", col.names = NA)

# dds object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countsData,
                                      colData = metaData,
                                      design = ~condition,
                                      tidy = TRUE)

# differential expression analysis
dds <- DESeq(dds)
# results
res <- results(dds)
# sort by p-value
res <- res[order(res$padj),]
# write results to file
write.table(res, file = "results.tsv", quote = FALSE, sep = "\t", col.names = NA)
