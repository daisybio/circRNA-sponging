#!/usr/bin/env Rscript
library("Rsubread", "DESeq2")

args = commandArgs(trailingOnly = TRUE)
# read sam file
sam_file <- args[1]
# genome version
genome_version <- args[2]
# gtf
gtf_file <- args[3]
# convert sam to countsFile
countsData <- featureCounts(sam_file, annot.inbuilt = genome_version, annot.ext = gtf_file,
                            isGTFAnnotationFile = TRUE)
# prepare data
dds <- DESeqDataSetFromMatrix(countData = countsData)
# differential expression analysis
dds <- DESeq(dds)
# save results
res <- results(dds)
