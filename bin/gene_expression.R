#!/usr/bin/env Rscript

library(Rsubread)

args = commandArgs(trailingOnly = TRUE)

# sam file
sam_file <- args[1]
# id 
id <- args[2]

# genome version
genome_version <- args[3]
# gtf
gtf_file <- args[4]
# read mode
is_single_end <- function(isT) {
  if (isT == "true") {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
# check for paired reads
isPaired <- is_single_end(args[5])

# convert sam to countsObject
countsObject <- Rsubread::featureCounts(sam_file, annot.inbuilt = genome_version, annot.ext = gtf_file,
                                        isGTFAnnotationFile = TRUE, isPairedEnd = isPaired)
# prepare deseq2 counts data
countsData <- data.frame(countsObject[["counts"]])
# sample id
id_c <- countsObject$targets
# reformat ids
countsData$ensembl_gene_id <- sapply(strsplit(as.character(rownames(countsData)), "\\."), "[", 1)
countsData <- aggregate(countsData[id_c], by=countsData['ensembl_gene_id'], sum)
rownames(countsData) <- countsData$ensembl_gene_id
countsData$ensembl_gene_id <- NULL
colnames(countsData) <- id

# save general gene expression of all samples
write.table(countsData, file = paste("gene_expression", "tsv", sep = "."), quote = FALSE, sep = "\t", col.names = NA)
