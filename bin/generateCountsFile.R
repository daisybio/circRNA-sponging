#!/usr/bin/env Rscript
library("Rsubread")

args = commandArgs(trailingOnly = TRUE)
# read sam file
sam_file <- args[1]
# genome version
genome_version <- args[2]
# gtf
gtf_file <- args[3]
# read mode
parse_boolean <- function(isT) {
  if (isT == "TRUE") {
    return (TRUE)
  } else {
    return (FALSE)
  }
}
isPaired <- parse_boolean(args[4])
# convert sam to countsFile

countsData <- Rsubread::featureCounts(sam_file, annot.inbuilt = genome_version, annot.ext = gtf_file,
                            isGTFAnnotationFile = TRUE, isPairedEnd = isPaired)
# write output file
write.table(countsData, file = "countsFile.tsv", quote = FALSE, sep = "\t", col.names = )
