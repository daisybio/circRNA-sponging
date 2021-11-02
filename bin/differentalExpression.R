#!/usr/bin/env Rscript
library("Rsubread")

args = commandArgs(trailingOnly = TRUE)
args <- c("/Users/leonschwartz/Desktop/Bioinformatik/local_data/fastqs/test.sam,/Users/leonschwartz/Desktop/Bioinformatik/local_data/fastqs/end.sam",
          "hg38",
          "/Users/leonschwartz/Desktop/Bioinformatik/local_data/references/gencode.v35.primary_assembly.annotation.gtf",
          "TRUE")
# read sam file(s)
sam_files <- strsplit(args[1], ",")[[1]]
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

countsData <- Rsubread::featureCounts(sam_files, annot.inbuilt = genome_version, annot.ext = gtf_file,
                            isGTFAnnotationFile = TRUE, isPairedEnd = isPaired)
# write output file
write.table(countsData[["counts"]], file = "countsFile.tsv", quote = FALSE, sep = "\t", col.names = NA)
