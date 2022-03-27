#!/usr/bin/env Rscript

library(Biostrings)
library(argparser)

args <- commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
parser <- add_argument(parser, "--circ_counts", help = "circRNA counts file")
parser <- add_argument(parser, "--dir", help = "Directory containing all samples with calculated abundances", default = ".")
parser <- add_argument(parser, "--keep_tmp", help = "Keep temporary files", default = T)

argv <- parse_args(parser, argv = args)


# read circRNA counts
circ.counts <- read.table(argv$circ_counts, header = T, sep = "\t", stringsAsFactors = F, check.names = F)

# determine annotation
annotation <- "circBaseID" %in% colnames(circ.counts)
# set search key
key <- ifelse(annotation, "circ", "chr")
# set row names
if (annotation) {
  rownames(circ.counts) <- circ.counts$circBaseID
} else {
  circ.counts$key <- paste0(circ.counts$chr, ":", circ.counts$start, "-", circ.counts$stop, "_", circ.counts$strand)
  circ.counts <- circ.counts[!duplicated(circ.counts$key),]
  circ.counts$key <- NULL
  rownames(circ.counts) <- paste0(circ.counts$chr, ":", circ.counts$start, "-", circ.counts$stop, "_", circ.counts$strand)
}

circ.quant <- circ.counts
mRNA.quant <- data.frame()

abundances <- list.files(path = argv$dir, pattern = "abundance.tsv", recursive = T, full.names = T)

for (path in abundances) {
  # get sample name
  sample <- basename(dirname(path))
  # read created abundances
  abundance <- read.table(normalizePath(path), sep = "\t", header = T)
  # split abundances accoording to circular and linear
  abundance$RNAtype <- ifelse(grepl(key, abundance$target_id), "circular", "linear")
  circ_or_linear <- split(abundance, abundance$RNAtype)
  abundance.circ <- circ_or_linear$circular
  abundance.mRNA <- circ_or_linear$linear
  
  # set counts to quantified levels
  circ.quant[abundance.circ$target_id, sample] <- abundance.circ$est_counts
  # save linear transcripts
  mRNA.quant[abundance.mRNA$target_id, sample] <- abundance.mRNA$est_counts
}

library(biomaRt)
# convert transcript ids to gene ids
mart <- NULL
not_done <- T
while(not_done){
  tryCatch(
    {
      mart <- useEnsembl(biomart = "ensembl", "hsapiens_gene_ensembl")
    },
    error=function(cond) {
      message(cond)
    },
    warning=function(cond) {
      message(cond)
    },
    finally={
      not_done = F
    }
  )
}

# convert row names
transcript.IDs <- sapply(strsplit(row.names(mRNA.quant), "\\."), "[", 1)
rownames(mRNA.quant) <- transcript.IDs
transcript2gene <- getBM(attributes=c("ensembl_transcript_id","external_gene_name","ensembl_gene_id"),
                         filters = "ensembl_transcript_id",
                         values = transcript.IDs, 
                         mart = mart)
conv <- merge(mRNA.quant, transcript2gene, by.x = 0, by.y = 1)
conv$Row.names <- NULL
conv$external_gene_name <- NULL
conv <- aggregate(conv[,-ncol(conv)], list(Gene=conv$ensembl_gene_id), FUN = sum)
row.names(conv) <- conv$Gene
conv$Gene <- NULL
# save converted samples
mRNA.quant <- conv
# remove tmp
if (!argv$keep_tmp) {
  unlink("tmp", recursive = T)
}
# write output to disk
linear.o <- "quant_linear_expression.tsv"
cat("writing mRNA output to", linear.o, "\n")
write.table(mRNA.quant, file = linear.o, sep = "\t", row.names = T)
o <- "quant_circ_expression.tsv"
cat("writing output file to ", o, "...\n")
write.table(circ.quant, file = o, sep = "\t", row.names = T)
print("done")

# extract counts
# circ.counts.c <- circ.counts[order(rownames(circ.counts)), -c(1:8)]
# circ.quant.c <- circ.quant[order(rownames(circ.quant)), -c(1:8)]
# calculate false positive rate
# tmp <- circ.counts.c - circ.quant.c
# fp <- length(tmp[tmp ==circ.counts.c])
# fp.rate <- fp/prod(dim(circ.counts.c))
# cat("False positive rate:", fp.rate, "according to psirc quantification\n")
