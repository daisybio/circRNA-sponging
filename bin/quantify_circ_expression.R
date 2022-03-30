#!/usr/bin/env Rscript

library(Biostrings)
library(argparser)
library(DESeq2)
require(reshape2)

args <- commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
parser <- add_argument(parser, "--circ_counts", help = "circRNA counts file")
parser <- add_argument(parser, "--samplesheet", help = "Samplesheet of pipeline")
parser <- add_argument(parser, "--dir", help = "Directory containing all samples with calculated abundances", default = ".")

argv <- parse_args(parser, argv = args)

norm <- function(data, samples) {
  meta <- data.frame(samples)
  row.names(meta) <- meta$samples 
  data <- as.matrix(data)
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(data + 1), colData = meta, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  
  return(DESeq2::counts(dds, normalized=T))
}


# read circRNA counts
print("reading circRNA raw counts")
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

n <- length(abundances)
c <- 0
# save tpm counts -> log2(tpm+1)
circ.tpm <- data.frame()
mRNA.tpm <- data.frame()
for (path in abundances) {
  # get sample name
  sample <- basename(dirname(path))
  cat("processing sample", sample, "|", c/n, "\r")
  # read created abundances
  abundance <- read.table(normalizePath(path), sep = "\t", header = T)
  # split abundances accoording to circular and linear
  abundance$RNAtype <- ifelse(grepl(key, abundance$target_id), "circular", "linear")
  circ_or_linear <- split(abundance, abundance$RNAtype)
  abundance.circ <- circ_or_linear$circular
  abundance.mRNA <- circ_or_linear$linear
  
  # save tpms
  circ.tpm[abundance.circ$target_id, sample] <- log2(abundance.circ$tpm + 1)
  mRNA.tpm[abundance.mRNA$target_id, sample] <- log2(abundance.mRNA$tpm + 1)
  
  # save circular transcripts
  circ.quant[abundance.circ$target_id, sample] <- abundance.circ$est_counts
  # save linear transcripts
  mRNA.quant[abundance.mRNA$target_id, sample] <- abundance.mRNA$est_counts
  c = c + 1
}
colnames(circ.tpm) <- c("tpm")
colnames(mRNA.tpm) <- c("tpm")

print("Converting transcripts to genes...")
library(biomaRt)
# convert transcript ids to gene ids
mart <- NULL
not_done <- T
while(not_done){
  tryCatch(
    {
      mart <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
    },
    error=function(cond) {
      message(cond)
    },
    warning=function(cond) {
      message(cond)
    },
    finally={
      not_done = F
      print("mart sucessfully created")
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
# summarizing transcripts that map to same gene
print("Summarizing transcripts of same gene...")
conv <- aggregate(conv[,-ncol(conv)], list(Gene=conv$ensembl_gene_id), FUN = sum)
row.names(conv) <- conv$Gene
conv$Gene <- NULL
# save converted samples
mRNA.quant <- conv

# aggregate tpms
mRNA.tpm <- merge(mRNA.tpm, transcript2gene, by.x = 0, by.y = 1)
mRNA.tpm$Row.names <- NULL
mRNA.tpm$external_gene_name <- NULL
mRNA.tpm <- aggregate(mRNA.tpm[,-ncol(mRNA.tpm)], list(Gene=mRNA.tpm$ensembl_gene_id), FUN = sum)
row.names(mRNA.tpm) <- mRNA.tpm$Gene
mRNA.tpm$Gene <- NULL
mRNA.tpm <- mRNA.tpm

# bind and save tpms
tpms <- rbind(circ.tpm, mRNA.tpm)
write.table(tpms, file = "TPM_map.tsv", sep = "\t", row.names = T)

# write output to disk
linear.o <- "quant_linear_expression.tsv"
cat("writing mRNA output to", linear.o, "\n")
write.table(mRNA.quant, file = linear.o, sep = "\t", row.names = T)
o <- "quant_circ_expression.tsv"
cat("writing output file to ", o, "...\n")
write.table(circ.quant, file = o, sep = "\t", row.names = T)
print("done")

print("Calculating quantification effects...")

# plot estimated counts vs junction reads
# est counts <- quantification
# junction reads <- CIRCexplorer2
samplesheet <- read.table(argv$samplesheet, sep = "\t", header = T)
samples <- samplesheet$sample

circ.norm <- circ.counts
rownames(circ.norm) <- paste(circ.norm$chr, circ.norm$start, circ.norm$stop, circ.norm$strand, circ.norm$gene_symbol, sep = "_")
circ.norm.quant <- circ.quant
rownames(circ.norm.quant) <- paste(circ.norm.quant$chr, circ.norm.quant$start, circ.norm.quant$stop, circ.norm.quant$strand, circ.norm.quant$gene_symbol, sep = "_")
# normalize
circ.norm <- norm(circ.norm[,samples], samples)
circ.norm.quant <- norm(circ.norm.quant[,samples], samples)

# compare counts
circ.norm.sums <- as.data.frame(rowSums(circ.norm))
circ.norm.quant.sums <- as.data.frame(rowSums(circ.norm.quant))

sums <- merge(circ.norm.sums, circ.norm.quant.sums, by = 0)
colnames(sums) <- c("id", "count.x", "count.y")
sums$key <- rownames(sums)
sums$chr <- sapply(gsub("chr", "", sapply(strsplit(sums$id, "_"), "[", 1)), "[", 1)
chrOrder <-c((1:22),"X","Y","M")
sums <- sums[order(factor(sums$chr, levels = chrOrder, ordered = T)),]
chromosomes <- unique(sums$chr)
middle.chromosome.pos <- table(sums$chr)/2

psirc.color <- "#3333FF"
circExplorer2.color <- "#009933"

# plot both
png("quant_effects.png", width = 1200, height = 800, units = "px")
plot(sums$key, log10(sums$count.y), type = "b", col = psirc.color, xlab = "chromosome", ylab = "circRNA counts over all samples (log10)", xaxt="n")
axis(1, at = round(seq(1, nrow(sums), nrow(sums)/length(chromosomes)))+middle.chromosome.pos, labels = F)
text(round(seq(1, nrow(sums), nrow(sums)/length(chromosomes))), par("usr")[3] - 0.2, labels = chromosomes, srt = 45, pos = 1, xpd = TRUE)
lines(sums$key, log10(sums$count.x), type = "b", lty = 2, col = circExplorer2.color, pch = 18)
legend("top", inset = c(-0.5, 0), legend = c("CircExplorer2", "psirc-quant"), col = c(circExplorer2.color, psirc.color), lty = c(1,1), cex=0.8)
dev.off()

# extract counts
# circ.counts.c <- circ.counts[order(rownames(circ.counts)), -c(1:8)]
# circ.quant.c <- circ.quant[order(rownames(circ.quant)), -c(1:8)]
# calculate false positive rate
# tmp <- circ.counts.c - circ.quant.c
# fp <- length(tmp[tmp ==circ.counts.c])
# fp.rate <- fp/prod(dim(circ.counts.c))
# cat("False positive rate:", fp.rate, "according to psirc quantification\n")
