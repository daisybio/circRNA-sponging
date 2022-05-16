#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(BSgenome, Biostrings, GenomicRanges, GenomicFeatures, seqinr)


args = commandArgs(trailingOnly = TRUE)

genome.version <- args[1]
genomes <- available.genomes(splitNameParts=TRUE)
select.genome <- genomes$pkgname[genomes$genome == genome.version]
if (length(select.genome) == 0) {
  stop("Error: Genome version not found in BSgenomes")
}
# set user library
dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
.libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
# select default genome
select.genome <- select.genome[!grepl("masked", select.genome)]
BiocManager::install(select.genome, ask = F)
library(package = select.genome, character.only = T)
fasta <- getBSgenome(select.genome)
# load circ expression
circ.expression <- read.table(args[2], sep = "\t", header = T)
# build Granges
ci.expression <- circ.expression[circ.expression$type == "ciRNA",]
circ.expression <- circ.expression[circ.expression$type == "circRNA",]
circ.ranges <- makeGRangesFromDataFrame(circ.expression[,1:4])
ci.ranges <- makeGRangesFromDataFrame(ci.expression[,1:4])
# set names according to annotation
if ("circBaseID" %in% colnames(circ.expression)) {
  circ_RNA_annotation <- ifelse(circ.expression$circBaseID != "None", 
                                circ.expression$circBaseID, 
                                paste0(circ.expression$chr, ":", circ.expression$start, ":", circ.expression$stop, ":", circ.expression$strand))
  ci_annotaion <- ifelse(ci.expression$circBaseID != "None", 
                         ci.expression$circBaseID, 
                         paste0(ci.expression$chr, ":", ci.expression$start, ":", ci.expression$stop, ":", ci.expression$strand))
  
} else {
  circ_RNA_annotation <- paste0(circ.expression$chr, ":", circ.expression$start, "-", circ.expression$stop, "_", circ.expression$strand)
  ci_annotaion <- paste0(ci.expression$chr, ":", ci.expression$start, "-", ci.expression$stop, "_", ci.expression$strand)
}
names(circ.ranges) <- circ_RNA_annotation
names(ci.ranges) <- ci_annotaion

# load gtf
gtf.pkg <- BiocManager::available(paste0(genome.version, ".knownGene"))
BiocManager::install(gtf.pkg, ask = F)
library(package = gtf.pkg, character.only = T)
gtf <- eval(parse(text=paste0(gtf.pkg,"::",gtf.pkg)))
exons <- exonicParts(gtf)

# splice sequences according to circ annotation
circ.exons <- findOverlapPairs(circ.ranges, exons, type = "any")
ids <- names(circ.exons@first)
circ.exons <- circ.exons@second
names(circ.exons) <- ids
circ.sequences <- getSeq(fasta, circ.exons)
# get intronic sequences
ci.sequences <- getSeq(fasta, ci.ranges)
# combine sequences
all.sequences <- c(circ.sequences, ci.sequences)
# concatenate exons from same ids
all.sequences.frame <- as.data.frame(all.sequences)
all.sequences.frame$id <- names(all.sequences)
all.sequences.frame <- aggregate(. ~ id, all.sequences.frame, function(x) paste0(x, collapse = ""))
all.sequences.combined <- DNAStringSet(all.sequences.frame$x)
names(all.sequences.combined) <- all.sequences.frame$id
# write sequences to file
writeXStringSet(all.sequences.combined, filepath = "circRNAs.fa")
