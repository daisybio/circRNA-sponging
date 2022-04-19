#!/usr/bin/env Rscript

library(BSgenome)
library(AnnotationHub)
library(argparser)
library(ggplot2)
library(biomaRt)
library(ensembldb)
library(Biostrings)

ah <- AnnotationHub(ask = F)

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
# ARGS
parser <- add_argument(parser, "--circ_targets", help = "circRNA targets as table in tsv format")
parser <- add_argument(parser, "--circ_fasta", help = "circRNA fasta file")
parser <- add_argument(parser, "--linear_targets", help = "mRNA 3UTR targets file as table in tsv format")
parser <- add_argument(parser, "--type", help = "Type of linear RNAs: CDS, 3UTR or 5UTR")
parser <- add_argument(parser, "--organism", help = "Organsim in three letter code e.g. hsa for Human")

argv <- parse_args(parser, argv = args)

plot_bindsites <- function(mRNA, circ, name){
  bindsites.to.length.plot <- ggplot() + 
    geom_point(data = mRNA, aes(x = length, y = bindsites, colour = "mRNA"), shape = 4) +
    geom_smooth(data = mRNA, aes(x = length, y = bindsites, colour = "mRNAregression"), method = "lm", col = "darkorchid3") +
    geom_point(data = circ, aes(x = length, y = bindsites, colour = "circRNA"), shape = 2) +
    geom_smooth(data = circ, aes(x = length, y = bindsites, colour = "circRNAregression"), method = "lm", col = "darkorange3") +
    scale_y_log10() + 
    scale_x_log10() + 
    labs(x = "length (log10)", y = "bindsites (log10)", colour = "Legend", title = name) +
    scale_colour_manual("",
                        breaks = c("mRNA", "circRNA", "mRNAregression", "circRNAregression"),
                        values = c("#0066CC", "#006033", "darkorchid3", "darkorange3")) +
    theme(text = element_text(size=20))
  png(paste(name, "png", sep = "."))
  plot(bindsites.to.length.plot)
  dev.off()
}

safe.mart.ensembl <- function(data.set){
  mart <- 0
  not_done=TRUE
  while(not_done)
  {
    tryCatch({
      mart <- useEnsembl(biomart = "ensembl", dataset = data.set)
      not_done=FALSE
    }, warning = function(w) {
      print("WARNING SECTION")
      print(w)
    }, error = function(e) {
      print("ERROR SECTION")
      print(e)
    }, finally = {
    })
  }
  print("SUCCESS")
  return(mart)
}

# load targets
print("reading circRNA targets")
circ.targets <- read.table(argv$circ_targets, sep = "\t", header = T)
# calculate total n of targets for circRNAs over all samples
circ.targets <- rowSums(circ.targets)
print("reading circRNA fasta")
circ.fasta <- readDNAStringSet(argv$circ_fasta)
circ.fasta <- data.frame(names(circ.fasta), circ.fasta@ranges@width, row.names = 1)

circ.all <- data.frame(row.names = 1, merge(circ.targets, circ.fasta, by = 0))
colnames(circ.all) <- c("bindsites", "length")

print("reading linear targets")
linear.targets <- read.table(argv$linear_targets, sep = "\t", header = T)
linear.targets <- rowSums(linear.targets)

organism <- argv$organism

# load latest EnsDb release for given organism
gtf <- tail(query(ah, pattern = c(organism, "EnsDb")), 1)[[1]]
# load genes
genes <- genes(gtf)
genes <- as.data.frame(genes)
# extract 3UTR regions
if (argv$type == "3UTR"){
  mRNA <- threeUTRsByTranscript(gtf)
} else if (argv$type == "5UTR"){
  mRNA <- fiveUTRsByTranscript(gtf)
} else if (argv$type == "CDS") {
  mRNA <- cdsBy(gtf)
} else {
  stop("invalid type of linear RNA given, use one of CDS, 3UTR or 5UTR")
}
mRNA <- as.data.frame(mRNA)
mRNA <- merge(mRNA, genes, by.x = "group_name", by.y = "canonical_transcript")
mRNA <- data.frame(row.names = 1, mRNA[!duplicated(mRNA$gene_id),c("gene_id", "width.x")])
mRNA <- data.frame(row.names = 1, colnames = c("bindsites", "length"), merge(linear.targets, mRNA, by = 0))
plot_bindsites(mRNA.3UTR.final, circ.all, argv$type)
