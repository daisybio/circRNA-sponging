#!/usr/bin/env Rscript

install.packages("pacman")
pacman::p_load(BSgenome, AnnotationHub, argparser, ggplot2, biomaRt, ensembldb, Biostrings)


ah <- AnnotationHub(ask = F)

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
# ARGS
parser <- add_argument(parser, "--circ_targets", help = "circRNA targets as table in tsv format")
parser <- add_argument(parser, "--circ_fasta", help = "circRNA fasta file")
parser <- add_argument(parser, "--linear_targets", help = "mRNA 3UTR targets file as table in tsv format")
parser <- add_argument(parser, "--miRWalk_data", help = "mRNA data directly from miRWalk")
parser <- add_argument(parser, "--type", help = "Type of linear RNAs: CDS, 3UTR or 5UTR")
parser <- add_argument(parser, "--organism", help = "Organsim in three letter code e.g. hsa for Human")
parser <- add_argument(parser, "--biomaRt", help = "BiomaRt dataset name")

argv <- parse_args(parser, argv = args)

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
# load biomaRt
mart <- safe.mart.ensembl(argv$biomaRt)
# load targets
print("reading circRNA targets")
circ.targets <- read.table(argv$circ_targets, sep = "\t", header = T)
# calculate total n of targets for circRNAs over all samples
circ.targets <- rowSums(circ.targets)
print("reading circRNA fasta")
circ.fasta <- readDNAStringSet(argv$circ_fasta)
circ.fasta <- data.frame(names(circ.fasta), circ.fasta@ranges@width, row.names = 1)

circ <- data.frame(row.names = 1, merge(circ.targets, circ.fasta, by = 0))
colnames(circ) <- c("bindsites", "length")

print("reading miRWalk2.0 targets")
linear.targets <- fread(argv$linear_targets, select = 1:2, sep = "\t", header = T)
# build table
linear.targets <- table(linear.targets$mRNA, linear.targets$miRNA)
# convert genebank names to ENSG
conv <- data.frame(
    getBM(attributes = c("refseq_mrna", "ensembl_gene_id"),
          filters = "refseq_mrna", values = rownames(linear.targets), mart = mart)
    , check.names = F, stringsAsFactors = F, row.names = 1)
rownames(linear.targets) <- conv[rownames(linear.targets),]
# aggregate duplicate rows
linear.targets <- as.data.frame.matrix(linear.targets) %>% group_by(rownames(linear.targets)) %>% summarise_all(funs(sum))
linear.targets <- data.frame(linear.targets[!is.na(linear.targets$`rownames(linear.targets)`),], row.names = 1, check.names = F, stringsAsFactors = F)
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
mRNA$bindsites <- linear.targets[rownames(mRNA)]
colnames(mRNA)[1] <- c("length")

# plot ratios
p <- ggplot() +
    geom_point(data = mRNA, aes(x = length, y = bindsites, colour = "mRNA"), shape = 4) +
    geom_smooth(data = mRNA, aes(x = length, y = bindsites, col = "mRNA regression"), method = "lm", col = "darkorchid3") +
    geom_point(data = circ, aes(x = length, y = bindsites, colour = "circRNA"), shape = 2) +
    geom_smooth(data = circ, aes(x = length, y = bindsites, col = "circRNA regression"), method = "lm", col = "darkorange3") +
    scale_y_log10() +
    scale_x_log10() +
    labs(x = "length (log10)", y = "bindsites (log10)", colour = "Legend", title = argv$type) +
    scale_colour_manual(values = c(mRNA="#0066CC", circRNA="#006033", "mRNA regression"="darkorchid3", "circRNA regression"="darkorange3")) +
    theme(text = element_text(size=20)) +
    guides(color = guide_legend(override.aes = list(size = 3, shape = c(16,16,16,16))))
# save as png
png(paste0(argv$type, ".png"), res = 200, height = 800, width = 1300)
plot(p)
dev.off()

# different plot
mRNA$type <- argv$type
circ$type <- "circRNA"
data <- rbind(mRNA, circ)
# plot binding sites with linear regression
bindings <- ggplot(data, aes(x=length, y=bindsites, col=type, shape = factor(type))) +
    geom_point(size = 1.75) +
    geom_smooth(method = "lm", se = F) +
    scale_x_log10() + scale_y_log10() +
    scale_color_manual(values = c("3UTR"="#0066CC", circRNA="#006033")) +
    scale_shape_manual(values = c(6,1))
