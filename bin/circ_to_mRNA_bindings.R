#!/usr/bin/env Rscript

library(BSgenome)
library(AnnotationHub)
library(argparser)
library(ggplot2)
library(biomaRt)

ah <- AnnotationHub(ask = F)

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for circRNA quantification", name = "quant_parser")
# ARGS
parser <- add_argument(parser, "--miranda", help = "Miranda targets file from pipeline")
parser <- add_argument(parser, "--tarpmir", help = "Tarpmir targets file as from pipeline but only pairs")
parser <- add_argument(parser, "--mRNA3UTR", help = "mRNA 3UTR targets file as table in tsv format")
parser <- add_argument(parser, "--mRNA5UTR", help = "mRNA 5UTR targets file as table in tsv format")
parser <- add_argument(parser, "--CDS", help = "mRNA CDS targets file as table in tsv format")
parser <- add_argument(parser, "--circ_annotation", help = "circRNA annotation file")
parser <- add_argument(parser, "--organism", help = "Organsim in three letter code e.g. hsa for Human")

argv <- parse_args(parser, argv = args)

gen_counts <- function(sums, mRNA.all){
  # merge with bindsites counts
  mirWalk.all <- merge(sums, mRNA.all, by.x = 0, by.y = "gene_id")
  mirWalk.all <- mirWalk.all[, c("Row.names", colnames(mirWalk.all)[2], "width.x")]
  colnames(mirWalk.all) <- c("gene_id", "bindsites", "length")
  tmp <- aggregate(length~gene_id, data=mirWalk.all, sum)
  tmp$bindsites <- mirWalk.all[!duplicated(mirWalk.all$gene_id), "bindsites"]
  mirWalk.all <- tmp
  mirWalk.all <- mirWalk.all[mirWalk.all$length != 0,]
  mirWalk.all$score <- mirWalk.all$bindsites / mirWalk.all$length
  return(mirWalk.all)
}

plot_bindsites <- function(mRNA, circ, name){
  # plot bindsites to length for mRNA3UTR and circRNA
  bindsites.to.length.plot <- ggplot() + 
    geom_point(data = mRNA, aes(x = length, y = bindsites, colour = "mRNA"), shape = 4) +
    geom_point(data = circ, aes(x = spliced.length, y = bindsites, colour = "circRNA"), shape = 2) + scale_y_log10()+     
    scale_x_log10() + 
    labs(x = "length (log10)", y = "bindsites (log10)", color = "Legend", title = name) +
    scale_colour_manual("", 
                        breaks = c("mRNA", "circRNA"),
                        values = c("#0066CC", "#006033")) +
    theme(text = element_text(size=20))
  png(paste(name, "png", sep = "."))
  plot(bindsites.to.length.plot)
  dev.off()
}

# load targets
print("reading miranda targets")
miranda.targets <- read.table(argv$miranda, sep = "\t", header = T)
miranda.targets <- as.data.frame.matrix(table(miranda.targets$Target, miranda.targets$miRNA))
miranda.scores <- rowSums(miranda.targets)
print("reading tarpmir targets")
tarpmir.targets <- read.table(argv$tarpmir, sep = "\t", header = F)
tarpmir.targets <- as.data.frame.matrix(table(tarpmir.targets$V2, tarpmir.targets$V1))
tarpmir.scores <- rowSums(tarpmir.targets)
# combine scores
circ.targets <- merge(miranda.scores, tarpmir.scores, by = 0)
colnames(circ.targets) <- c("circRNA.ID", "miranda_targets", "tarpmir_targets")

# TODO: save converted counts to disk
print("reading mirWalk targets")
mirWalk.targets <- read.table(argv$mRNA3UTR, sep = "\t", header = T)
mRNA.3UTR.sums <- rowSums(mirWalk.targets)
mRNA.5UTR.sums.df <- read.table(argv$mRNA5UTR, sep = "\t", header = T)
mRNA.5UTR.sums.df <- as.data.frame.matrix(table(mRNA.5UTR.sums.df$mRNA, mRNA.5UTR.sums.df$miRNA))
mRNA.5UTR.sums.df <- data.frame(rowSums(mRNA.5UTR.sums.df))
mRNA.CDS.sums.df <- read.table(argv$CDS, sep = "\t", header = T)
mRNA.CDS.sums.df <- as.data.frame.matrix(table(mRNA.CDS.sums.df$mRNA, mRNA.CDS.sums.df$miRNA))
mRNA.CDS.sums.df <- data.frame(rowSums(mRNA.CDS.sums.df))

# convert genebank to ENSIDs
mart <- 0
not_done=TRUE
while(not_done)
{
  tryCatch({
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
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
gene.ens.all <- biomaRt::getBM(attributes = c("ensembl_gene_id", "refseq_mrna", "refseq_ncrna"),
                               filter = "refseq_mrna",
                               values=unique(c(names(mRNA.5UTR.sums.df), names(mRNA.CDS.sums.df))),
                               mart=mart,
                               checkFilters = T)
gene.ens.all <- gene.ens.all[!duplicated(gene.ens.all$ensembl_gene_id),]
gene.ens.all <- gene.ens.all[!is.null(gene.ens.all$refseq_mrna),]
# replace ids
matches5 <- merge(gene.ens.all, mRNA.5UTR.sums.df, by.x = "refseq_mrna", by.y = 0)
colnames(matches5)[5] <- "bindsites"
matches5 <- aggregate(matches5$bindsites, list(ensembl_gene_id=matches5$ensembl_gene_id), FUN=sum)
rownames(matches5) <- matches5$ensembl_gene_id
mRNA.5UTR.sums.df <- matches5[, 2, drop = F]
matchesCDS <- merge(gene.ens.all, mRNA.CDS.sums.df, by.x = "refseq_mrna", by.y = 0)
colnames(matchesCDS)[5] <- "bindsites"
matchesCDS <- aggregate(matchesCDS$bindsites, list(ensembl_gene_id=matchesCDS$ensembl_gene_id), FUN=sum)
rownames(matchesCDS) <- matchesCDS$ensembl_gene_id
mRNA.CDS.sums.df <- matchesCDS[, 2, drop = F]

# load annotations
print("reading circ annotations")
circ.annotations <- read.table(argv$circ_annotation, sep = "\t", header = T)
print("loading GTF")
organism <- argv$organism

# load latest EnsDb release for given organism
gtf <- tail(query(ah, pattern = c(organism, "EnsDb")), 1)[[1]]
# load genes
genes <- genes(gtf)
genes.df <- as.data.frame(genes)
# extract 3UTR regions
mRNA.3UTR <- threeUTRsByTranscript(gtf)
mRNA.3UTR.df <- as.data.frame(mRNA.3UTR)
# extract 5UTR regions
mRNA.5UTR <- fiveUTRsByTranscript(gtf)
mRNA.5UTR.df <- as.data.frame(mRNA.5UTR)
# extract cds regions
mRNA.CDS.regions <- cdsBy(gtf)
mRNA.CDS.regions.df <- as.data.frame(mRNA.CDS.regions)
# get gene id for transcripts
mRNA.3UTR.all <- merge(mRNA.3UTR.df, genes.df, by.x = "group_name", by.y = "canonical_transcript")
mRNA.5UTR.all <- merge(mRNA.5UTR.df, genes.df, by.x = "group_name", by.y = "canonical_transcript")
mRNA.CDS.regions.all <- merge(mRNA.CDS.regions.df, genes.df, by.x = "group_name", by.y = "canonical_transcript")
# get counts
mRNA.3UTR.final <- gen_counts(mRNA.3UTR.sums, mRNA.3UTR.all)
mRNA.5UTR.final <- gen_counts(mRNA.5UTR.sums.df, mRNA.5UTR.all)
mRNA.CDS.final <- gen_counts(mRNA.CDS.sums.df, mRNA.CDS.regions.all)

# get scores for circRNAs
circ.annotations <- merge(circ.annotations, circ.targets)
# exclude tarpmir matches
circ.annotations$bindsites <- circ.annotations$miranda_targets + circ.annotations$tarpmir_targets
# circ.annotations$bindsites <- circ.annotations$miranda_targets
circ.annotations$bindscores <- circ.annotations$bindsites / circ.annotations$spliced.length
circ.score <- mean(circ.annotations$bindscores[!is.na(circ.annotations$bindscores)])
# get scores for mRNAs
mRNA.3UTR.score <- mean(mRNA.3UTR.final$score[!is.na(mRNA.3UTR.final$score)])
mRNA.5UTR.score <- mean(mRNA.5UTR.final$score[!is.na(mRNA.5UTR.final$score)])
mRNA.CDS.score <- mean(mRNA.CDS.final$score[!is.na(mRNA.CDS.final$score)])

# calc ratio
ratio <- circ.score / mRNA.3UTR.score
cat("mean circRNA bindings per sequence length:", circ.score, "\nmean mRNA bindings per sequence length", mRNA.3UTR.score, 
    "\ncircRNAbpl/mRNAbpl =", ratio, "\n")

# barplot
barplot <- ggplot(data = data.frame(name=c("circRNA", "3UTR", "5UTR", "CDS"), value = c(circ.score, mRNA.3UTR.score, mRNA.5UTR.score, mRNA.CDS.score)), aes(x=name, y=value)) + 
  geom_bar(stat = "identity", fill = c("#006633", "#003399", "#CC3333", "#CC6633")) + xlab("RNA type") + ylab("mean of bindsites/length") +
  theme(text = element_text(size=20))
png("circToLinearBindings.png")
plot(barplot)
dev.off()

# plot individual mRNA types
plot_bindsites(mRNA = mRNA.3UTR.final, circ = circ.annotations, name = "3UTR")
plot_bindsites(mRNA = mRNA.5UTR.final, circ = circ.annotations, name = "5UTR")
plot_bindsites(mRNA = mRNA.CDS.final, circ = circ.annotations, name = "CDS")
# plot all mRNA against circRNA
all <- rbind(mRNA.3UTR.final, mRNA.5UTR.final, mRNA.CDS.final)
plot_bindsites(mRNA = all, circ = circ.annotations, name = "ALL")

