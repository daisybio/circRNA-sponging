#!/usr/bin/env Rscript

library(SPONGE)

args = commandArgs(trailingOnly = TRUE)

# load SPONGE workspace
load(args[1])
# differentially expressed circRNAs
differential.circ <- read.table(args[2], sep = "\t", header = T)
# PARAMS
FDR <- 0.1
L2FC <- 1

# get significantly differentially expressed circRNAs
signif.hits <- differential.circ[!is.na(differential.circ$padj) &
                                   differential.circ$padj<as.double(FDR) &
                                   abs(differential.circ$log2FoldChange)>=as.numeric(L2FC),]

# get all differentially expressed circRNAs that act as ceRNAs
differential.circ.in.ce.network.A <- merge(ceRNA_interactions_fdr, signif.hits, by.x = "geneA", by.y = "ENS_ID")
differential.circ.in.ce.network.B <- merge(ceRNA_interactions_fdr, signif.hits, by.x = "geneB", by.y = "ENS_ID")
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network.A, differential.circ.in.ce.network.B, all = T)

differential.circ.plot <- sponge_plot_network(differential.circ.in.ce.network, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
differential.circ.plot$x$edges$label <- paste("mscor:", round(differential.circ.in.ce.network$mscor, 2))
differential.circ.plot
