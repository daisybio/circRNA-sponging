#!/usr/bin/env Rscript

library(SPONGE)

args = commandArgs(trailingOnly = TRUE)

# load SPONGE workspace
load(args[1])
# differentially expressed circRNAs
signif.hits <- read.table(args[2], sep = "\t", header = T)


# get all differentially expressed circRNAs that act as ceRNAs
differential.circ.in.ce.network.A <- merge(ceRNA_interactions_fdr, signif.hits, by.x = "geneA", by.y = "X")
differential.circ.in.ce.network.B <- merge(ceRNA_interactions_fdr, signif.hits, by.x = "geneB", by.y = "X")
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network.A, differential.circ.in.ce.network.B, all = T)

differential.circ.plot <- sponge_plot_network(differential.circ.in.ce.network, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
differential.circ.plot$x$edges$label <- paste("mscor:", round(differential.circ.in.ce.network$mscor, 2))
differential.circ.plot
