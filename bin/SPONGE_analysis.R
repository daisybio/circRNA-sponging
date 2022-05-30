#!/usr/bin/env Rscript

library(SPONGE)
library(visNetwork)
library(MetBrewer)

args = commandArgs(trailingOnly = TRUE)

# load SPONGE workspace
load(args[1])
# differentially expressed circRNAs
signif.hits.circ <- read.table(args[2], sep = "\t", header = T)
signif.hits.mRNA <- read.table(args[3], sep = "\t", header = T)
gtf <- args[4]

signif.hits <- rbind(signif.hits.circ, signif.hits.mRNA)

differential.circ.in.ce.network <- circ.mRNA.only
# convert ENSG to hgnc
gtf <- rtracklayer::readGFF(gtf)
gene.ens.all <- unique(gtf[!is.na(gtf$transcript_id),c("gene_id", "gene_name", "gene_biotype")])
colnames(gene.ens.all) <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
rownames(gene.ens.all) <- gene.ens.all$gene_id

differential.circ.plot <- sponge_plot_network(differential.circ.in.ce.network, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
differential.circ.plot$x$edges$label <- paste("mscor:", round(differential.circ.in.ce.network$mscor, 2))

hgncs <- merge(signif.hits, gene.ens.all, by.x = "X", by.y = "ensembl_gene_id")
nodes <- differential.circ.plot$x$nodes
nodes <- merge(nodes, gene.ens.all, by = 1, all.x = T)
# change to hgnc
nodes[!is.na(nodes$hgnc_symbol),"id"] <- nodes[!is.na(nodes$hgnc_symbol),"hgnc_symbol"]
# add circRNA as biotype
nodes[is.na(nodes$gene_biotype),"gene_biotype"] <- "circRNA"
biotypes <- unique(nodes$gene_biotype)
style <- data.frame(groupname=biotypes, 
                    shape=seq(1,length(biotypes)), 
                    color = as.vector(met.brewer("Juarez", length(biotypes))))
nodes$color <- NULL
nodes[nodes$id %in% hgncs$hgnc_symbol | nodes$id %in% signif.hits$X,"color"] <- "#CC3333"
nodes$shape <- ifelse(grepl("c", nodes$id), "rectangle", "triangle")
nodes$group <- ifelse(grepl("c", nodes$id), "circRNA", "mRNA")
nodes$group[grep("c", nodes$id)] <- "circRNA"
nodes[nodes$id %in% hgncs$hgnc_symbol | nodes$id %in% signif.hits$X,"group"] <- "DE"

graph <- visNetwork(nodes = nodes, edges = differential.circ.plot$x$edges) %>%
  visIgraphLayout(type = "full", physics = F) %>%
  visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
  visGroups(groupname = "mRNA", shape = "triangle", color = "#0066CC") %>% 
  visGroups(groupname = "DE", color = "#CC3333", shape = "rectangle") %>%
  visLegend()
graph

visNetwork::visSave(graph, file = "DE_SPONGE.html")
