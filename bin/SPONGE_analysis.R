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
# convert geneA and geneB
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network, gene.ens.all, by.x = "geneA", by.y = 1, all.x = T)
differential.circ.in.ce.network[!is.na(differential.circ.in.ce.network$hgnc_symbol),"geneA"] <- differential.circ.in.ce.network$hgnc_symbol[!is.na(differential.circ.in.ce.network$hgnc_symbol)]
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network, gene.ens.all, by.x = "geneB", by.y = 1, all.x = T)
differential.circ.in.ce.network[!is.na(differential.circ.in.ce.network$hgnc_symbol.y),"geneB"] <- differential.circ.in.ce.network$hgnc_symbol.y[!is.na(differential.circ.in.ce.network$hgnc_symbol.y)]

differential.circ.plot <- sponge_plot_network(differential.circ.in.ce.network, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
differential.circ.plot$x$edges$label <- paste("mscor:", round(differential.circ.in.ce.network$mscor, 2))

hgncs <- merge(signif.hits, gene.ens.all, by.x = "X", by.y = "ensembl_gene_id")
nodes <- differential.circ.plot$x$nodes
nodes <- merge(nodes, gene.ens.all[,c("ensembl_gene_id", "gene_biotype")], by = 1, all.x = T)

# add circRNA as biotype
nodes[is.na(nodes$gene_biotype),"gene_biotype"] <- "circRNA"
biotypes <- unique(nodes$gene_biotype)
style <- data.frame(groupname=biotypes, 
                    shape=seq(1,length(biotypes)), 
                    color = as.vector(met.brewer("Juarez", length(biotypes))))
# remove preset color and shape
nodes <- nodes[,-c(3,4)]
# change to group
colnames(nodes)[7] <- "group"
# change edges
edges <- differential.circ.plot$x$edges

graph <- visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(type = "full", physics = F) %>%
  visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
  visGroups(groupname = "mRNA", shape = "triangle", color = "#0066CC") %>% 
  visGroups(groupname = "DE", color = "#CC3333", shape = "rectangle") %>%
  visLegend()
graph

visNetwork::visSave(graph, file = "DE_SPONGE.html")
