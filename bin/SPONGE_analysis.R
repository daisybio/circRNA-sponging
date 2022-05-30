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

# break types of RNA down to lncRNA and coding
gene.ens.all$gene_biotype[gene.ens.all$gene_biotype != "protein_coding" & gene.ens.all$gene_biotype != "lncRNA"] <- "other_RNA"

hgncs <- merge(signif.hits, gene.ens.all, by.x = "X", by.y = "ensembl_gene_id")
nodes <- differential.circ.plot$x$nodes
nodes <- merge(nodes, gene.ens.all, by = 1, all.x = T)
# change to hgnc
nodes[!is.na(nodes$hgnc_symbol),1:2] <- nodes[!is.na(nodes$hgnc_symbol),"hgnc_symbol"]

# add circRNA as biotype
nodes[is.na(nodes$gene_biotype),"gene_biotype"] <- "circRNA"
biotypes <- unique(nodes$gene_biotype)
style <- data.frame(groupname=biotypes, 
                    shape=seq(1,length(biotypes)), 
                    color = as.vector(met.brewer("Juarez", length(biotypes))))
# remove preset color and shape
nodes <- nodes[,-c(3,4)]
# change to group
colnames(nodes)[8] <- "group"
# add new shape and color
# nodes <- merge(nodes, style, by.x = "group", by.y = "groupname", all.x = T)
# mark differentially expressed RNAs
nodes[nodes$id %in% hgncs$hgnc_symbol | nodes$id %in% signif.hits$X,"group"] <- "DE"
# change edges
edges <- differential.circ.plot$x$edges
# convert geneA and geneB
edges <- merge(edges, gene.ens.all, by.x = "from", by.y = 1, all.x = T)
edges[!is.na(edges$hgnc_symbol),"from"] <- edges$hgnc_symbol[!is.na(edges$hgnc_symbol)]
edges <- merge(edges, gene.ens.all, by.x = "to", by.y = 1, all.x = T)
edges[!is.na(edges$hgnc_symbol.y),"to"] <- edges$hgnc_symbol.y[!is.na(edges$hgnc_symbol.y)]

graph <- visNetwork(nodes = nodes, edges = edges) %>%
  visIgraphLayout(type = "full", physics = F) %>%
  visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
  visGroups(groupname = "protein_coding", shape = "triangle", color = "#0066CC") %>%
  visGroups(groupname = "lncRNA", shape = "square", color = "#F8766D") %>%
  visGroups(groupname = "other_RNA", color = "#DC71FA", shape = "diamond") %>%
  visGroups(groupname = "DE", color = "#CC3333") %>%
  visLegend()
graph

visNetwork::visSave(graph, file = "DE_SPONGE.html")
