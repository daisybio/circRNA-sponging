#!/usr/bin/env Rscript

library(SPONGE)
library(biomaRt)
library(visNetwork)

args = commandArgs(trailingOnly = TRUE)

# load SPONGE workspace
load(args[1])
# differentially expressed circRNAs
signif.hits.circ <- read.table(args[2], sep = "\t", header = T)
signif.hits.mRNA <- read.table(args[3], sep = "\t", header = T)

signif.hits <- rbind(signif.hits.circ, signif.hits.mRNA)

# load biomaRt
mart <- NULL
not_done <- T
while(not_done){
  tryCatch(
    {
      mart <- useEnsembl(biomart = "ensembl", "hsapiens_gene_ensembl")
    },
    error=function(cond) {
      message(cond)
    },
    warning=function(cond) {
      message(cond)
    },
    finally={
      not_done = F
    }
  )
}
print("SUCCESS")

# get all differentialy expressed circRNAs that act as ceRNAs
# build network for all DE circRNAs
de.circ.network <- c()
for (hit in signif.hits$X) {
  circ.network <- subnetwork(ceRNA_interactions_fdr, pattern = hit)
  de.circ.network <- rbind(de.circ.network, circ.network)
}
differential.circ.in.ce.network <- de.circ.network
# find gene symbols for ensg
targets <- unique(c(differential.circ.in.ce.network$geneA, differential.circ.in.ce.network$geneB))
targets <- targets[!grepl(ifelse(annotation, "circ", "chr"), targets)]
# start query
gene.ens.all <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                               filter = "ensembl_gene_id",
                               values = targets,
                               mart = mart,
                               checkFilters = T)
# merge for geneA and geneB
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network, gene.ens.all, by.x = "geneA", by.y = 1, all.x = T)
differential.circ.in.ce.network[!is.na(differential.circ.in.ce.network$hgnc_symbol),"geneA"] <- differential.circ.in.ce.network$hgnc_symbol[!is.na(differential.circ.in.ce.network$hgnc_symbol)]
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network, gene.ens.all, by.x = "geneB", by.y = 1, all.x = T)
differential.circ.in.ce.network[!is.na(differential.circ.in.ce.network$hgnc_symbol.y),"geneB"] <- differential.circ.in.ce.network$hgnc_symbol.y[!is.na(differential.circ.in.ce.network$hgnc_symbol.y)]

differential.circ.plot <- sponge_plot_network(differential.circ.in.ce.network, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
differential.circ.plot$x$edges$label <- paste("mscor:", round(differential.circ.in.ce.network$mscor, 2))
differential.circ.plot

hgncs <- merge(signif.hits, gene.ens.all, by.x = "X", by.y = "ensembl_gene_id")
nodes <- differential.circ.plot$x$nodes
nodes$color <- NULL
nodes[nodes$id %in% hgncs$hgnc_symbol | nodes$id %in% signif.hits$X,"color"] <- "#CC3333"
nodes$shape <- ifelse(grepl("c", nodes$id), "rectangle", "triangle")
nodes$group <- ifelse(grepl("c", nodes$id), "circRNA", "mRNA")
nodes[nodes$id %in% hgncs$hgnc_symbol | nodes$id %in% signif.hits$X,"group"] <- "DE"

graph <- visNetwork(nodes = nodes, edges = differential.circ.plot$x$edges) %>%
  visIgraphLayout(type = "full", physics = F) %>%
  visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
  visGroups(groupname = "mRNA", shape = "triangle", color = "#0066CC") %>% 
  visGroups(groupname = "DE", color = "#CC3333", shape = "rectangle") %>%
  visLegend()
graph

visNetwork::visSave(graph, file = "DE_SPONGE.html")
