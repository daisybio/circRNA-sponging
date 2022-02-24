#!/usr/bin/env Rscript

library(SPONGE)
library(biomaRt)

args = commandArgs(trailingOnly = TRUE)

# load SPONGE workspace
load(args[1])
# differentially expressed circRNAs
signif.hits <- read.table(args[2], sep = "\t", header = T)

# load biomaRt
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

# get all differentially expressed circRNAs that act as ceRNAs
differential.circ.in.ce.network.A <- merge(ceRNA_interactions_all_circ, signif.hits, by.x = "geneA", by.y = "X")
differential.circ.in.ce.network.B <- merge(ceRNA_interactions_all_circ, signif.hits, by.x = "geneB", by.y = "X")
differential.circ.in.ce.network <- merge(differential.circ.in.ce.network.A, differential.circ.in.ce.network.B, all = T)

# find gene symbols for ensg
targets <- c(differential.circ.in.ce.network$geneA, differential.circ.in.ce.network$geneB)
targets <- unique(targets[!grepl(ifelse(annotation, "circ", "chr"), targets)])
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
# build network for all DE circRNAs
de.circ.network <- c()
for (hit in signif.hits$X) {
  circ.network <- subnetwork(ceRNA_interactions_all_circ, pattern = hit)
  de.circ.network <- rbind(de.circ.network, circ.network)
}
differential.circ.in.ce.network <- de.circ.network

differential.circ.plot <- sponge_plot_network(differential.circ.in.ce.network, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
differential.circ.plot$x$edges$label <- paste("mscor:", round(differential.circ.in.ce.network$mscor, 2))
differential.circ.plot

nodes <- differential.circ.plot$x$nodes
nodes$color <- NULL
nodes$shape <- NULL
nodes$group <- ifelse(grepl("circ", nodes$id), "circRNA", "mRNA")
test <- visNetwork(nodes = nodes, edges = differential.circ.plot$x$edges) %>%
  visIgraphLayout(type = "full", physics = F) %>%
  visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1))) %>%
  visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
  visGroups(groupname = "mRNA", shape = "triangle", color = "#0066CC") %>% 
  visLegend()
test
