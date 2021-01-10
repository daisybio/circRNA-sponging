#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied", call.=FALSE)
}
output_dir = args[1]

file <- paste0(output_dir, "/results/binding_sites/output/circRNA_bind_sites_results.txt")
plot_folder <- paste0(output_dir, "results/binding_sites/plots/")
raw_bindSites <- read.table(file, header = T, sep = '\t', stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]

allBindSites <- count(bindSites, Target, name="freq")
distinct <- distinct(bindSites)
distinctBindSites <- count(distinct, Target, name="freq")

write.table(raw_bindSites, file = paste0(output_dir, "/results/binding_sites/output/bindsites_raw.tsv"), sep = "\t", quote = F, row.names = F)
write.table(allBindSites, file = paste0(output_dir, "/results/binding_sites/output/bindsites_per_circRNA.tsv"), sep = "\t", quote = F, row.names = F)

# miRanda score distribution
scores <- raw_bindSites[,c(1,3)]

# fiter 25% worst scores
filtered_scores <- raw_bindSites[raw_bindSites$Score > quantile(raw_bindSites$Score, 0.25),]
write.table(filtered_scores, file = paste0(output_dir, "/results/binding_sites/output/bindsites_25%_filtered.tsv"), sep = "\t", quote = F, row.names = F)
filtered_scores <- filtered_scores[,c(1,2)]

# redo plots after filtering
allBindSites <- count(filtered_scores, Target, name="freq")
distinct <- distinct(filtered_scores)
distinctBindSites <- count(distinct, Target, name="freq")
