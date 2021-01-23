#!/usr/bin/env Rscript
library(ggplot2)
library(plyr)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=1) {
  stop("One argument must be supplied", call.=FALSE)
}
raw_bindSites_file = args[1]

raw_bindSites <- read.table(raw_bindSites_file, header = T, sep = '\t', stringsAsFactors = F)
bindSites <- raw_bindSites[,c(1,2)]

allBindSites <- count(bindSites, Target, name="freq")
distinct <- distinct(bindSites)
distinctBindSites <- count(distinct, Target, name="freq")

write.table(allBindSites, file = paste0("bindsites_per_circRNA.tsv"), sep = "\t", quote = F, row.names = F)

# miRanda score distribution
scores <- raw_bindSites[,c(1,3)]

# fiter 25% worst scores
filtered_scores <- raw_bindSites[raw_bindSites$Score > quantile(raw_bindSites$Score, 0.25),]
write.table(filtered_scores, file = paste0("bindsites_25%_filtered.tsv"), sep = "\t", quote = F, row.names = F)
