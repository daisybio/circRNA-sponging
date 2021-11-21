#!/usr/bin/env Rscript
library(biomaRt)

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=5) {
  stop("Five arguments must be supplied", call.=FALSE)
}
expression_norm_path = args[1]
output_dir = args[2]
samples_percentage = as.numeric(args[3]) # default 0.2, minimum percentage of samples, a circRNA has to be expressed in is to pass filtering
read_cutoff = as.numeric(args[4]) # default 5, minimum number of reads, a circRNA is required to have to pass filtering
organism = args[5] # organism in three letter code to provide further gene annotation

# define organism three letter codes
org_codes <- list("ebv" = c("Epstein Barr virus", ""),
                  "oar" = c("Ovis aries", ""),
                  "hcmv" = c("Human cytomegalovirus", ""),
                  "cel" = c("Caenorhabditis elegans", ""),
                  "cfa" = c("Canis familiaris", ""),
                  "gga" = c("Gallus gallus", ""),
                  "cgr" = c("Cricetulus griseus", ""),
                  "tgu" = c("Taeniopygia guttata", ""),
                  "xla" = c("Yenips laevis", ""),
                  "ola" = c("Oryzias latipes", ""),
                  "dme" = c("Drosophila melanogaster", "dmelanogaster_gene_ensembl"),
                  "bmo" = c("Bombyx mori", ""),
                  "mmu" = c("Mus musculus", "mmusculus_gene_ensembl"),
                  "rno" = c("Rattus norvegicus", "rnorvegicus_gene_ensembl"),
                  "dre" = c("Dano rerio", ""),
                  "hsa" = c("Homo sapiens", "hsapiens_gene_ensembl"),
                  "kshv" = c("Kaposi sarcoma-associated herpesvirus", "Herpes"),
                  "osa" = c("Oryza sativa", ""),
                  "ssc" = c("Sus scrofa", ""),
                  "bta" = c("Bos taurus", "btaurus_gene_ensembl"),
                  "ath" = c("Arabidopsis thaliana", ""),
                  "xtr" = c("Xenopus tropicalis", ""))
# init organism
org_data <- org_codes[organism][[1]]

expression_norm <- read.table(expression_norm_path, sep = "\t", header=T, stringsAsFactors = F, check.names = F)

samples <- colnames(expression_norm)[-c(1:5)]

# filter data: counts > 5 in at least 20% of samples
if(length(samples) < 5){
  stop("Cannot perform filtering on less than 5 samples")
}
sample_nr_cutoff <- ceiling(samples_percentage *length(samples))

rows_to_keep <- c()
for (i in 1:nrow(expression_norm)){
  number_of_samples_containing_this_circRNA <- 0
  for (j in 6:ncol(expression_norm)){
    if(expression_norm[i,j] >= read_cutoff){
      number_of_samples_containing_this_circRNA <- number_of_samples_containing_this_circRNA + 1
    }
  }
  if(number_of_samples_containing_this_circRNA >= sample_nr_cutoff){
    rows_to_keep <- append(rows_to_keep, i)
  }
}
filtered_data <- expression_norm[rows_to_keep,]
# annotate filtered data
targets <- filtered_data$gene_symbol
targets <- targets[!duplicated(targets)]
# init mart
mart <- biomaRt::useDataset(org_data[2], useMart("ensembl"))
# annotate data: gene symbol + ENSG
gene.ens.all <- biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values=targets,
                      mart=mart,
                      checkFilters = T)
gene.ens.all <- gene.ens.all[!duplicated(gene.ens.all$hgnc_symbol),]
# append data
filtered_data <- merge(filtered_data, gene.ens.all, by.x = "gene_symbol", by.y = "hgnc_symbol", all.x = T)
# write final output
write.table(filtered_data, paste(output_dir, "circRNA_counts_filtered.tsv", sep = "/"), quote = F, sep = "\t", row.names = F)

