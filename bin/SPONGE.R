# install necessary packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SPONGE")
args = commandArgs(trailingOnly = TRUE)
# imports
# library("SPONGE")
library(biomaRt, "tidyverse")

args = c("/Users/leonschwartz/Desktop/Bioinformatik/local_data/pipeline_expample_output/output/results/differential_expression/gene_expression.tsv",
         "/Users/leonschwartz/Desktop/Bioinformatik/local_data/pipeline_expample_output/output/results/miRNA/miRNA_counts_filtered.tsv",
         "/Users/leonschwartz/Desktop/Bioinformatik/local_data/miRNA/miRTarBase_MTI.csv",
         "hsa")

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

# create count matrix
count_matrix <- function(miR, g_ids, frame) {
  df <- data.frame(matrix(ncol=length(miR),nrow=length(g_ids), dimnames=list(g_ids, NULL)))
  colnames(df) <- miR
  for (id in g_ids) {
    matches <- frame[frame$Target.Gene == id,]$miRNA
    for (match in matches) {
      if (!is.na(df[id, match])) {
        df[id, match] = df[id, match] + 1
      } else {
        df[id, match] = 1
      }
    }
  }
  df[is.na(df)] <- 0
  return(df)
}

gene_expr <- as.data.frame(read.table(file = args[1], header = TRUE, sep = "\t"))
mi_rna_expr <- t(as.data.frame(read.table(file = args[2], header = TRUE, sep = "\t")))

target_scan_symbols <- as.data.frame(read.csv(file = args[3], header = TRUE))
org_data <- org_codes[args[4]][[1]]
# filter by organism
target_scan_symbols <- target_scan_symbols[target_scan_symbols$Species..miRNA. == org_data[1],]

genes <- gene_expr$X
simple_genes <- c()
for (i in seq_along(genes)) {
  simple_genes[i] = sub("\\..*", "", genes[i])
}
gene_expr$X <- simple_genes
genes <- simple_genes

# set up mart
mart <- useDataset(org_data[2], useMart("ensembl"))
gene_names <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = genes, mart = mart)
# converting gene ids and reformatting file
gene_expr <- merge(gene_expr, gene_names, by.x = "X", by.y = "ensembl_gene_id")
gene_expr$X <- NULL
gene_expr <- gene_expr %>%
  relocate(hgnc_symbol.x)
gene_expr <- t(gene_expr)
# extract all unique miRNAs
miRNAs <- target_scan_symbols$miRNA
miRNAs <- miRNAs[!duplicated(miRNAs)]
# extract all unique gene ids
gene_ids <- target_scan_symbols$Target.Gene
gene_ids <- gene_ids[!duplicated(gene_ids)]
# create target scan counts
target_scan_symbols <- count_matrix(miR = miRNAs, g_ids = gene_ids, frame = target_scan_symbols)


