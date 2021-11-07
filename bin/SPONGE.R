# install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SPONGE")
# imports
library("SPONGE")

# gene expression
gene_expr <- t(read.table(file = args[1], header = TRUE, sep = "\t"))
# miRNA expression
mir_expr <- t(read.table(file = args[2], header = TRUE, sep = "\t"))

print(1)
