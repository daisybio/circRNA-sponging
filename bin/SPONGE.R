# install necessary packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SPONGE")
# imports
library("SPONGE")

print(1)