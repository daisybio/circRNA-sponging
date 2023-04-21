args <- commandArgs(trailingOnly = T)

psiPerIsoform <- read.csv(args[1], sep = "\t", check.names = F)
meta <- read.csv(args[3], sep = "\t")

# split psi into cell types
sapply(unique(meta$condition), function(ct) {
  s <- meta[meta$condition==ct,"sample"]
  write.table(psiPerIsoform[,s], file = paste0(ct, ".psi"),
              sep = "\t", quote = F)
})