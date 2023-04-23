args <- commandArgs(trailingOnly = T)

tpms <- read.csv(args[1], sep = "\t", check.names = F)
meta <- read.csv(args[2], sep = "\t")

# split psi into cell types
sapply(unique(meta$condition), function(ct) {
  s <- meta[meta$condition==ct,"sample"]
  write.table(tpms[,s], file = paste0(ct, ".tsv"),
              sep = "\t", quote = F)
})