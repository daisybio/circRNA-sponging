#!/usr/bin/env Rscript

args = commandArgs(trailingOnly = TRUE)

sample_files <- list.files(path = args[1], full.names = T, pattern = "gene_expression\\.tsv$")
out_file <- args[2]

print(sample_files)

total_expr <- data.frame()
for (idx in seq_along(sample_files)) {
  file <- sample_files[idx]
  gene_expr <- data.frame(read.table(file = file, header = T, sep = "\t"))
  if (idx == 1) {
    total_expr <- gene_expr
  } else {
    total_expr <- merge(total_expr, gene_expr, by = "X", all = T)
  }
}
write.table(total_expr, file = out_file, quote = F, sep = "\t")
