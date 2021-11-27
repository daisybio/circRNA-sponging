#!/usr/bin/env Rscript

library("reshape2", "ggplot")

library(DESeq2, ensembldb, "tximport")

args = commandArgs(trailingOnly = TRUE)

# create ouptut data and plots
create_outputs <- function(d, results, marker, out) {
  # create dirs in cwd
  dir.create(out, showWarnings = FALSE)
  # write data to disk
  write.table(cbind(ENS_ID=rownames(results), results), file = file.path(out, paste(out, "tsv", sep = ".")), quote = FALSE, sep = "\t", col.names = NA)
  # PCA
  # variance stabilizing transformation
  deseq_vst <- DESeq2::vst(d, blind = FALSE)
  PCA_plot <- DESeq2::plotPCA(deseq_vst, intgroup = marker)
  pca_name <- paste(out, "pca", sep = "_")
  png(filename = file.path(out, paste(pca_name, "png", sep = ".")))
  plot(PCA_plot)
  # HEAT MAP
  # filter
  ntd <- normTransform(d)
  select <- order(rowMeans(counts(d,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  df <- as.data.frame(colData(d)[,c("condition", "smallRNA")])
  # set output file loc
  heatmap_name <- paste(out, "HMAP", sep = "_")
  # plot heatmap
  pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
           cluster_cols=FALSE, annotation_col=df, 
           filename = file.path(out, paste(heatmap_name, "png", sep = ".")),
           height = 25, width = 25)
  # close device
  dev.off
}

# read gene expression counts
samples_loc <- args[1]

# metadata
samplesheet <- read.table(file = args[2], sep = "\t", header = TRUE)

# get quant files
quant.files <- file.path(samples_loc, samplesheet$sample, "salmon", "quant.sf")
names(quant.files) <- samplesheet$sample

# read input gtf
txdb <- ensembldb::ensDbFromGtf(args[3])
tx <- ensembldb::EnsDb(txdb)
tx2gene <- ensembldb::transcripts(tx, return.type="DataFrame")
tx2gene <- tx2gene[, c("tx_name", "gene_id")]
# txi object
txi <- tximport::tximport(quant.files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T)
# write total gene expression over samples to file
counts <- txi$counts
colnames(counts) <- samplesheet$sample
write.table(counts, file = "gene_expression.tsv")
# dds object
dds <- DESeqDataSetFromTximport(txi,
                                colData = samplesheet,
                                design = ~ condition)

# differential expression analysis
dds <- DESeq2::DESeq(dds)
# results
res <- DESeq2::results(dds)
# sort by p-value
res <- res[order(res$padj),]
# create summary
DESeq2::summary(res)

# WRITE OUTPUTS
# total_RNA
create_outputs(d = dds, results = res, marker = "condition", out = "total_rna")
# circRNA only
circ_RNAs <- read.table(file = args[4], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds, results = filtered_res, marker = "condition", out = "circ_rna")
# save R image
save.image(file = "DESeq2.RData")
