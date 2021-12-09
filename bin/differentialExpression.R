#!/usr/bin/env Rscript

library(ggplot2)
library(ensembldb)
library(tximport)
library(pheatmap)
library(DESeq2)

args = commandArgs(trailingOnly = TRUE)

# create ouptut data and plots
create_outputs <- function(d, results, marker, out, filteredRows) {
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
  # group
  df <- as.data.frame(colData(d))
  df <- df[order(as.character(df$condition)),]
  d@colData@listData <- as.list(df)
  d@colData@rownames <- rownames(df)
  df <- df[,c("sample", "condition")]
  
  # filter
  if (is.null(filteredRows)) {
    top_genes <- head(order(rowVars(assay(deseq_vst)), decreasing = T), 20)
    filtered <- assay(deseq_vst)[top_genes,]
  } else {
    filtered <- assay(deseq_vst)
    filtered <- filtered[row.names(filtered) %in% filteredRows,]
    top_genes <- head(order(rowVars(filtered), decreasing = T), 20)
    filtered <- filtered[top_genes,]
  }

  # set output file loc
  heatmap_name <- paste(out, "HMAP", sep = "_")
  # plot heatmap
  pheatmap::pheatmap(filtered, cluster_rows=T, show_rownames=T,
           cluster_cols=T, annotation_col=df,
           filename = file.path(out, paste(heatmap_name, "png", sep = ".")),
           height = 15, width = 25, legend = F)
}

txi <- readRDS(args[1])

# metadata
samplesheet <- read.table(file = args[2], sep = "\t", header = TRUE)

# dds object
dds <- DESeq2::DESeqDataSetFromTximport(txi,
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
create_outputs(d = dds, results = res, marker = "condition", out = "total_rna", filteredRows = NULL)
# circRNA only
circ_RNAs <- read.table(file = args[3], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
dds_filtered <- dds
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds_filtered, results = filtered_res, marker = "condition", out = "circ_rna", filteredRows = ens_ids)
# save R image
save.image(file = "DESeq2.RData")
