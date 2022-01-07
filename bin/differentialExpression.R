#!/usr/bin/env Rscript

library(ggplot2)
library(ensembldb)
library(tximport)
library(pheatmap)
library(DESeq2)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

# create ouptut data and plots
create_outputs <- function(d, results, marker, out, nsub=1000, n = 20, padj = 0.1, log2FC = 1) {
  # create dirs in cwd
  dir.create(out, showWarnings = FALSE)
  # write data to disk
  write.table(cbind(ENS_ID=rownames(results), results), file = file.path(out, paste(out, "tsv", sep = ".")), quote = FALSE, sep = "\t", col.names = NA)
  # group
  df <- as.data.frame(colData(d))
  df <- df[,c("sample", "condition")]
  df <- df[order(df$condition, decreasing = T),]
  # PCA
  # variance stabilizing transformation
  deseq_vst <- DESeq2::vst(d, blind = FALSE, nsub = nsub)
  PCA_plot <- DESeq2::plotPCA(deseq_vst, intgroup = marker)
  pca_name <- paste(out, "pca", sep = "_")
  png(filename = file.path(out, paste(pca_name, "png", sep = ".")), res = 200, width = 1024, height = 800)
  plot(PCA_plot)
  # HEAT MAP
  
  # filter for significant differential expression
  signif.hits <- results[!is.na(results$padj) &
                      results$padj<as.double(padj) &
                      abs(results$log2FoldChange)>=as.numeric(log2FC),]
  
  signif.top <- head(signif.hits, n)
  # select all
  # PSEUDOCOUNTS
  selected <- rownames(signif.hits)
  counts <- counts(d,normalized=T)[rownames(d) %in% selected,]+1e-3
  filtered <- as.data.frame(log2(counts))
  filtered <- filtered[, df$sample]
  # select top n
  selected.top <- rownames(signif.top)
  counts.top <- counts(d,normalized=TRUE)[rownames(d) %in% selected.top,]+1e-6
  filtered.top <- as.data.frame(log2(counts.top))
  filtered.top <- filtered.top[, df$sample]

  # set output file loc
  heatmap_name <- paste(out, "HMAP", sep = "_")
  # plot heatmap
  pheatmap::pheatmap(filtered, cluster_rows=T, show_rownames=T,
           cluster_cols=F, annotation_col=df,
           filename = file.path(out, paste(heatmap_name, "png", sep = ".")),
           height = 15, width = 25, legend = T)
  # top n
  heatmap_name_top <- paste(out, "HMAP_top", sep = "_")
  pheatmap::pheatmap(filtered.top, cluster_rows=T, show_rownames=T,
                     cluster_cols=F, annotation_col=df,
                     filename = file.path(out, paste(heatmap_name_top, "png", sep = ".")),
                     height = 15, width = 25, legend = T)
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

# WRITE OUTPUTS GENE EXPRESSION
# total_RNA
create_outputs(d = dds, results = res, marker = "condition", out = "total_rna")
# circRNA gene expression
circ_RNAs <- read.table(file = args[3], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_id
dds_filtered <- dds
filtered_res <- res[rownames(res) %in% ens_ids,]
create_outputs(d = dds_filtered, results = filtered_res, marker = "condition", out = "circ_rna_GE")

# DIFFERENTIAL CIRCRNA EXPRESSION
# use given annotation if possible
if ("circBaseID" %in% colnames(circ_RNAs)) {
  circ_RNA_annotation <- data.table(circRNA.ID=circ_RNAs$circBaseID)
  # cut table and annotate rownames
  circ_expr <- circ_RNAs[,-c(1:8)]
} else {
  circ_RNA_annotation <- data.table(circRNA.ID=paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand))
  # cut table and annotate row names
  circ_expr <- circ_RNAs[,-c(1:7)]
}
rownames(circ_expr) <- circ_RNA_annotation$circRNA.ID
# format miRNA ids
colnames(circ_expr) <- sapply(gsub("\\.", "-", colnames(circ_expr)), "[", 1)
circ_expr <- as.matrix(circ_expr)

dds.circ <- DESeq2::DESeqDataSetFromMatrix(countData = round(circ_expr),
                                           colData = samplesheet,
                                           design = ~ condition)
dds.circ <- DESeq2::DESeq(dds.circ)
res.circ <- DESeq2::results(dds.circ)
# sort by p-value
res.circ <- res.circ[order(res.circ$padj),]
# create summary
DESeq2::summary(res.circ)

create_outputs(dds.circ, res.circ, marker = "condition", out = "circ_rna_DE", nsub = 100)

# save R image
save.image(file = "DESeq2.RData")
