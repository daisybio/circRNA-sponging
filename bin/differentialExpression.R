#!/usr/bin/env Rscript

library(ggplot2)
library(ensembldb)
library(tximport)
library(pheatmap)
library(DESeq2)
library(data.table)

args = commandArgs(trailingOnly = TRUE)

# create ouptut data and plots
create_outputs <- function(d, results, marker, out, filteredRows, nsub=1000) {
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
  png(filename = file.path(out, paste(pca_name, "png", sep = ".")))
  plot(PCA_plot)
  # HEAT MAP
  
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
  # filter for significant differential expression
  signif.hits <- results[!is.na(results$padj) &
                      results$padj<0.10 &
                      abs(results$log2FoldChange)>=1,]
  
  selected <- rownames(signif.hits);selected
  filtered <- as.data.frame(log2(counts(d,normalized=TRUE)[rownames(d) %in% selected,]))
  filtered <- filtered[, df$sample]

  # set output file loc
  heatmap_name <- paste(out, "HMAP", sep = "_")
  # plot heatmap
  pheatmap::pheatmap(filtered, cluster_rows=T, show_rownames=T,
           cluster_cols=F, annotation_col=df,
           filename = file.path(out, paste(heatmap_name, "png", sep = ".")),
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
create_outputs(d = dds, results = res, marker = "condition", out = "total_rna", filteredRows = NULL)
# circRNA gene expression
circ_RNAs <- read.table(file = args[3], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
dds_filtered <- dds
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds_filtered, results = filtered_res, marker = "condition", out = "circ_rna", filteredRows = ens_ids)

# circRNA differential expression
circ_RNA_annotation <- 0
# use annotation if possible
if (length(args)==4) {
  circ_RNA_annotation <- data.frame(read.table(args[4], sep = "\t", header = T))
  circ_RNA_annotation <- data.table(pos=paste0(circ_RNA_annotation$chr, ":", circ_RNA_annotation$start, ":", circ_RNA_annotation$stop, ":", circ_RNA_annotation$strand),
                                    circRNA.ID=circ_RNA_annotation$circRNA.ID)
  tmp <- data.table(pos=paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand))
  circ_RNAs <- circ_RNAs[rownames(merge(tmp, circ_RNA_annotation, by = "pos")),]
} else {
  circ_RNA_annotation <- data.table(circRNA.ID=paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand))
}
circ_RNAs$circRNA.ID <- circ_RNA_annotation$circRNA.ID
circ_expr <- circ_RNAs[,-c(1:6)]
rownames(circ_expr) <- circ_expr$circRNA.ID
circ_expr$circRNA.ID <- NULL
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

create_outputs(dds.circ, res.circ, marker = "condition", out = "circ_rna_DE", filteredRows = NULL, nsub = 100)

# save R image
save.image(file = "DESeq2.RData")
