#!/usr/bin/env Rscript

library(ggplot2)
library(ensembldb)
library(tximport)
library(pheatmap)
library(DESeq2)
library(data.table)
library(EnhancedVolcano)

args = commandArgs(trailingOnly = TRUE)

# create output data and plots
create_outputs <- function(d, results, marker, out, nsub=1000, n = 20, padj = 0.1, log2FC = 1, pseudocount = 1e-3) {
  # create dirs in cwd
  dir.create(out, showWarnings = FALSE)
  results <- results(d,
                     contrast = c('condition', 'Non-Seminoma', 'Seminoma'))
  results <- lfcShrink(d,
                       contrast = c('condition', 'Non-Seminoma', 'Seminoma'), res=results, type = 'normal')
  # write data to disk
  write.table(cbind(ENS_ID=rownames(results), results), file = file.path(out, paste(out, "tsv", sep = ".")), quote = FALSE, sep = "\t", col.names = NA)
  # col data
  df <- as.data.frame(colData(d))
  df <- df[,c("sample", "condition")]
  
  # PCA
  # variance stabilizing transformation
  deseq_vst <- DESeq2::vst(d, blind = FALSE, nsub = nsub)
  PCA_plot <- DESeq2::plotPCA(deseq_vst, intgroup = marker)
  png(filename = file.path(out, paste("pca", "png", sep = ".")), res = 200, width = 1024, height = 800)
  plot(PCA_plot)
  dev.off()
  # VOLCANO PLOT
  volcano <- EnhancedVolcano(results,
                  lab = rownames(results),
                  x = 'log2FoldChange',
                  y = 'pvalue')
  png(filename = file.path(out, paste("volcano", "png", sep = ".")), res = 200, width = 1024, height = 800)
  plot(volcano)
  dev.off()
  # HEAT MAP
  
  # filter for significant differential expression
  signif.hits <- results[!is.na(results$padj) &
                      results$padj<as.double(padj) &
                        abs(results$log2FoldChange) > log2FC,]
  
  signif.top <- head(signif.hits[order(signif.hits$padj),], n)
  # select all
  # PSEUDOCOUNTS
  selected <- rownames(signif.hits)
  counts <- counts(d, normalized = T)[rownames(d) %in% selected,]+pseudocount
  filtered <- as.data.frame(log10(counts))
  filtered <- filtered[, df$sample]
  # select top n
  selected.top <- rownames(signif.top)
  counts.top <- counts(d, normalized = T)[rownames(d) %in% selected.top,]+pseudocount
  filtered.top <- as.data.frame(log2(counts.top))
  filtered.top <- filtered.top[, df$sample]

  # set output file loc
  # plot heatmap
  d <- dist(t(filtered))
  m <- 5 / max(d)
  # set colors
  colors <- c(colorRampPalette(c("blue", "orange"))(100), colorRampPalette(c("orange", "red"))(100))
  annotation.colors <- list(condition = c("Seminoma" = "#339300", "Non-Seminoma" = "#CC0000"))
  row.names(df) <- df$sample
  df <- df[, "condition", drop = F]
  pheatmap::pheatmap(filtered, cluster_rows=T, show_rownames=F,
           cluster_cols=T, annotation_col=df,
           filename = file.path(out, paste("HMAP", "png", sep = ".")),
           height = 15, width = 25, legend = T, annotation_legend = T,
           show_colnames = F, color = colors, annotation_names_col = F, main = out,
           treeheight_row = 0, treeheight_col = 0, annotation_colors = annotation.colors)
  # top n
  pheatmap::pheatmap(filtered.top, cluster_rows=T, show_rownames=T,
                     cluster_cols=T, annotation_col=df,
                     filename = file.path(out, paste("HMAP_top", "png", sep = ".")),
                     height = 15, width = 25, legend = F,
                     color = colors, annotation_colors = annotation.colors)
}
# load total gene expression of samples
txi <- readRDS(args[1])

# metadata
samplesheet <- read.table(file = args[2], sep = "\t", header = T)

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

# circRNAs filtered
circ_RNAs <- read.table(file = args[3], sep = "\t", header = T, stringsAsFactors = F, check.names = F)
filtered.circs <- paste0(circ_RNAs$chr, ":", circ_RNAs$start, "-", circ_RNAs$stop, ":", circ_RNAs$strand)
# raw circRNA expression
circ.raw <- read.table(file = args[4], sep = "\t", header = T, stringsAsFactors = F, check.names = F)
circ.raw$key <- paste0(circ.raw$chr, ":", circ.raw$start, "-", circ.raw$stop, ":", circ.raw$strand)

ens_ids <- circ_RNAs$ensembl_gene_id
dds_filtered <- dds
filtered_res <- res[rownames(res) %in% ens_ids,]
create_outputs(d = dds_filtered, results = filtered_res, marker = "condition", out = "circ_rna_GE")

annotation <- "circBaseID" %in% colnames(circ_RNAs)
# get samples
samples <- samplesheet$sample

# DIFFERENTIAL CIRCRNA EXPRESSION
# use given annotation if possible
if (annotation) {
  circ_RNA_annotation <- data.table(circRNA.ID=circ_RNAs$circBaseID)
} else {
  circ_RNA_annotation <- data.table(circRNA.ID=paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand))
}
# get raw expression values for filtered circRNAs and samples
circ_expr <- circ.raw[circ.raw$key %in% filtered.circs, samples]
rownames(circ_expr) <- circ_RNA_annotation$circRNA.ID
# add pseudo counts
circ_expr <- as.matrix(circ_expr) + 1

dds.circ <- DESeq2::DESeqDataSetFromMatrix(countData = round(circ_expr),
                                           colData = samplesheet,
                                           design = ~ condition)
dds.circ <- DESeq2::DESeq(dds.circ)
res.circ <- DESeq2::results(dds.circ)
# sort by p-value
res.circ <- res.circ[order(res.circ$padj),]
# create summary
DESeq2::summary(res.circ)

create_outputs(dds.circ, res.circ, marker = "condition", out = "circ_rna_DE", nsub = 100, pseudocount = 0, log2FC = 0)

# save R image
save.image(file = "DESeq2.RData")
