#!/usr/bin/env Rscript

install.packages("pacman")
pacman::p_load("ggplot2", "ensembldb", "pheatmap", "DESeq2", "data.table", "EnhancedVolcano", "argparser", "MetBrewer") 

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for differenial expression analysis", name = "DE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--samplesheet", help = "Meta data for expressions in tsv")
parser <- add_argument(parser, "--circ_filtered", help = "circRNA filtered expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--circ_raw", help = "circRNA raw expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--tpm_map", help = "TPM map of circular and linear transcripts provided by pipeline")
# FLAGS
parser <- add_argument(parser, "--tpm", help = "Use TPM instead of counts", flag = T)

argv <- parse_args(parser, argv = args)

# create output data and plots
create_outputs <- function(d, results, marker, out, nsub=1000, n = 20, padj = 0.1, log2FC = 0, pseudocount = 0, filter = NULL, isLogFransformed = F, palette = "Renoir") {
  # create dirs in cwd
  dir.create(out, showWarnings = FALSE)
  # col data
  df <- as.data.frame(colData(d))
  df <- df[,c("sample", marker)]
  
  conditions <- split(df, df[,marker])
  
  conditions <- unique(df[,marker])
  
  # results <- results(d,
  #                   contrast = c(marker, as.vector(conditions)))
  # results <- lfcShrink(d,
  #                     contrast = c(marker, as.vector(conditions)), res=results, type = 'normal')
  # write data to disk
  write.table(results, file = file.path(out, paste(out, "tsv", sep = ".")), quote = FALSE, sep = "\t", col.names = NA)
  
  # PCA
  # variance stabilizing transformation
  deseq_vst <- DESeq2::vst(d, blind = FALSE, nsub = nsub)
  PCA_plot <- DESeq2::plotPCA(deseq_vst, intgroup = marker)
  png(filename = file.path(out, paste("pca", "png", sep = ".")), res = 200, width = 1024, height = 800)
  plot(PCA_plot)
  dev.off()
  # HEAT MAP
  
  # filter for significant differential expression
  signif.hits <- results[!is.na(results$padj) &
                      results$padj<as.double(padj) &
                        abs(results$log2FoldChange) > log2FC,]
  # filter for specific RNAs
  if (!is.null(filter)){
    cat("using specific filtering for:", filter, "\n")
    signif.hits <- signif.hits[rownames(signif.hits) %in% filter,]
  }
  cat(nrow(signif.hits), "of", nrow(results), "hits survived filtering of padj:", padj, "and log2FC:", log2FC, "\n")
  # VOLCANO PLOT
  volcano <- EnhancedVolcano(signif.hits,
                             lab = rownames(signif.hits),
                             x = 'log2FoldChange',
                             y = 'pvalue')
  png(filename = file.path(out, paste("volcano", "png", sep = ".")), res = 200, width = 1024, height = 800)
  plot(volcano)
  dev.off()
  
  write.table(signif.hits, file = file.path(out, paste(out, "signif", "tsv", sep = ".")), quote = FALSE, sep = "\t", col.names = NA)
  
  # select all
  # PSEUDOCOUNTS
  selected <- rownames(signif.hits)
  counts <- counts(d, normalized = T)[rownames(d) %in% selected,]+pseudocount
  if (!isLogFransformed) {
    counts <- log2(counts)
  }
  filtered <- as.data.frame(counts)
  filtered <- filtered[, df$sample]

  # set colors
  colors <- c(colorRampPalette(c("blue", "orange"))(100), colorRampPalette(c("orange", "red"))(100))
  #colors <- hcl.colors(101, rev = T)
  #annotation.colors <- hcl.colors(length(conditions), palette = hcl.pals(type = "diverging")[12])
  annotation.colors <- met.brewer(palette, n = length(conditions), type = "continuous")
  
  names(annotation.colors) <- conditions
  annotation.colors <- list(condition = annotation.colors)
  
  # plot total counts per sample
  png(filename = file.path(out, "hits.per.sample.png"), res = 200, width = 1300, height = 800)
  cons <- split(df, df[,marker])
  counts.per.condition <- as.data.frame(lapply(cons, function(x) colSums(counts(d)[,x[["sample"]]] != 1)))
  par(mar=c(5,4,4,5), xpd = T)
  matplot(counts.per.condition, type = "l", xaxt="n", yaxt="n", ylab = NA, main = "unique RNAs per sample", lty = 1, col = annotation.colors, lwd = 2)
  axis(1, at=1:nrow(counts.per.condition), labels = rownames(counts.per.condition), las = 2)
  axis(2, las = 2, at = seq(min(counts.per.condition), max(counts.per.condition), 1000))
  legend("topright", legend = colnames(counts.per.condition), col = annotation.colors, lty = c(1,1), cex=0.8, lwd = 2, inset = c(-0.2, 0))
  par(xpd=F)
  abline(v = 1:nrow(counts.per.condition), lty = 2, col = "grey")
  abline(h = seq(min(counts.per.condition), max(counts.per.condition), 1000), lty = 2, col = "grey")
  dev.off()
  
  row.names(df) <- df$sample
  df <- df[, marker, drop = F]
  # plot heatmap
  pheatmap::pheatmap(filtered, cluster_rows=T, show_rownames=F,
           cluster_cols=T, annotation_col=df,
           filename = file.path(out, paste("HMAP", "png", sep = ".")),
           height = 15, width = 25, legend = T, annotation_legend = T,
           show_colnames = F, color = colors, annotation_names_col = F, main = out,
           treeheight_row = 0, treeheight_col = 0, 
           annotation_colors = annotation.colors, 
           fontsize = 25)
}

# read gene expression and add pseudocount
gene_expression <- as.matrix(read.table(file = argv$gene_expr, header = T, sep = "\t", stringsAsFactors = F, check.names = F))

# metadata
samplesheet <- read.table(file = argv$samplesheet, sep = "\t", header = T)

# circRNAs filtered
circ_RNAs <- read.table(file = argv$circ_filtered, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
filtered.circs <- paste0(circ_RNAs$chr, ":", circ_RNAs$start, "-", circ_RNAs$stop, ":", circ_RNAs$strand)
# raw circRNA expression
circ.raw <- read.table(file = argv$circ_raw, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
circ.raw$key <- paste0(circ.raw$chr, ":", circ.raw$start, "-", circ.raw$stop, ":", circ.raw$strand)

annotation <- "circBaseID" %in% colnames(circ_RNAs)
# get samples
samples <- samplesheet$sample

# DIFFERENTIAL circRNA EXPRESSION
# use given annotation if possible
if (annotation) {
  circ_RNA_annotation <- ifelse(circ_RNAs$circBaseID != "None", 
                                circ_RNAs$circBaseID, 
                                paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand))
} else {
  circ_RNA_annotation <- paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand)
}

# get raw expression values for filtered circRNAs and samples
circ_expr <- circ.raw[circ.raw$key %in% filtered.circs, samples]
rownames(circ_expr) <- circ_RNA_annotation

# use tpms instead of counts
if (argv$tpm) {
  # convert circRNA and linear expression
  TPM.map <- read.table(argv$tpm_map, header = T, sep = "\t", stringsAsFactors = F, check.names = F)
  gene_expression_tpm <- as.matrix(TPM.map[rownames(gene_expression), colnames(gene_expression)])
  circ_expr_tpm <- as.matrix(TPM.map[rownames(circ_expr), colnames(circ_expr)])
  
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(gene_expression_tpm),
                                        colData = samplesheet,
                                        design = ~ condition)
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  # sort by p-value
  res <- res[order(res$padj),]
  DESeq2::summary(res)
  
  create_outputs(d = dds, results = res, marker = "condition", out = "total_rna_tpm", nsub = 100, isLogFransformed = argv$tpm)
  # CIRCULAR RNA
  
  dds.circ <- DESeq2::DESeqDataSetFromMatrix(countData = round(circ_expr_tpm),
                                             colData = samplesheet,
                                             design = ~ condition)
  dds.circ <- DESeq2::DESeq(dds.circ)
  res.circ <- DESeq2::results(dds.circ)
  # sort by p-value
  res.circ <- res.circ[order(res.circ$padj),]
  # create summary
  DESeq2::summary(res.circ)
  create_outputs(dds.circ, res.circ, marker = "condition", out = "circ_rna_DE_tpm", nsub = 100, isLogFransformed = argv$tpm)
}
pseudocount = 1
gene_expression <- gene_expression + pseudocount
circ_expr <- as.matrix(circ_expr) + pseudocount

# TOTAL RNA
dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(gene_expression),
                                      colData = samplesheet,
                                      design = ~ condition)
dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds)
# sort by p-value
res <- res[order(res$padj),]
DESeq2::summary(res)

create_outputs(d = dds, results = res, marker = "condition", out = "total_rna", nsub = 100)
# CIRCULAR RNA

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
