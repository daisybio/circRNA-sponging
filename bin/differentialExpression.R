#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load("ggplot2", "ensembldb", "pheatmap", "DESeq2", "data.table", "EnhancedVolcano", "argparser", "MetBrewer") 

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for differenial expression analysis", name = "DE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--samplesheet", help = "Meta data for expressions in tsv")
parser <- add_argument(parser, "--circ_filtered", help = "circRNA filtered expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--circ_raw", help = "circRNA raw expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--tpm_map", help = "TPM map of circular and linear transcripts provided by pipeline")
# parameter settings
parser <- add_argument(parser, "--fdr", help = "FDR threshold", default = 0.01)
parser <- add_argument(parser, "--log2fc", help = "Log2FoldChange threshold", default = 0)
parser <- add_argument(parser, "--palette", help = "Palette to use for MetBrewer", default = "Renoir")

argv <- parse_args(parser, argv = args)

# create output data and plots
create_outputs <- function(d, results, marker, out, nsub=1000, n = 20, padj = 0.1, log2FC = 0, pseudocount = 0, filter = NULL, isLogFransformed = F, palette = "Renoir") {
  # create dirs in cwd
  dir.create(out, showWarnings = FALSE)
  # col data
  df <- as.data.frame(colData(d))
  df <- df[,c("sample", marker)]
  
  conditions <- unique(df[,marker])
  
  # results <- results(d,
  #                   contrast = c(marker, as.vector(conditions)))
  # results <- lfcShrink(d,
  #                     contrast = c(marker, as.vector(conditions)), res=results, type = 'normal')
  # write data to disk
  write.table(results, file = file.path(out, paste(out, "tsv", sep = ".")), quote = FALSE, sep = "\t", col.names = NA)
  
  # HEAT MAP
  
  # filter for significant differential expression
  signif.hits <- results[!is.na(results$padj) &
                           results$padj<as.double(padj) &
                           abs(results$log2FoldChange) > log2FC,]
  
  # DE vs not
  png(filename = file.path(out, "nDE.png"), res = 200, width = 1300, height = 800)
  nDE <- c(rest=(nrow(results)-nrow(signif.hits)), DE=nrow(signif.hits))
  col <- met.brewer(palette, n = length(nDE))
  pie(nDE, labels = nDE, col = col)
  par(mar = c(5, 4, 4, 8), xpd = T)
  legend("topright", legend = names(nDE), inset = c(-0.05, 0), cex = 0.75, fill = col)
  dev.off()
  
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
  png(filename = file.path(out, paste("volcano", "png", sep = ".")), res = 200, width = 1200, height = 1200)
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
  annotation.colors <- met.brewer(palette, n = length(conditions))
  
  names(annotation.colors) <- conditions
  
  # PCA
  # variance stabilizing transformation
  deseq_vst <- DESeq2::vst(d, blind = FALSE, nsub = nsub)
  PCA_data <- DESeq2::plotPCA(deseq_vst, intgroup = marker, returnData = T)
  percentVar <- round(100 * attr(PCA_data, "percentVar"))
  png(filename = file.path(out, paste("pca", "png", sep = ".")), res = 200, width = 1024, height = 800)
  ggplot(PCA_data, aes(x = PC1, y = PC2, color = group)) + 
    geom_point(size = 3) +
    scale_color_manual(values = annotation.colors) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
    ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
    ggtitle(paste0("PCA plot of ", out))
  dev.off()
  
  # reformat meta data for heatmap annotation
  row.names(df) <- df$sample
  df <- df[, marker, drop = F]
  
  # plot heatmap
  pheatmap::pheatmap(filtered, cluster_rows=T, show_rownames=F,
                     cluster_cols=T, annotation_col=df,
                     filename = file.path(out, paste("HMAP", "png", sep = ".")),
                     height = 15, width = 25, legend = T, annotation_legend = T,
                     show_colnames = F, color = colors, annotation_names_col = F, main = out,
                     treeheight_row = 0, treeheight_col = 0, 
                     annotation_colors = list(condition = annotation.colors), 
                     fontsize = 25)
}

# read gene expression and add pseudocount
gene_expression <- as.matrix(read.table(file = argv$gene_expr, header = T, sep = "\t", stringsAsFactors = F, check.names = F))

# metadata
samplesheet <- read.table(file = argv$samplesheet, sep = "\t", header = T)

# plot ratio of conditions
png(filename = "conditions.png", res = 200, width = 1300, height = 800)
condition.occurences <- table(samplesheet$condition)
label <- paste(round(prop.table(condition.occurences)*100), "%", sep = "")
cond.col <- met.brewer(argv$palette, n = length(condition.occurences))
pie(condition.occurences, col = cond.col, labels = label)
par(mar = c(5, 4, 4, 8), xpd = T)
legend("topright", legend = names(condition.occurences), fill = cond.col, inset = c(-0.05, 0), cex = 0.75)
dev.off()

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

# plot ratio if database annotation
png(filename = "db_annotation.png", res = 200, width = 1300, height = 800)
db.rate <- table(grepl("circ", circ_RNA_annotation))
label <- paste(round(prop.table(db.rate)*100), "%", sep = "")
cond.col <- met.brewer(argv$palette, n = length(db.rate))
pie(db.rate, col = cond.col, labels = label, main = "circRNA database annotation")
par(mar = c(5, 4, 4, 8), xpd = T)
legend("topright", legend = names(db.rate), fill = cond.col, inset = c(-0.05, 0), cex = 0.75)
dev.off()

# get raw expression values for filtered circRNAs and samples
samples <- intersect(colnames(circ.raw), samples)
circ_expr <- circ.raw[circ.raw$key %in% filtered.circs, samples]
rownames(circ_expr) <- circ_RNA_annotation

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

create_outputs(d = dds, results = res, marker = "condition", out = "total_rna", nsub = 100,
               padj = argv$fdr, log2FC = argv$log2fc, palette = argv$palette)
# CIRCULAR RNA

samplesheet <- samplesheet[samplesheet$sample %in% samples,]
dds.circ <- DESeq2::DESeqDataSetFromMatrix(countData = round(circ_expr),
                                           colData = samplesheet,
                                           design = ~ condition)
dds.circ <- DESeq2::DESeq(dds.circ)
res.circ <- DESeq2::results(dds.circ)
# sort by p-value
res.circ <- res.circ[order(res.circ$padj),]
# create summary
DESeq2::summary(res.circ)
create_outputs(dds.circ, res.circ, marker = "condition", out = "circ_rna_DE", nsub = 100,
               padj = argv$fdr, log2FC = argv$log2fc, palette = argv$palette)

# save R image
save.image(file = "DESeq2.RData")
