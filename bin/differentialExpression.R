#!/usr/bin/env Rscript

library("Rsubread", "DESeq2", "reshape2", "ggplot")

args = commandArgs(trailingOnly = TRUE)

# read gene expression counts
countsData <- data.frame(read.table(file = args[1], header = T, sep = "\t"))

# metadata
metaData <- read.table(file = args[2], sep = "\t", header = TRUE)

# create ouptut data and plots
create_outputs <- function(d, results, marker, file_name) {
  # write data to disk
  write.table(cbind(ENS_ID=rownames(results), results), file = paste(file_name, "tsv", sep = "."), quote = FALSE, sep = "\t", col.names = NA)
  # PCA
  vsdata <- DESeq2::vst(d, blind = FALSE)
  PCA_plot <- DESeq2::plotPCA(vsdata, intgroup = marker)
  pca_name <- paste(file_name, "pca", sep = "_")
  png(filename = paste(pca_name, "png", sep = "."))
  plot(PCA_plot)
  # MA
  ma_plot <- DESeq2::plotMA(res)
  ma_name <- paste(file_name, "MA", sep = "_")
  png(filename = paste(ma_name, "png", sep = "."))
  plot(ma_plot)
  # HEAT MAP
  # variance stabilizing transformation
  deseq_vst <- DESeq2::vst(d)
  # convert to data frame
  deseq_vst <- assay(deseq_vst)
  deseq_vst <- as.data.frame(deseq_vst)
  deseq_vst$Gene <- rownames(deseq_vst)
  # filter for significantly differntiated genes log2FolChnage > 3
  deseq_res_df <- as.data.frame(results)
  signif_genes <- rownames(deseq_res_df[deseq_res_df$padj <= .05 &
  abs(deseq_res_df$log2FoldChange) > 3,])
  deseq_vst <- deseq_vst[deseq_vst$Gene %in% signif_genes,]
  # make heatmap
  heatmap <- ggplot(deseq2VST, aes(x=variable, y=Gene, fill=value)) + 
  geom_raster() + scale_fill_viridis(trans="sqrt") + 
  theme(axis.text.x=element_text(angle=65, hjust=1), 
  axis.text.y=element_blank(), axis.ticks.y=element_blank())
  # set output file loc
  heatmap_name <- paste(file_name, "HMAP", sep = "_")
  png(filename = paste(heatmap_name, "png", sep = "."))
  plot(heatmap)
  # close device
  dev.off
}

# dds object
dds <- DESeq2::DESeqDataSetFromMatrix(countData = countsData,
                                      colData = metaData,
                                      design = ~condition)

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
output_loc <- "/nfs/home/students/l.schwartz/test/total_rna"
create_outputs(d = dds, results = res, marker = "condition", file_name = output_loc)
# circRNA only
circ_RNAs <- read.table(file = args[3], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds, results = filtered_res, marker = "condition", file_name = "circ_rna_only")
