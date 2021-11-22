#!/usr/bin/env Rscript

library("reshape2", "ggplot")
# HUMAN ONLY
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

suppressWarnings(BiocManager::install(c("DESeq2", "ensembldb", "tximport")))
library(DESeq2, ensembldb, tximport)

args = commandArgs(trailingOnly = TRUE)

# create ouptut data and plots
create_outputs <- function(d, results, marker, out) {
  # create dirs in cwd
  dir.create(out, showWarnings = FALSE)
  # write data to disk
  write.table(cbind(ENS_ID=rownames(results), results), file = paste(out, "tsv", sep = "."), quote = FALSE, sep = "\t", col.names = NA)
  # PCA
  vsdata <- DESeq2::vst(d, blind = FALSE)
  PCA_plot <- DESeq2::plotPCA(vsdata, intgroup = marker)
  pca_name <- paste(out, "pca", sep = "_")
  png(filename = paste(pca_name, "png", sep = "."))
  plot(PCA_plot)
  # MA
  ma_plot <- DESeq2::plotMA(res)
  ma_name <- paste(out, "MA", sep = "_")
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
  heatmap_name <- paste(out, "HMAP", sep = "_")
  png(filename = paste(heatmap_name, "png", sep = "."))
  plot(heatmap)
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
tx2gene <- transcripts(tx, return.type="DataFrame")
tx2gene <- tx2gene[, c("tx_name", "gene_id")]
# txi object
txi <- tximport::tximport(quant.files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T)
# write total gene expression over samples to file
counts <- txi$counts
colnames(counts) <- samplesheet$sample
write.table(counts, file = file.path(output_dir, "gene_expression.tsv"), sep = "\t")
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
output_loc <- args[5]
create_outputs(d = dds, results = res, marker = "condition", out = paste("total_rna", output_loc, sep = "/"))
# circRNA only
circ_RNAs <- read.table(file = args[4], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds, results = filtered_res, marker = "condition", out = paste("circ_rna", output_loc, sep = "/"))
