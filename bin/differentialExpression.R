#!/usr/bin/env Rscript

library(ggplot2)
library(ensembldb)
library(tximport)
library(pheatmap)
library(DESeq2)

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
  # group
  df <- as.data.frame(colData(d))
  df <- df[order(as.character(df$condition)),]
  d@colData@listData <- as.list(df)
  d@colData@rownames <- rownames(df)
  df <- df[,c("sample", "condition")]
  # normalization
  ntd <- DESeq2::normTransform(d)
  select <- order(rowMeans(counts(d,normalized=TRUE)),
                  decreasing=TRUE)[1:20]
  # filter
  top_genes <- head(order(rowVars(assay(deseq_vst)), decreasing = T), 20)
  match <- assay(deseq_vst)[top_genes,]
  # set output file loc
  heatmap_name <- paste(out, "HMAP", sep = "_")
  # plot heatmap
  pheatmap::pheatmap(match, cluster_rows=T, show_rownames=T,
           cluster_cols=T, annotation_col=df,
           filename = file.path(out, paste(heatmap_name, "png", sep = ".")),
           height = 15, width = 25, legend = F)
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
tx2gene <- ensembldb::transcripts(tx, return.type = "data.frame", columns = c("gene_name", "gene_id"))
failures <- rownames(tx2gene[is.na(tx2gene$gene_name),])
tx2gene[failures, "gene_name"] <- tx2gene[failures, "gene_id"]
tx2gene <- tx2gene[, c("tx_id", "gene_name")]
# txi object
txi <- tximport::tximport(quant.files, type="salmon", tx2gene=tx2gene, ignoreTxVersion = T)
# write total gene expression over samples to file
counts <- as.data.frame(txi$counts)
colnames(counts) <- samplesheet$sample
write.table(counts, file = "gene_expression.tsv", sep = "\t")
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
create_outputs(d = dds, results = res, marker = "condition", out = "total_rna")
# circRNA only
circ_RNAs <- read.table(file = args[4], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
dds_filtered <- dds
filtered_dds <- dds_filtered@assays@data@listData
for(i in seq_along(filtered_dds)) { 
  filtered <- filtered_dds[i][[1]]
  filtered_dds[i][[1]] <- filtered[row.names(filtered) %in% ens_ids,]
}
dds_filtered@assays@data@listData <- filtered_dds
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds_filtered, results = filtered_res, marker = "condition", out = "circ_rna")
# save R image
save.image(file = "DESeq2.RData")
