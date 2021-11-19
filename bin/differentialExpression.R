#!/usr/bin/env Rscript
# install required packages if they do not already exist
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "Rsubread"))

# start script
library("Rsubread", "DESeq2", "reshape2", "ggplot")
# ARGS: path_to_samples[file] path_to_metadata[file] genome_version[str] gtf[file] is_single_end[boolean][true/false], OPT: path_to_circRNAs[file]

args = commandArgs(trailingOnly = TRUE)

# read sam file(s)
sam_files <- list.files(args[1], pattern = "\\.sam$", full.names = TRUE, recursive = TRUE)

# metadata
metaData <- read.table(file = args[2], sep = "\t", header = TRUE)
# genome version
genome_version <- args[3]
# gtf
gtf_file <- args[4]
# read mode
is_single_end <- function(isT) {
  if (isT == "true") {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
# check for paired reads
isPaired <- is_single_end(args[5])
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
# convert sam to countsObject
countsObject <- Rsubread::featureCounts(sam_files, annot.inbuilt = genome_version, annot.ext = gtf_file,
                            isGTFAnnotationFile = TRUE, isPairedEnd = isPaired)
# prepare deseq2 counts data
countsData <- countsObject[["counts"]]
# save general gene expression of all samples
write.table(countsData, file = paste("gene_expression", "tsv", sep = "."), quote = FALSE, sep = "\t", col.names = NA)
# reformat ids
new_names <- list()
for (name in row.names(countsData)) {
  general_name <- strsplit(name, "\\.")[[1]][1]
  new_names <- c(new_names, general_name)
}
row.names(countsData) <- new_names

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
circ_RNAs <- read.table(file = args[6], sep = "\t", header = TRUE)
ens_ids <- circ_RNAs$ensembl_gene_ID
filtered_res <- res[row.names(res) %in% ens_ids,]
create_outputs(d = dds, results = filtered_res, marker = "condition", file_name = "circ_rna_only")
