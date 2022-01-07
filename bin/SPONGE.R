#!/usr/bin/env Rscript

library(SPONGE)
library(argparser)
library(data.table)
library(dplyr)
library(ggplot2)

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for SPONGE analysis", name = "SPONGE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--circ_rna", help = "Path to circ_rna detection results containing gene ids in tsv")
parser <- add_argument(parser, "--mirna_expr", help = "miRNA expression file in tsv format")
parser <- add_argument(parser, "--organism", help = "Organism given in three letter code")
# add all target scan symbol options to be included -> will generate final target scan symbols
parser <- add_argument(parser, "--output_dir", help = "Output directory", default = getwd())
parser <- add_argument(parser, "--fdr", help = "FDR rate for ceRNA networks", default = 0.01)
parser <- add_argument(parser, "--target_scan_symbols", help = "Contingency matrix of target scan symbols provided as tsv", default = "null")
parser <- add_argument(parser, "--miRTarBase_loc", help = "MiRTarBase contingency data location in csv format", default = "null")
parser <- add_argument(parser, "--miranda_data", help = "miRanda default output in tsv", default = "null")
parser <- add_argument(parser, "--tarpmir_data", help = "default tarpmir output file in tsv", default = "null")
parser <- add_argument(parser, "--TargetScan_data", help = "TargetScan contingency data location in tsv format", default = "null")
parser <- add_argument(parser, "--lncBase_data", help = "LncBase contingency data location in tsv format", default = "null")
parser <- add_argument(parser, "--miRDB_data", help = "miRDB contingency data location in tsv format", default = "null")

argv <- parse_args(parser, argv = args)

print("-------------- PARAMETERS ----------------")
print(argv)

out <- argv$output_dir

# create plots in cwd
dir.create(file.path(out, "plots"), showWarnings = F)

# define organism three letter codes
org_codes <- list("ebv" = c("Epstein Barr virus", ""),
                  "oar" = c("Ovis aries", ""),
                  "hcmv" = c("Human cytomegalovirus", ""),
                  "cel" = c("Caenorhabditis elegans", ""),
                  "cfa" = c("Canis familiaris", ""),
                  "gga" = c("Gallus gallus", ""),
                  "cgr" = c("Cricetulus griseus", ""),
                  "tgu" = c("Taeniopygia guttata", ""),
                  "xla" = c("Yenips laevis", ""),
                  "ola" = c("Oryzias latipes", ""),
                  "dme" = c("Drosophila melanogaster", "dmelanogaster_gene_ensembl"),
                  "bmo" = c("Bombyx mori", ""),
                  "mmu" = c("Mus musculus", "mmusculus_gene_ensembl"),
                  "rno" = c("Rattus norvegicus", "rnorvegicus_gene_ensembl"),
                  "dre" = c("Dano rerio", ""),
                  "hsa" = c("Homo sapiens", "hsapiens_gene_ensembl"),
                  "kshv" = c("Kaposi sarcoma-associated herpesvirus", "Herpes"),
                  "osa" = c("Oryza sativa", ""),
                  "ssc" = c("Sus scrofa", ""),
                  "bta" = c("Bos taurus", "btaurus_gene_ensembl"),
                  "ath" = c("Arabidopsis thaliana", ""),
                  "xtr" = c("Xenopus tropicalis", ""))

# FUNCTIONS ----------------------------------------------------
remove_ext <- function(g) {
  simple_genes <- c()
  for (i in seq_along(g)) {
    simple_genes[i] = sub("\\..*", "", g[i])
  }
  return(simple_genes)
}

det_strand <- function(x) {
  if (x == "+") {
    return("1")
  }
  if (x == "-") {
    return("-1")
  }
  if (x == 1) {
    return("+")
  }
  if (x == -1) {
    return("-")
  }
}

notset <- function(x) {
  if (x == "null" || is.na(x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

adj_fdr <- function(interactions, x){
  tmp = interactions
  while (nrow(tmp) > 10000 && x < 1){
    tmp = interactions[which(interactions$p.adj < x),]
    x = x+0.1
  }
}

split_encoding <- function(coded_vector) {
  split1 <- strsplit(coded_vector, "_")
  strand <- det_strand(split1[[1]][2])
  split2 <- strsplit(split1[[1]][1], ":")
  chr <- substr(split2[[1]][1], 4, 5)
  split3 <- strsplit(split2[[1]][2], "-")
  start <- split3[[1]][1]
  end <- split3[[1]][2]
  return(paste(c(chr, start, end, strand), collapse = ":"))
}

# use pipeline outputs to create target scan symbols
create_target_scan_symbols <- function(merged_data, miRTarBase, miranda, tarpmir, TargetScan, lncBase, miRDB, org_data) {
  print("CREATING TARGET SCAN SYMBOLS")
  # check targets
  start_file <- 0
  data <- list(merged_data, miRTarBase, miranda, TargetScan, lncBase, miRDB)
  
  # check files
  if (length(Filter(file.exists, data)) == 0) {
    stop("Error: At least one target symbol file has to be supplied")
  }
  
  # process given targets
  merged_data_targets <- NULL
  if (file.exists(merged_data)) {
    print("using given targets")
    merged_data_targets <- data.frame(read.table(merged_data, header = T, sep = "\t"))
  }
  # process miRTarBase
  miRTarBase_targets <- NULL
  if (file.exists(miRTarBase)) {
    print("processing miRTarBase targets")
    miRTarBase_targets <- data.frame(read.csv(miRTarBase, header = T))
    # filter by organism
    miRTarBase_targets <- miRTarBase_targets[miRTarBase_targets$Species..miRNA. == org_data[1],]
    # create target scan counts matrix
    miRTarBase_targets <- as.data.frame.matrix(table(miRTarBase_targets$Target.Gene, miRTarBase_targets$miRNA))
  }
  # process targets from miranda
  miranda_targets <- NULL
  if (file.exists(miranda)) {
    print("processing miranda targets")
    miranda.bp <- data.frame(read.table(miranda, header = T, sep = "\t"))
    miranda_targets <- as.data.frame.matrix(table(miranda.bp$Target, miranda.bp$miRNA))
  }
  # process tarpmir data
  tarpmir_targets <- NULL
  if (file.exists(tarpmir)) {
    print("processing TarPmiR data")
    tarpmir_data <- read.table(tarpmir, header = F, sep = "\t")
    tarpmir_targets <- as.data.frame.matrix(table(tarpmir_data$V2, tarpmir_data$V1))
  }
  # process TargetScan data
  target_scan_targets <- NULL
  if (file.exists(TargetScan)) {
    print("processing TargetScan data")
    target_scan_data <- data.frame(read.table(TargetScan, header = T, sep = "\t"))
    target_scan_targets <- as.data.frame.matrix(table(target_scan_data$Gene.Symbol, target_scan_data$miR.Family))
  }
  # process lncBase data
  lncBase_targets <- NULL
  if (file.exists(lncBase)) {
    print("processing lncBase data")
    lncBase_targets <- as.data.frame(data.frame(read.table(lncBase, sep = "\t", header = T)))
  }
  
  # MERGE DATA
  merged.targets <- NULL
  targets_data <- list(merged_data_targets, miRTarBase_targets, miranda_targets, tarpmir_targets, target_scan_targets, lncBase_targets)
  for (target in targets_data) {
    # append data if present
    if (!is.null(target)) {
      colnames(target) <- sapply(gsub("\\.", "-", colnames(target)), "[", 1)
      # first data set
      if (is.null(merged.targets)) {
        merged.targets <- target
      } else {
        # merge by rowname
        merged.targets <- merge(merged.targets, target, by = 0, all = T)
        # redo columns
        colnames(merged.targets) <- sapply(strsplit(colnames(merged.targets), "\\."), "[", 1)
        # reformat rows
        rownames(merged.targets) <- merged.targets$Row
        # drop Row
        merged.targets$Row <- NULL
        # remove NAs
        merged.targets[is.na(merged.targets)] <- 0
        # combine tables
        merged.targets <- do.call(cbind,lapply(split(seq_len(ncol(merged.targets)),names(merged.targets)),function(x) rowSums(merged.targets[x])))
      }
    }
  }
  # return contingency table
  return(merged.targets)
}

# PROCESS INPUTS ------------------------------------------------
# check required inputs
if (notset(argv$gene_expr) || notset(argv$mirna_expr) || notset(argv$organism)) {
  stop("One or more mandatory arguments are not given")
}

# choose organism
org_data <- org_codes[argv$organism][[1]]

# SET TARGET SCAN SYMBOLS
target_scan_symbols_counts <- create_target_scan_symbols(merged_data = argv$target_scan_symbols,
                                                        miRTarBase = argv$miRTarBase_loc,
                                                        miranda = argv$miranda_data,
                                                        TargetScan = argv$TargetScan_data,
                                                        lncBase = argv$lncBase_data,
                                                        miRDB = argv$miRDB_data,
                                                        org_data = org_data)
# SET MIRNA EXPRESSION
print("reading miRNA expression...")
mi_rna_expr <- t(data.frame(read.table(file = argv$mirna_expr, header = T, sep = "\t"), row.names = 1))
# SET GENE EXPRESSION
print("reading gene expression...")
gene_expr <- as.data.frame(t(read.table(file = argv$gene_expr, header = T, sep = "\t")))

# READ CIRC_RNA EXPRESSION AND COMBINE THEM
print("adding circRNA expression...")
circ_RNAs <- as.data.frame(read.table(file = argv$circ_rna, header = T, sep = "\t"))
# use given annotation if possible
if ("circBaseID" %in% colnames(circ_RNAs)) {
  circ_RNA_annotation <- data.table(circRNA.ID=circ_RNAs$circBaseID)
  # cut table and annotate rownames
  circ_filtered <- circ_RNAs[,-c(1:8)]
} else {
  circ_RNA_annotation <- data.table(circRNA.ID=paste0(circ_RNAs$chr, ":", circ_RNAs$start, ":", circ_RNAs$stop, ":", circ_RNAs$strand))
  # cut table and annotate row names
  circ_filtered <- circ_RNAs[,-c(1:7)]
}
rownames(circ_filtered) <- circ_RNA_annotation$circRNA.ID

circ_filtered <- circ_filtered[complete.cases(circ_filtered),]
circ_filtered <- as.data.frame(t(circ_filtered))

gene_expr <- cbind(gene_expr, circ_filtered)

# filter for matching samples
print("gene_expr samples:")
dim(gene_expr)
print("miRNA expr samples:")
dim(mi_rna_expr)
print("target scan symbols samples:")
dim(target_scan_symbols_counts)
gene_expr <- gene_expr[rownames(mi_rna_expr),]
print("using gene_expr samples:")
dim(gene_expr)
# transform for sponge
gene_expr[is.na(gene_expr)] <- 0
gene_expr <- as.matrix(gene_expr)
mi_rna_expr[is.na(mi_rna_expr)] <- 0
mi_rna_expr <- as.matrix(mi_rna_expr)
target_scan_symbols_counts <- as.matrix(target_scan_symbols_counts)

print("Gene expression:")
print(gene_expr[1:5, 1:5])
print("miRNA expression:")
print(mi_rna_expr[1:5, 1:5])
print("target scan symbols:")
print(target_scan_symbols_counts[1:5, 1:5])

# ----------------------------- SPONGE -----------------------------
# SET UP CLUSTER
library(doParallel)
library(foreach)

logging.file <- ".sponge.log"

num.of.cores <- 20

cl <- makeCluster(num.of.cores, outfile=logging.file) 
registerDoParallel(cl)

print("calculating gene-miRNA interactions...")
# (A) gene-miRNA interactions
genes_miRNA_candidates <- SPONGE::sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr,
  mir_expr = mi_rna_expr,
  mir_predicted_targets = target_scan_symbols_counts,
  log.level = "INFO")
print("calculating ceRNA interactions...")
# (B) ceRNA interactions
ceRNA_interactions <- SPONGE::sponge(gene_expr = gene_expr,
                             mir_expr = mi_rna_expr,
                             mir_interactions = genes_miRNA_candidates,
                             log.level = "INFO")

save.image(file = file.path(out, "sponge.RData"))

print("building null model...")
# (C) Null-model-based p-value computation
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, number_of_samples = nrow(gene_expr))
# simulation plot
sim_plot <- sponge_plot_simulation_results(mscor_null_model)
png(file.path(out, "plots/simulation.png"))
plot(sim_plot)
# ceRNA interaction signs
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, 
                                                   null_model = mscor_null_model)
print("building ceRNA network...")
# (D) ceRNA interaction network
fdr <- as.double(argv$fdr)
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
if (nrow(ceRNA_interactions_fdr)==0) {
  print("Warning: fdr setting too strict, no significant interactions detected; min of padj is:")
  print(min(ceRNA_interactions_sign$p.adj))
  print("using pvalue and selecting top 100 samples")
  ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.val < fdr),]
  ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$p.val),]
  ceRNA_interactions_fdr <- head(ceRNA_interactions_fdr, 1000)
}
if (nrow(ceRNA_interactions_fdr)>10000){
  print("Warning: fdr setting too loose, generated over 10000 significant hits; adjusting")
  ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$p.adj),]
  ceRNA_interactions_fdr <- head(ceRNA_interactions_fdr, 5000)
}
# GENERAL NETWORK
ceRNA_network_plot <- sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)
ceRNA_network_plot <- sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
ceRNA_network_plot$x$edges$label <- paste("mscor:", round(ceRNA_interactions_fdr$mscor, 2))
visNetwork::visSave(ceRNA_network_plot, file = "plots/network.html")
# MOST SIGNIFICANT SPONGES
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
ceRNA_interactions_weight <- ceRNA_interactions_fdr
ceRNA_interactions_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_weight)
weighed_network_plot <- sponge_plot_network_centralities(weighted_network_centralities, top = 3)

# STRONGEST LINEAR
ceRNA_strongest <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$mscor, decreasing = T),]
ceRNA_strongest <- head(ceRNA_strongest, 100)
ceRNA_strongest_plot <- sponge_plot_network(ceRNA_strongest, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
ceRNA_strongest_plot$x$edges$label <- paste("mscor:", round(ceRNA_strongest$mscor, 2))

# CIRC RNA ONLY
ceRNA_interactions_circ <- ceRNA_interactions_fdr[grepl("circ", ceRNA_interactions_fdr$geneA) | grepl("circ", ceRNA_interactions_fdr$geneB), ]
not_circ_A <- ceRNA_interactions_circ$geneA[!grepl("circ", ceRNA_interactions_circ$geneA)]
not_circ_B <- ceRNA_interactions_circ$geneB[!grepl("circ", ceRNA_interactions_circ$geneB)]
not_circ <- c(not_circ_A, not_circ_B)
ceRNA_interactions_w_circ <- 
interactions_w_circ <- ceRNA_interactions_fdr[ceRNA_interactions_fdr$geneA %in% not_circ | ceRNA_interactions_fdr$geneB %in% not_circ,]
ceRNA_interactions_all_circ <- merge(ceRNA_interactions_circ, interactions_w_circ, all = T)

# add scores
circRNA_network_plot <- sponge_plot_network(ceRNA_interactions_all_circ, genes_miRNA_candidates, ) %>%
                        visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
circRNA_network_plot$x$edges$label <- paste("mscor:", round(ceRNA_interactions_all_circ$mscor, 2))

visNetwork::visSave(circRNA_network_plot, file = "plots/circ_network.html")

# NETWORK ANALYSIS
ceRNA_interactions_circ_weight <- ceRNA_interactions_circ
ceRNA_interactions_circ_weight$weight <- -log10(ceRNA_interactions_circ$p.val)
weighted_network_centralities_circ <- sponge_node_centralities(ceRNA_interactions_circ_weight)
# plot top n samples
n = 3
top_network_plot <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n)
png(file = "plots/circ_top_network.png")
plot(top_network_plot)

dev.off()
stopCluster(cl) # stop cluster
# save R objects
save.image(file = file.path(out, "sponge.RData"))
