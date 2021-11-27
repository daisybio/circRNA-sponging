#!/usr/bin/env Rscript

library(SPONGE, biomaRt, argparser, data.table, ggplot2, dplyr)
# create dirs in cwd
dir.create("plots", showWarnings = FALSE)

# TODO: add transpose argument or check if it is necessary
parser <- arg_parser("Argument parser for SPONGE analysis", name = "SPONGE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--circ_rna", help = "Path to circ_rna detection results containing gene ids in tsv")
parser <- add_argument(parser, "--mirna_expr", help = "miRNA expression file in tsv format")
parser <- add_argument(parser, "--organism", help = "Organism given in three letter code")
parser <- add_argument(parser, "--target_scan_symbols", help = "Matrix of target scan symbols provided as tsv")
parser <- add_argument(parser, "--fdr", help = "FDR rate for ceRNA networks")
# add all target scan symbol options to be included -> will generate final target scan symbols
parser <- add_argument(parser, "--miRTarBase_loc", help = "MiRTarBase data location in csv format")
parser <- add_argument(parser, "--miranda_data", help = "Miranda output data location in tsv format")
parser <- add_argument(parser, "--TargetScan_data", help = "TargetScan data location")
parser <- add_argument(parser, "--lncBase_data", help = "LncBase data location")
# FLAGS
parser <- add_argument(parser, "-circ_annotated", help = "CircRNA file is annotated with circBaseId or default CircExplorer output", flag = T)

argv <- parse_args(parser)

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

det_strand <- function(x) {
  if (x == "+") {
    return("1")
  }
  if (x == "-") {
    return("-1")
  }
  if (x == "1") {
    return("+")
  }
  if (x == "-1") {
    return("-")
  }
}

notset <- function(x) {
  if (x == "null" || (x)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

split_encoding <- function(coded_vector) {
  d <- c()
  for (i in seq_along(coded_vector)) {
    split1 <- strsplit(coded_vector[i], "_")
    strand <- det_strand(split1[[1]][2])
    split2 <- strsplit(split1[[1]][1], ":")
    chr <- substr(split2[[1]][1], 4, 5)
    split3 <- strsplit(split2[[1]][2], "-")
    start <- split3[[1]][1]
    end <- split3[[1]][2]
    d <- c(d, paste(c(chr, start, end, strand), collapse = ":"))
  }
  return(d)
}

# get gene names for chromosomal regions in miranda output
annotate_miranda <- function(miRTarBase_loc, ensembl_mart) {
  df = data.frame(read.table(miranda_loc, sep = "\t", header = T))
  targets = split_encoding(df$Target)
  
  print("Start mapping miranda ouptut to gene symbols")
  for (i in seq_along(targets)) {
    gene.Symbol <- getBM(attributes = "hgnc_symbol",
                         filters = c("chromosomal_region"),
                         values=targets[i],
                         mart=ensembl_mart)
    df[i, "gene_symbol"] <- gene.Symbol
  }
  print("Finished mapping miranda ouptut")
  return(df)
}

# use pipeline outputs to create target scan symbols
create_target_scan_symbols <- function(miRTarBase, miranda, TargetScan, org_name, ensembl_mart) {
  # process miRTarBase
  targets <- 0
  if (!notset(miRTarBase)) {
    miRTarBase_targets <- data.frame(read.csv(miRTarBase, header = T))
    # filter by organism
    miRTarBase_targets <- miRTarBase_targets[miRTarBase_targets$Species..miRNA. == org_name,]
    # create target scan counts matrix
    targets <- miRTarBase_targets[, c("miRNA", "Target.Gene")]
  }
  # process targets from miranda
  if (!notset(miranda)) {
    miranda_file <- data.frame(read.table(argv$miranda_data, header = T, sep = "\t"))
    annotate_miranda(miranda_file, ensembl_mart = mart)
    miranda_data <- miranda_file[, c("miRNA", "gene_symbol")]
    colnames(miranda_data) <- c("miRNA", "Target.Gene")
    if (targets == 0) {
      targets <- miranda_data
    } else {
      targets <- merge(targets, miranda_data, by = c("miRNA", "Target.Gene"), all = T)
    }
  }
  # TODO: process TargetScan data
  if (!notset(TargetScan)) {
    target_scan_file <- data.frame(read.table(argv$TargetScan_data, header = T, sep = "\t"))
    # make gene names simple
  }
  # remove NA rows
  targets[complete.cases(targets), ]
  # return contingency table
  return(as.data.frame.matrix(table(target_scan_symbols$miRNA, target_scan_symbols$Target.Gene)))
}

# check required inputs
if (notset(argv$gene_expr) || notset(argv$mirna_expr) || notset(argv$organism)) {
  stop("One or more mandatory arguments are not given")
}
# choose organism
org_data <- org_codes[argv$organism][[1]]
# set up mart
mart <- 0
not_done=TRUE
while(not_done)
{
  tryCatch({
    mart <- useDataset(org_data[2], useMart("ensembl"))
    not_done=FALSE
  }, warning = function(w) {
    print("WARNING SECTION")
    print(w)
  }, error = function(e) {
    print("ERROR SECTION")
    print(e)
  }, finally = {
  })
}

# SET TARGET SCAN SYMBOLS
target_scan_symbols_counts <- 0
# use given target scan symbols for SPONGE
if (!notset(argv$target_scan_symbols)) {
  target_scan_symbols_counts <- data.frame(read.table(file = argv$target_scan_symbols, sep = "\t", header = T))
} else {
  # check for at least one given data set
  if (notset(argv$miRTarBase_loc) && notset(argv$miranda_data) && notset(argv$TargetScan_data)) {
    stop("At least one target scan data source has to be provided")
  }
  # use data from pipeline and build target scan symbols
  target_scan_symbols_counts <- create_target_scan_symbols(miRTarBase = argv$miRTarBase_loc,
                                                           miranda = argv$miranda_data,
                                                           TargetScan = argv$TargetScan_data,
                                                           lncBase = argv$TlncBase_data,
                                                           org_name = org_data[1],
                                                           ensembl_mart = mart)
}
# SET GENE EXPRESSION
gene_expr <- t(as.data.frame(read.table(file = argv$gene_expr, header = TRUE, sep = "\t")))
# READ CIRC_RNA EXPRESSION AND COMBINE THEM
circ_rna_expression <- as.data.frame(read.table(file = argv$circ_rna, header = T, sep = "\t"))
circ_rna_table <- 0
# use db_annotation
if (argv$circ_annotated) {
  circ_filtered <- circ_rna_expression[, c("circRNA.ID", "ensembl_gene_id")]
  circ_filtered <- circ_filtered[complete.cases(circ_filtered),]
  circ_rna_table <- as.data.frame.matrix(table(circ_filtered))
} else {
  compact_raw <- data.table(circRNA.ID=paste0(circ_rna_expression$chr, ":", circ_rna_expression$start, "-", circ_rna_expression$stop,"_", circ_rna_expression$strand), ensembl_gene_id = circ_rna_expression$ensembl_gene_id)
  circ_rna_table <- as.data.frame.matrix(table(compact_raw))
}
# combine expressions
gene_expr <- bind_rows(gene_expr, circ_rna_table)
# SET MIRNA EXPRESSION
mi_rna_expr <- t(as.data.frame(read.table(file = argv$mirna_expr, header = TRUE, sep = "\t")))
# Make genes simple -> ENSG0000001.1 -> ENSG00000001
genes <- gene_expr$X
simple_genes <- c()
for (i in seq_along(genes)) {
  simple_genes[i] = sub("\\..*", "", genes[i])
}
gene_expr$X <- simple_genes
genes <- simple_genes
# get gene names for ENS ids
gene_names <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "hgnc_symbol"), values = genes, mart = mart)
# converting gene ids and reformatting file
gene_expr <- merge(gene_expr, gene_names, by.x = "X", by.y = "ensembl_gene_id")
gene_expr$X <- NULL
colidx <- grep("hgnc_symbol", names(gene_expr))
gene_expr <- gene_expr[, c(colidx, (1:ncol(gene_expr))[-colidx])]
gene_expr <- t(gene_expr)
gene_expr[complete.cases(gene_expr), ]

# ----------------------------- SPONGE -----------------------------
# (A) gene-miRNA interactions
genes_miRNA_candidates <- SPONGE::sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr,
  mir_expr = mi_rna_expr,
  mir_predicted_targets = target_scan_symbols_counts)
# (B) ceRNA interactions
ceRNA_interactions <- SPONGE::sponge(gene_expr = gene_expr,
                             mir_expr = mi_rna_expr,
                             mir_interactions = genes_miRNA_candidates)
# (C) Null-model-based p-value computation
mscor_null_model <- sponge_build_null_model(number_of_datasets = 100, number_of_samples = nrow(gene_expr))
# simulation plot
sim_plot <- sponge_plot_simulation_results(mscor_null_model)
png("plots/simulation.png")
plot(sim_plot)
# ceRNA interaction signs
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, 
                                                   null_model = mscor_null_model)
# (D) ceRNA interaction network
fdr <- as.double(argv$fdr)
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
ceRNA_network_plot <- sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)
png("plots/ceRNA_network.png")
plot(ceRNA_network_plot)
# NETWORK ANALYSIS
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
# different weights for ceRNA interactions
ceRNA_interactions_fdr_weight <- ceRNA_interactions_fdr
ceRNA_interactions_fdr_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
# plot top n samples
n = 3
top_network_plot <- sponge_plot_network_centralities(weighted_network_centralities, top = n)
png("plots/top_ceRNA_network.png")
plot(top_network_plot)
# save R objects
save.image(file = "sponge.RData")
