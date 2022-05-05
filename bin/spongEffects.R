#!/usr/bin/env Rscript

install.packages("pacman")
pacman::p_load(SPONGE, doParallel, foreach, dplyr, argparser)

args = commandArgs(trailingOnly = TRUE)
# TODO add fine tuning parameters
parser <- arg_parser("Argument parser for differenial expression analysis", name = "DE_parser")
parser <- add_argument(parser, "--spongeData", help = "SPONGE Rdata containing all results, e.g. gene expression, miRNA expression, sponge centralities etc.")
parser <- add_argument(parser, "--meta", help = "Meta data for samples in tsv")
parser <- add_argument(parser, "--fdr", help = "False discovery rate to use for filtering SPONGE results", default = 0.05)
parser <- add_argument(parser, "--cpus", help = "Number of cores to use for backend", default = 4)

argv <- parse_args(parser, argv = args)

# backend
num.of.cores <- argv$cpus
cl <- makeCluster(num.of.cores) 
registerDoParallel(cl)

# meta data
meta <- read.csv(file = argv$meta, sep = "\t")

# gene expression
train_gene_expr <- gene_expr

# miRNA expression
train_mirna_expr <- mi_rna_expr

# ceRNA interactions
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < argv$fdr),]
# network centralities
network_centralities  <- sponge_node_centralities(ceRNA_interactions_fdr)


# spongEffects filtering
filtered_network_centralities <- filter_ceRNA_network(sponge_effects = train_ceRNA_interactions, 
                                                   Node_Centrality = train_network_centralities,
                                                   add_weighted_centrality=T, 
                                                   mscor.threshold = 0.01, 
                                                   padj.threshold = 0.1)


