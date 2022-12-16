#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(SPONGE, argparser, data.table, dplyr, ggplot2, reshape2, stringr, VennDiagram, Biostrings, MetBrewer)

args = commandArgs(trailingOnly = T)

parser <- arg_parser("Argument parser for SPONGE analysis", name = "SPONGE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--circ_expr", help = "circRNA expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--circ_annotation", help = "circRNA expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--mirna_expr", help = "miRNA expression file in tsv format")
parser <- add_argument(parser, "--mir_fasta", help = "miRNA fasta file")
parser <- add_argument(parser, "--meta", help = "Meta data describing the samples")
parser <- add_argument(parser, "--circ_bs", help = "circRNA-miRNA binding sites")
parser <- add_argument(parser, "--target_scan_symbols", help = "Linear RNA-miRNA binding sites")
parser <- add_argument(parser, "--output_dir", help = "Output directory", default = getwd())
parser <- add_argument(parser, "--fdr", help = "FDR rate for ceRNA networks", default = 0.01)
parser <- add_argument(parser, "--pseudocount", help = "Pseudo-counts to use", default = 1e-3)
parser <- add_argument(parser, "--cpus", help = "Number of cores to use for parallelization", default = 25)
parser <- add_argument(parser, "--tpm_map", help = "TPM map of circular and linear transcripts provided by pipeline")
parser <- add_argument(parser, "--total_bindings", help = "Option to pass a compatible complete mRNA-miRNA binding site tsv file including circRNAs", default = NULL)
parser <- add_argument(parser, "--f_test_p", help = "Option to filter with different f test p values", default = 0.05)
parser <- add_argument(parser, "--coef_t", help = "Option to filter with different coefficient thresholds", default = -0.05)

# FLAGS
parser <- add_argument(parser, "--normalize", help = "Normalize given gene expression before analysis", flag = T)
parser <- add_argument(parser, "--tpm", help = "Use TPM instead of counts", flag = T)
parser <- add_argument(parser, "--no_net", help = "Option to filter with elastic net", flag = T)

argv <- parse_args(parser, argv = args)

print("-------------- PARAMETERS ----------------")
print(argv)

out <- argv$output_dir

# create plots in cwd
dir.create(file.path(out, "plots"), showWarnings = F)

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
create_target_scan_symbols <- function(merged_data, majority, total_bindings) {
  if(!is.null(total_bindings)) {
    print("using manually provided total RNA-miRNA bindings")
    return(total_bindings)
  }
  merged_data_targets <- read.table(merged_data, header = T, sep = "\t", check.names = F, stringsAsFactors = F)

  print("Combining linear and circular miRNA-binding sites...")
  targets_data <- list(merged_data_targets, majority)
  # remove null elements from list
  targets_data[sapply(targets_data, is.null)] <- NULL
  # bind targets
  merged.targets <- bind_rows(targets_data)
  # merge
  merged.targets <- data.frame(dcast(merged.targets, Var1 ~ Var2), row.names = 1, check.names = F, stringsAsFactors = F)
  # set NAs to 0
  merged.targets[is.na(merged.targets)] <- 0
  # convert back to matrix
  merged.targets <- as.matrix(merged.targets)
  print("finished")
  # return contingency table
  return(merged.targets)
}

# TODO: proper sub-network with all extensions
# filter ceRNA networks
subnetwork <- function(interactions, pattern){
  subnetwork <- interactions[grepl(pattern, interactions$geneA) | grepl(pattern, interactions$geneB), ]
  not_circ_A <- subnetwork$geneA[!grepl(pattern, subnetwork$geneA)]
  not_circ_B <- subnetwork$geneB[!grepl(pattern, subnetwork$geneB)]
  not_circ <- c(not_circ_A, not_circ_B)
  interactions_w_circ <- interactions[interactions$geneA %in% not_circ | interactions$geneB %in% not_circ,]
  return(merge(subnetwork, interactions_w_circ, all = T))
}

circ.mRNA.subnetwork <- function(interactions, pattern) {
  return(interactions[(grepl(pattern, interactions$geneA) & grepl("EN", interactions$geneB)) | 
                        (grepl(pattern, interactions$geneB) & grepl("EN", interactions$geneA)),])
}

# Normalize
normalize.data <- function(data){
  samples <- colnames(data)
  meta <- data.frame(samples)
  row.names(meta) <- meta$samples 
  data <- as.matrix(data)
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(data + 1), colData = meta, design = ~ 1)
  dds <- DESeq2::estimateSizeFactors(dds)
  return(DESeq2::counts(dds, normalized=T))
}

# PROCESS INPUTS ------------------------------------------------
majority <- NULL
total_bindings <- NULL
# read precomputed total bindings of linear and cirular RNA
if (file.exists(argv$total_bindings)){
  total_bindings <- read.table(argv$total_bindings, sep = "\t", stringsAsFactors = F, check.names = F)
# read majority vote
} else if (file.exists(argv$circ_bs)) {
  majority <- read.table(argv$circ_bs, sep = "\t", check.names = F, header = T)
} else {
  stop("Please provide circular RNA-miRNA binding sites (--circ_bs) or precomputed RNA-miRNA binding sites (--total_bindings)")
}

# SET TARGET SCAN SYMBOLS
target_scan_symbols_counts <- create_target_scan_symbols(merged_data = argv$target_scan_symbols,
                                                         majority = majority,
                                                         total_bindings = total_bindings)
# SET MIRNA EXPRESSION
print("reading miRNA expression...")
mi_rna_expr <- read.table(file = argv$mirna_expr, header = T, check.names = F, row.names = 1)
# SET GENE EXPRESSION
print("reading gene expression...")
gene_expr <- read.table(file = argv$gene_expr, header = T, check.names = F)

if (argv$normalize & !argv$tpm) {
  # normalize expressions if not already done
  print("normalizing gene expression")
  gene_expr <- normalize.data(gene_expr)
}

# READ CIRC_RNA EXPRESSION AND COMBINE THEM
print("adding circRNA expression...")
circ_filtered_raw <- read.table(file = argv$circ_expr, header = T, check.names = F)
print("reading samplesheet...")
meta <- read.csv(argv$meta, sep = "\t", check.names = F)

# filter for expressions only
circ_filtered <- circ_filtered_raw[,meta$sample]
# remove NAs
circ_filtered <- circ_filtered[complete.cases(circ_filtered),]

# combine linear and circular expressions
gene_expr <- rbind(gene_expr, circ_filtered)

# filter for matching samples
print("gene_expr samples:")
dim(gene_expr)
print("miRNA expr samples:")
dim(mi_rna_expr)
print("target scan symbols samples:")
dim(target_scan_symbols_counts)
# shared samples
shared.samples <- intersect(colnames(gene_expr), colnames(mi_rna_expr))
all.samples <- unique(c(colnames(gene_expr), colnames(mi_rna_expr)))
if (length(shared.samples) < length(all.samples)) {
  cat("Discarding samples", all.samples[!all.samples %in% shared.samples], "because they are missing in either gene or miRNA expression\n")
}
gene_expr <- gene_expr[,shared.samples]
mi_rna_expr <- mi_rna_expr[,shared.samples]
dim(gene_expr)

# transform for sponge
gene_expr[is.na(gene_expr)] <- 0
mi_rna_expr[is.na(mi_rna_expr)] <- 0

# use tpms instead of counts
if (argv$tpm) {
  # convert circRNA and linear expression
  print("using TPMs instead of counts for gene and miRNA expression")
  TPM.map <- read.table(argv$tpm_map, header = T, sep = "\t")
  gene_expr <- TPM.map[rownames(gene_expr), colnames(gene_expr)]
  # convert miRNA expression
  mir_fasta <- readDNAStringSet(argv$mir_fasta)
  mi_rna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% names(mir_fasta),]
  lengths <- mir_fasta[names(mir_fasta)%in%rownames(mi_rna_expr)]@ranges@width
  mi_rna_expr <- mi_rna_expr/lengths
  mi_rna_expr <- log2(t(t(mi_rna_expr)*1e6/colSums(mi_rna_expr)) + argv$pseudocount)
  
  # transform for sponge
  gene_expr[is.na(gene_expr)] <- 0
  mi_rna_expr[is.na(mi_rna_expr)] <- 0
}

# annotate circRNAs if possible
if("circBaseID" %in% colnames(circ_filtered_raw)) {
    # get all circBase IDs for row names in the circRNA expression file
    IDs <- rownames(gene_expr)
    IDs <- merge(IDs, circ_filtered_raw[,"circBaseID", drop = F], by = 0, all.x = T)
    # set NAs to None keyword
    IDs[is.na(IDs$y),"y"] <- "None"
    # only change names that are present in annotation
    IDs[IDs[,"y"]!="None","x"] <- IDs[IDs[,"y"]!="None","y"]
    rownames(gene_expr) <- IDs$x
}

# transpose for SPONGE
mi_rna_expr <- as.matrix(t(mi_rna_expr))
gene_expr <- as.matrix(t(gene_expr))
# cast to matrix
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
# number of cores from args
num.of.cores <- argv$cpus

cl <- makeCluster(num.of.cores, outfile=logging.file) 
registerDoParallel(cl)

print("calculating gene-miRNA interactions...")
# (A) gene-miRNA interactions
genes_miRNA_candidates <- SPONGE::sponge_gene_miRNA_interaction_filter(
  gene_expr = gene_expr,
  mir_expr = mi_rna_expr,
  mir_predicted_targets = target_scan_symbols_counts,
  log.level = "INFO",
  elastic.net = !argv$no_net,
  F.test.p.adj.threshold = argv$f_test_p,
  coefficient.threshold = argv$coef_t)

save.image(file = file.path(out, "sponge.RData"))
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
min.interactions <- 5000
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
if (nrow(ceRNA_interactions_fdr)<min.interactions && nrow(ceRNA_interactions_sign)>min.interactions) {
  print("Warning: fdr setting too strict, no significant interactions detected; min of padj is:")
  print(min(ceRNA_interactions_sign$p.adj))
  print("adjusting...")
  fdr <- min(ceRNA_interactions_sign$p.adj) * 1.01
  ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
  while (nrow(ceRNA_interactions_fdr)<min.interactions) {
    fdr <- fdr * 1.01
    ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
  }
  cat("adjusted fdr to :", fdr, "to allow for a minimum", min.interactions, "interactions", "\n")
  cat("current fdr relevant interactions:", nrow(ceRNA_interactions_fdr))
  ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$p.adj),]
}
# cut samples if too many are detected TODO: choose cutoff
cutoff <- 8500
if (nrow(ceRNA_interactions_fdr)>10000){
  print("Warning: fdr setting too loose, generated over 10000 significant hits; adjusting")
  ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$p.adj),]
  ceRNA_interactions_fdr <- head(ceRNA_interactions_fdr, cutoff)
}

dir.create("circRNA/plots", recursive = T)
dir.create("total/plots", recursive = T)

# save R objects
save.image(file = file.path(out, "sponge.RData"))

# GENERAL NETWORK
ceRNA_network_plot <- sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)
ceRNA_network_plot <- sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates, )
ceRNA_network_plot$x$edges$label <- paste("mscor:", round(ceRNA_interactions_fdr$mscor, 2))
visNetwork::visSave(ceRNA_network_plot, file = "total/plots/network.html")

# MOST SIGNIFICANT SPONGES
network_centralities <- sponge_node_centralities(ceRNA_interactions_fdr)
ceRNA_interactions_weight <- ceRNA_interactions_fdr
ceRNA_interactions_weight$weight <- -log10(ceRNA_interactions_fdr$p.adj)
weighted_network_centralities <- sponge_node_centralities(ceRNA_interactions_weight)
weighed_network_plot <- sponge_plot_network_centralities(weighted_network_centralities, top = 3)
png("total/plots/centralities.png")
plot(weighed_network_plot)
dev.off()

# STRONGEST LINEAR
ceRNA_strongest <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$mscor, decreasing = T),]
ceRNA_strongest <- head(ceRNA_strongest, 100)
ceRNA_strongest_plot <- sponge_plot_network(ceRNA_strongest, genes_miRNA_candidates, ) %>%
  visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
ceRNA_strongest_plot$x$edges$label <- paste("mscor:", round(ceRNA_strongest$mscor, 2))
visNetwork::visSave(ceRNA_strongest_plot, file = "total/plots/strongest_network.html")
write.table(ceRNA_strongest, "total/strongest_linear_ceRNAs.tsv", sep = "\t", row.names = F)

# CIRC RNA ONLY
ceRNA_interactions_all_circ <- subnetwork(ceRNA_interactions_fdr, pattern = "c")
# CIRC RNA MRNA ONLY
circ.mRNA.only <- circ.mRNA.subnetwork(ceRNA_interactions_fdr, "c")

# add scores
circRNA_network_plot <- sponge_plot_network(ceRNA_interactions_all_circ, genes_miRNA_candidates, )
circRNA_network_plot$x$edges$label <- paste("mscor:", round(ceRNA_interactions_all_circ$mscor, 2))

visNetwork::visSave(circRNA_network_plot, file = "circRNA/plots/circ_network.html")
# save significant circRNAs that act as ceRNAs
write.table(ceRNA_interactions_all_circ, "circRNA/circRNAs_as_ceRNAs.tsv", sep = "\t", row.names = F)

# NETWORK ANALYSIS CIRC
ceRNA_interactions_circ_weight <- ceRNA_interactions_all_circ
ceRNA_interactions_circ_weight$weight <- -log10(ceRNA_interactions_all_circ$p.val)
weighted_network_centralities_circ <- sponge_node_centralities(ceRNA_interactions_circ_weight)
# plot top n samples
n = 10
# betweeness
top_network_plot_btw <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n, measure = "btw")
png(file = "circRNA/plots/circ_btw.png")
plot(top_network_plot_btw)
dev.off()
# eigenvector
top_network_plot_ev <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n, measure = "ev")
png(file = "circRNA/plots/circ_ev.png")
plot(top_network_plot_ev)
dev.off()
# counts
top_network_plot_c <- sponge_plot_network_centralities(weighted_network_centralities_circ, top = n, measure = "count")
png(file = "circRNA/plots/circ_counts.png")
plot(top_network_plot_c)
dev.off()
stopCluster(cl) # stop cluster
# save R objects
save.image(file = file.path(out, "sponge.RData"))
