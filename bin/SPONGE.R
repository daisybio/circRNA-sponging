#!/usr/bin/env Rscript

library(SPONGE)
library(argparser)
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(VennDiagram)
library(Biostrings)

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for SPONGE analysis", name = "SPONGE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--circ_rna", help = "circRNA expression file in tsv format as given by pipeline")
parser <- add_argument(parser, "--mirna_expr", help = "miRNA expression file in tsv format")
parser <- add_argument(parser, "--mir_fasta", help = "miRNA fasta file")
# add all target scan symbol options to be included -> will generate final target scan symbols
parser <- add_argument(parser, "--output_dir", help = "Output directory", default = getwd())
parser <- add_argument(parser, "--fdr", help = "FDR rate for ceRNA networks", default = 0.01)
parser <- add_argument(parser, "--cpus", help = "Number of cores to use for paralellization", default = 25)
parser <- add_argument(parser, "--target_scan_symbols", help = "Contingency matrix of target scan symbols provided as tsv", default = "null")
parser <- add_argument(parser, "--miranda_data", help = "miRanda default output in tsv", default = "null")
parser <- add_argument(parser, "--tarpmir_data", help = "default tarpmir output file in tsv", default = "null")
parser <- add_argument(parser, "--pita_data", help = "Default PITA output", default = "null")
parser <- add_argument(parser, "--majority_matcher", help = "Majority match setting, choose between (start, end, complete)", default = "end")
parser <- add_argument(parser, "--tpm_map", help = "TPM map of circular and linear transcripts provided by pipeline")
# FLAGS
parser <- add_argument(parser, "--normalize", help = "Normalize given gene expression before analysis", flag = T)
parser <- add_argument(parser, "--tpm", help = "Use TPM instead of counts", flag = T)

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

# use only binding sites that 2/3 methods approve of; match criteria: complete(start and end pos), start(start pos), end(end pos)
majority_vote <- function(miranda, tarpmir, pita, match) {
  data <- list(miranda, tarpmir, pita)
  # check files
  if (length(Filter(file.exists, data)) != 3) {
    print("not all three predictions given, no majority vote available")
    return(NULL)
  }
  print("calculating majority vote of circRNA-miRNA binding sites")
  
  print("process miranda targets")
  miranda.bp <- data.frame(read.table(miranda, header = T, sep = "\t"))
  miranda.bp[,c("start", "end")] <- str_split_fixed(miranda.bp$Subject.Al.Start.End., " ", 2)

  print("process tarpmir targets")
  tarpmir_data <- read.table(tarpmir, header = F, sep = "\t")
  tarpmir_data <- tarpmir_data[,c(1:3)]
  tarpmir_data[, c("start", "end")] <- str_split_fixed(tarpmir_data$V3, ",", 2)
  colnames(tarpmir_data) <- c("miRNA", "mRNA", "pos", "start", "end")
  
  print("process PITA targets")
  pita_data <- read.table(pita, header = T, sep = "\t")
  # switch pita output labels
  colnames(pita_data) <- c("UTR", "microRNA", "end", "start")
  
  print("creating keys...")
  cat("using", match, "for matching binding sites\n")
  
  if (match=="complete") {
    miranda.keys <- paste(miranda.bp$miRNA, miranda.bp$Target, miranda.bp$start, miranda.bp$end, sep="|")
    tarpmir.keys <- paste0(tarpmir_data$miRNA, tarpmir_data$mRNA, tarpmir_data$start, tarpmir_data$end, sep="|")
    pita.keys <- paste0(pita_data$microRNA, pita_data$UTR, pita_data$start, pita_data$end, sep="|")
  } else if (match=="start") {
    miranda.keys <- paste(miranda.bp$miRNA, miranda.bp$Target, miranda.bp$start, sep="|")
    tarpmir.keys <- paste(tarpmir_data$miRNA, tarpmir_data$mRNA, tarpmir_data$start, sep="|")
    pita.keys <- paste(pita_data$microRNA, pita_data$UTR, pita_data$start, sep="|")
  } else if (match=="end") {
    miranda.keys <- paste(miranda.bp$miRNA, miranda.bp$Target, miranda.bp$end, sep="|")
    tarpmir.keys <- paste(tarpmir_data$miRNA, tarpmir_data$mRNA, tarpmir_data$end, sep="|")
    pita.keys <- paste(pita_data$microRNA, pita_data$UTR, pita_data$end, sep="|")
  } else {
    stop("Wrong --majority_matcher argument given, use one of 'complete', 'start', 'end'")
  }
  
  # search for intersections
  miranda_x_tarpmir <- intersect(miranda.keys, tarpmir.keys)
  miranda_x_pita <- intersect(miranda.keys, pita.keys)
  tarpmir_x_pita <- intersect(tarpmir.keys, pita.keys)
  # apply majority vote: only use binding sites where 2/3 approve
  majority.vote <- c(miranda_x_tarpmir, miranda_x_pita, tarpmir_x_pita)
  # build table
  splitter <- ifelse(match=="complete", 4, 3)
  majority.vote <- str_split_fixed(majority.vote, "\\|", splitter)
  majority.vote <- as.data.frame.matrix(table(majority.vote[,2], majority.vote[,1]))
  # output VENN diagram
  library(RColorBrewer)
  myCol <- brewer.pal(3, "Pastel2")
  venn.diagram(x = list(miranda.keys, tarpmir.keys, pita.keys), 
               category.names = c("miRanda", "TarPmiR", "PITA"), 
               filename = file.path(out, "plots/binding_sites.png"), output = T,
               imagetype="png",
               height = 800, 
               width = 1200, 
               resolution = 300,
               compression = "lzw",
               
               # Circles
               lwd = 2,
               lty = 'blank',
               fill = myCol,
               
               # Numbers
               cex = .6,
               fontface = "bold",
               fontfamily = "sans",
               
               # Set names
               cat.cex = 0.6,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(-27, 27, 135),
               cat.dist = c(0.055, 0.055, 0.085),
               cat.fontfamily = "sans",
               rotation = 1)
  return(majority.vote)
}

# use pipeline outputs to create target scan symbols
create_target_scan_symbols <- function(merged_data, majority, miranda, tarpmir, pita) {
  print("CREATING TARGET SCAN SYMBOLS")
  print("using given targets")
  merged_data_targets <- data.frame(read.table(merged_data, header = T, sep = "\t"))
  miranda_targets <- NULL
  tarpmir_targets <- NULL
  pita_targets <- NULL
  
  # add individually if no majority vote could be formed
  if (is.null(majority)) {
    # process targets from miranda
    if (file.exists(miranda)) {
      print("processing miranda targets")
      miranda.bp <- data.frame(read.table(miranda, header = T, sep = "\t"))
      miranda_targets <- as.data.frame.matrix(table(miranda.bp$Target, miranda.bp$miRNA))
    }
    # process tarpmir data
    if (file.exists(tarpmir)) {
      print("processing TarPmiR data")
      tarpmir_data <- read.table(tarpmir, header = F, sep = "\t")
      tarpmir_targets <- as.data.frame.matrix(table(tarpmir_data$V2, tarpmir_data$V1))
    }
    # process PITA data
    if (file.exists(pita)) {
      print("processsing PITA data")
      pita_data <- read.table(pita, header = T, sep = "\t")
      pita_targets <- as.data.frame.matrix(table(pita_data$RefSeq, pita_data$microRNA))
    }
  }
  
  # MERGE DATA
  merged.targets <- NULL
  targets_data <- list(merged_data_targets, majority, miranda_targets, tarpmir_targets, pita_targets)
  for (target in targets_data) {
    # append data if present
    if (!is.null(target)) {
      colnames(target) <- gsub("\\.", "-", colnames(target))
      # first data set
      if (is.null(merged.targets)) {
        merged.targets <- target
      } else {
        # merge by row names
        merged.targets <- merge(merged.targets, target, by = 0, all = T)
        # reformat rows
        rownames(merged.targets) <- merged.targets$Row.names
        # drop Row
        merged.targets$Row.names <- NULL
        # remove NAs
        merged.targets[is.na(merged.targets)] <- 0
        # combine tables
        colnames(merged.targets) <- sapply(strsplit(colnames(merged.targets), "\\."), "[", 1)
        merged.targets <- do.call(cbind,lapply(split(seq_len(ncol(merged.targets)),names(merged.targets)),function(x) rowSums(merged.targets[x])))
      }
    }
  }
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
# check required inputs
if (notset(argv$gene_expr) || notset(argv$mirna_expr)) {
  stop("One or more mandatory arguments are not given")
}
majority <- NULL
n.targets <- length(sapply(list(argv$miranda, argv$tarpmir, argv$pita), notset))
if (n.targets == 0) {
  stop("Supply at least one ")
} else if (n.targets == 3) {
  majority <- majority_vote(argv$miranda, argv$tarpmir, argv$pita, argv$majority_matcher)
}

# SET TARGET SCAN SYMBOLS
target_scan_symbols_counts <- create_target_scan_symbols(merged_data = argv$target_scan_symbols,
                                                        majority = majority,
                                                        miranda = argv$miranda,
                                                        tarpmir = argv$tarpmir,
                                                        pita = argv$pita)
# SET MIRNA EXPRESSION
print("reading miRNA expression...")
mi_rna_expr <- data.frame(read.table(file = argv$mirna_expr, header = T, sep = "\t"), row.names = 1)
# SET GENE EXPRESSION
print("reading gene expression...")
gene_expr <- as.data.frame(read.table(file = argv$gene_expr, header = T, sep = "\t"))

# READ CIRC_RNA EXPRESSION AND COMBINE THEM
print("adding circRNA expression...")
circ_RNAs <- as.data.frame(read.table(file = argv$circ_rna, header = T, sep = "\t"))

# use given annotation if possible
annotation <- "circBaseID" %in% colnames(circ_RNAs)
if (annotation) {
  circ_RNA_annotation <- ifelse(circ_RNAs$circBaseID != "None", 
                                circ_RNAs$circBaseID, 
                                paste(circ_RNAs$chr, circ_RNAs$start, circ_RNAs$stop, circ_RNAs$strand, sep = ":"))
  # cut table and annotate rownames
  circ_filtered <- circ_RNAs[,-c(1:8)]
} else {
  circ_RNA_annotation <- paste(circ_RNAs$chr, circ_RNAs$start, circ_RNAs$stop, circ_RNAs$strand, sep = ":")
  # cut table and annotate row names
  circ_filtered <- circ_RNAs[,-c(1:7)]
}
rownames(circ_filtered) <- paste0(circ_RNAs$chr, ":", circ_RNAs$start, "-", circ_RNAs$stop, "_", circ_RNAs$strand)

circ_filtered <- circ_filtered[complete.cases(circ_filtered),]

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
cat("gene expression and miRNA expression share", length(shared.samples), "samples\n")
gene_expr <- gene_expr[,shared.samples]
mi_rna_expr <- mi_rna_expr[,shared.samples]
dim(gene_expr)
# transform for sponge
gene_expr[is.na(gene_expr)] <- 0
mi_rna_expr[is.na(mi_rna_expr)] <- 0

target_scan_symbols_counts <- as.matrix(target_scan_symbols_counts)

# use tpms instead of counts
if (argv$tpm) {
  # convert circRNA and linear expression
  TPM.map <- read.table(argv$tpm_map, header = T, sep = "\t")
  gene_expr <- TPM.map[rownames(gene_expr), colnames(gene_expr)]
  # convert miRNA expression
  mir_fasta <- readDNAStringSet(argv$mir_fasta)
  mi_rna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% names(mir_fasta),]
  lengths <- mir_fasta[names(mir_fasta)%in%rownames(mi_rna_expr)]@ranges@width
  mi_rna_expr <- mi_rna_expr/lengths
  mi_rna_expr <- log2(t(t(mi_rna_expr)*1e6/colSums(mi_rna_expr))+1e-3)
}

save.image(file = file.path(out, "sponge.RData"))

# normalize expressions if not already done
if (argv$normalize) {
  print("normalizing gene expression")
  gene_expr <- normalize.data(gene_expr)
}

# annotate circRNAs
if (annotation){
  rownames(gene_expr)[grepl("c", rownames(gene_expr))] <- circ_RNA_annotation
}

# transpose for SPONGE
mi_rna_expr <- as.matrix(t(mi_rna_expr))
gene_expr <- as.matrix(t(gene_expr))

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
min.interactions <- 5000
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
if (nrow(ceRNA_interactions_fdr)<min.interactions) {
  print("Warning: fdr setting too strict, no significant interactions detected; min of padj is:")
  print(min(ceRNA_interactions_sign$p.adj))
  print("adjusting...")
  ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
  while (nrow(ceRNA_interactions_fdr)<min.interactions) {
    fdr <- fdr * 1.01
    ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
  }
  cat("adjusted fdr to :", fdr, "to allow for a minimum", min.interactions, "interactions", "\n")
  cat("current fdr relevant interactions:", nrow(ceRNA_interactions_fdr))
  ceRNA_interactions_fdr <- ceRNA_interactions_fdr[order(ceRNA_interactions_fdr$p.adj),]
  # ceRNA_interactions_fdr <- head(ceRNA_interactions_fdr, 8000)
}
# cut samples if too many are detected TODO: choose cutoff
cutoff <- 5000
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
