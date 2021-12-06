#!/usr/bin/env Rscript

library(SPONGE)
library(biomaRt)
library(argparser)
library(data.table)
# library(ggplot2)
args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for SPONGE analysis", name = "SPONGE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--circ_rna", help = "Path to circ_rna detection results containing gene ids in tsv")
parser <- add_argument(parser, "--mirna_expr", help = "miRNA expression file in tsv format")
parser <- add_argument(parser, "--organism", help = "Organism given in three letter code")
# add all target scan symbol options to be included -> will generate final target scan symbols
parser <- add_argument(parser, "--output_dir", help = "Output directory", default = getwd())
parser <- add_argument(parser, "--fdr", help = "FDR rate for ceRNA networks", default = 0.01)
parser <- add_argument(parser, "--target_scan_symbols", help = "Matrix of target scan symbols provided as tsv", default = "null")
parser <- add_argument(parser, "--miRTarBase_loc", help = "MiRTarBase data location in csv format", default = "null")
parser <- add_argument(parser, "--miranda_data", help = "Miranda output data location in tsv format", default = "null")
parser <- add_argument(parser, "--TargetScan_data", help = "TargetScan data location in tsv format", default = "null")
parser <- add_argument(parser, "--lncBase_data", help = "LncBase data location in tsv format", default = "null")
parser <- add_argument(parser, "--miRDB_data", help = "miRDB data location in tsv format", default = "null")
parser <- add_argument(parser, "--gz_targets", help = "Target scan symbols file is compressed", flag = T)
parser <- add_argument(parser, "--circ_annotation", help = "Path to circRNA annotation file containing circBaseIDs and genomic position of circRNAs in tsv format", default = "null")

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
  if (x == "1") {
    return("+")
  }
  if (x == "-1") {
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
annotate_miranda <- function(miranda_loc, org_data) {
  # set up mart
  mart <- 0
  print("SETTING UP ENSEMBL MART")
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
create_target_scan_symbols <- function(merged_data, miRTarBase, miranda, TargetScan, lncBase, miRDB, org_data) {
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
    # unzip given file if it is
    if (argv$gz_targets) {
      untar(merged_data, exdir = out)
      split = strsplit(merged_data, "/")[[1]]
      merged_data <- split[length(split)]
      merged_data <- paste0(strsplit(merged_data, "\\.")[[1]][1], ".tsv")
      merged_data <- file.path(out, merged_data)
    }
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
    miRTarBase_targets <- as.data.frame.matrix(table(miRTarBase_targets$miRNA, miRTarBase_targets$Target.Gene))
  }
  # process targets from miranda
  miranda_targets <- NULL
  if (file.exists(miranda)) {
    print("processing miranda targets")
    miranda_data <- annotate_miranda(miranda_loc = miranda_loc, org_data = org_data)
    miranda_targets <- as.data.frame.matrix(table(miranda_data$miRNA, miranda_data$gene_symbol))
  }
  # process TargetScan data
  target_scan_targets <- NULL
  if (file.exists(TargetScan)) {
    print("processing TargetScan data")
    target_scan_data <- data.frame(read.table(TargetScan, header = T, sep = "\t"))
    target_scan_targets <- as.data.frame.matrix(table(target_scan_data$miR.Family, target_scan_data$Gene.Symbol))
  }
  # process lncBase data
  lncBase_targets <- NULL
  if (file.exists(lncBase)) {
    print("processing lncBase data")
    lncBase_targets <- as.data.frame(t(data.frame(read.table(lncBase, sep = "\t", header = T))))
  }
  
  # MERGE DATA
  merged.targets <- NULL
  targets_data <- list(merged_data_targets, miRTarBase_targets, miranda_targets, target_scan_targets, lncBase_targets)
  for (target in targets_data) {
    # append data if present
    if (!is.null(target)) {
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
  merged.targets[is.na(merged.targets)] <- 0
  rownames(merged.targets) <- merged.targets$Gene
  merged.targets$Gene <- NULL
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
mi_rna_expr <- as.data.frame(t(read.table(file = argv$mirna_expr, header = T, sep = "\t")))
# SET GENE EXPRESSION
print("reading gene expression...")
gene_expr <- as.data.frame(t(read.table(file = argv$gene_expr, header = T, sep = "\t")))

# READ CIRC_RNA EXPRESSION AND COMBINE THEM
print("adding circRNA expression...")
circ_rna_expression <- as.data.frame(read.table(file = argv$circ_rna, header = T, sep = "\t"))
circ_filtered <- 0
# use db_annotation
if (file.exists(argv$circ_annotation)) {
  print("using circBase IDs")
  circ_rna_annotation <- as.data.frame(read.table(file = argv$circ_annotation, header = T, sep = "\t"))
  # use circBase ID
  circ_annotation <- circ_rna_annotation[, c("chr", "start", "stop", "strand", "circRNA.ID")]
  circ_filtered <- merge(circ_rna_expression, circ_annotation, by = c("chr", "start", "stop", "strand"), all = T)
} else {
  print("using circRNA position as ID")
  circ_filtered <- data.table(circRNA.ID=paste0(circ_rna_expression$chr, ":", circ_rna_expression$start, "-", circ_rna_expression$stop,"_", circ_rna_expression$strand))
  circ_rna_expression$circRNA.ID <- circ_filtered$circRNA.ID
  circ_filtered <- circ_rna_expression
}
circ_filtered <- circ_filtered[, 7:length(circ_filtered)]
circ_filtered <- circ_filtered[complete.cases(circ_filtered),]
rownames(circ_filtered) <- circ_filtered$circRNA.ID
circ_filtered$circRNA.ID <- NULL
circ_filtered <- as.data.frame(t(circ_filtered))

gene_expr <- cbind(gene_expr, circ_filtered)

# filter for matching samples
print("gene_expr samples:")
nrow(gene_expr)
print("miRNA expr samples:")
nrow(mi_rna_expr)
gene_expr <- gene_expr[rownames(mi_rna_expr),]

# transform for sponge
gene_expr[is.na(gene_expr)] <- 0
mi_rna_expr[is.na(mi_rna_expr)] <- 0

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
png(file.path(out, "plots/simulation.png"))
plot(sim_plot)
# ceRNA interaction signs
ceRNA_interactions_sign <- sponge_compute_p_values(sponge_result = ceRNA_interactions, 
                                                   null_model = mscor_null_model)
# (D) ceRNA interaction network
fdr <- as.double(argv$fdr)
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < fdr),]
ceRNA_network_plot <- sponge_plot_network(ceRNA_interactions_fdr, genes_miRNA_candidates)
png(file.path(out, "plots/ceRNA_network.png"))
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
png(file.path(out, "plots/top_ceRNA_network.png"))
plot(top_network_plot)
# save R objects
save.image(file = file.path(out, "sponge.RData"))
