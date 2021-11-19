# install necessary packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SPONGE")
args = commandArgs(trailingOnly = TRUE)
# library("SPONGE")
library(biomaRt, argparser)

# TODO: add transpose argument or check if it is necessary
parser <- arg_parser("Argument parser for SPONGE analysis", name = "SPONGE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")
parser <- add_argument(parser, "--mirna_expr", help = "miRNA expression file in tsv format")
parser <- add_argument(parser, "--organism", help = "Organism given in three letter code")
parser <- add_argument(parser, "--target_scan_symbols", help = "Matrix of target scan symbols provided as tsv")
# add all target scan symbol options to be included -> will generate final target scan symbols
parser <- add_argument(parser, "--miRTarBase_loc", help = "MiRTarBase data location in csv format")
parser <- add_argument(parser, "--miranda_data", help = "Miranda output data location in tsv format")
parser <- add_argument(parser, "--TargetScan_data", help = "TargetScan data location")

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
  } else {
    return("-1")
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
  
  gene.Symbols <- getBM(attributes = c("hgnc_symbol", "mirbase_accession", "mirbase_id"),
                        filters = c("chromosomal_region"),
                        values=targets,
                        mart=ensembl_mart)
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
  if (!is.na(miRTarBase)) {
    miRTarBase_targets <- data.frame(read.csv(miRTarBase, header = T))
    # filter by organism
    miRTarBase_targets <- miRTarBase_targets[miRTarBase_targets$Species..miRNA. == org_name,]
    # create target scan counts matrix
    targets <- miRTarBase_targets[, c("miRNA", "Target.Gene")]
  }
  # process targets from miranda
  if (!is.na(miranda)) {
    miranda_file <- data.frame(read.table(args[5], header = T, sep = "\t"))
    annotate_miranda(miranda_file, ensembl_mart = mart)
    miranda_data <- miranda_file[, c("miRNA", "gene_symbol")]
    colnames(miranda_data) <- c("miRNA", "Target.Gene")
    targets <- merge(target_scan_data, miranda_data, by = c("miRNA", "Target.Gene"), all = T)
  }
  # TODO: process TargetScan data
  if (!is.na(TargetScan)) {
    target_scan_file <- data.frame(read.table(TargetScan, header = T, sep = "\t"))
  }
  # remove NA rows
  targets[complete.cases(targets), ]
  # return contingency table
  return(as.data.frame.matrix(table(target_scan_symbols$miRNA, target_scan_symbols$Target.Gene)))
}

# check required inputs
if (is.na(argv$gene_expr) || is.na(argv$mirna_expr) || is.na(argv$organism)) {
  stop("One or more mandatory arguments are not given")
}
# choose organism
org_data <- org_codes[argv$organism][[1]]
# set up mart
mart <- useDataset(org_data[2], useMart("ensembl"))
# SET TARGET SCAN SYMBOLS
target_scan_symbols_counts <- 0
# use given target scan symbols for SPONGE
if (!is.na(argv$target_scan_symbols)) {
  target_scan_symbols_counts <- data.frame(read.table(file = argv$target_scan_symbols, sep = "\t", header = T))
} else {
  # check for at least one given data set
  if (is.na(argv$miRTarBase_loc) && is.na(argv$miranda_data) && is.na(argv$TargetScan_data)) {
    stop("At least one target scan data source has to be provided")
  }
  # use data from pipeline and build target scan symbols
  target_scan_symbols_counts <- create_target_scan_symbols(miRTarBase = argv$miRTarBase_loc,
                                                           miranda = argv$miranda_data,
                                                           TargetScan = argv$TargetScan_data,
                                                           org_name = org_data[1],
                                                           ensembl_mart = mart)
}
# SET GENE EXPRESSION
gene_expr <- as.data.frame(read.table(file = args[1], header = TRUE, sep = "\t"))
# SET MIRNA EXPRESSION
mi_rna_expr <- t(as.data.frame(read.table(file = args[2], header = TRUE, sep = "\t")))
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
# extract all unique miRNAs
miRNAs <- target_scan_symbols$miRNA
miRNAs <- miRNAs[!duplicated(miRNAs)]

# INIT SPONGE


