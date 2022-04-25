#!/usr/bin/env Rscript

library(SPONGE)
library(doParallel)
library(foreach)
library(dplyr)

args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for differenial expression analysis", name = "DE_parser")
parser <- add_argument(parser, "--gene_expr", help = "Gene expression file in tsv format as given by DeSeq2")

argv <- parse_args(parser, argv = args)

# backend
num.of.cores <- 4
cl <- makeCluster(num.of.cores) 
registerDoParallel(cl)