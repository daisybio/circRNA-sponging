#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pacman::p_load(char=args)
errors <- args[!args %in% intersect(args, rownames(installed.packages()))]
if (length(errors)>0) stop("Could not properly install package(s): ", paste(errors, collapse = ","))