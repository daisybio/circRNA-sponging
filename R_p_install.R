#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pacman::p_load(args)
sapply(args, function(pkg) if(!require(pkg, character.only = T)) stop("package ", pkg, " could not be properly installed"))
