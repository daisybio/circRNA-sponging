#!/usr/bin/env Rscript

# majority vote
library(BiocGenerics)
library(data.table)
library(stringr)
library(RColorBrewer)
library(VennDiagram)
library(argparser)

# parse arguments
args = commandArgs(trailingOnly = T)
parser <- arg_parser("Majority vote parser", name = "majority")
parser <- add_argument(parser, "--miranda_data", help = "miRanda default output in tsv", default = "null")
parser <- add_argument(parser, "--tarpmir_data", help = "default tarpmir output file in tsv", default = "null")
parser <- add_argument(parser, "--pita_data", help = "Default PITA output", default = "null")
parser <- add_argument(parser, "--majority_matcher", help = "Majority match setting, choose between (start, end, complete)", default = "end")

argv <- parse_args(parser, argv = args)

majority_vote <- function(miranda, tarpmir, pita, match, out = "./") {
  data <- list(miranda, tarpmir, pita)
  # check files
  if (length(Filter(file.exists, data)) != 3) {
    print("not all three predictions given, no majority vote available")
    return(NULL)
  }
  print("calculating majority vote of circRNA-miRNA binding sites")

  print("process miRanda targets")
  miranda.bp <- data.frame(read.table(miranda, header = T, sep = "\t"))
  miranda.bp[,c("start", "end")] <- str_split_fixed(miranda.bp$Subject.Al.Start.End., " ", 2)

  print("process TarPmiR targets")
  tarpmir_data <- fread(tarpmir, header = F, sep = "\t", check.names = F, stringsAsFactors = F, select = 1:3, fill = T)
  tarpmir_data <- cbind(tarpmir_data, str_split_fixed(tarpmir_data$V3, ",", 2))
  colnames(tarpmir_data) <- c("miRNA", "mRNA", "pos", "start", "end")

  print("process PITA targets")
  pita_data <- fread(pita, header = T, sep = "\t", check.names = F, stringsAsFactors = F, select = 1:4, fill = T)
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

  # calculate proportions
  majority.vote <- c(miranda.keys, tarpmir.keys, pita.keys)
  majority.vote <- table(majority.vote)
  # apply majority vote: only use binding sites where 2/3 approve
  majority.vote <- names(majority.vote[majority.vote != 1])
  
  # build table
  splitter <- ifelse(match=="complete", 4, 3)
  majority.vote <- str_split_fixed(majority.vote, "\\|", splitter)
  majority.vote <- data.frame(table(majority.vote[,2], majority.vote[,1]))
  # output venn diagram
  myCol <- brewer.pal(3, "Pastel2")
  venn.diagram(x = list(miranda.keys, tarpmir.keys, pita.keys),
               category.names = c("miRanda", "TarPmiR", "PITA"),
               filename = file.path(out, "binding_sites.png"), output = T,
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
# write majority vote to disk
write.table(majority_vote(miranda = argv$miranda_data,
              tarpmir = argv$tarpmir_data,
              pita = argv$pita_data,
              match = "end"), file = gzfile("majority.tsv.gz"), sep = "\t", row.names = T, quote = F)
