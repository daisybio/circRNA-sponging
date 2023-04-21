library(rtracklayer)
library(GenomicRanges)
library(biomaRt)

#### CLI ARGS ####
args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for differenial splicing analysis", name = "DE_parser")
parser <- add_argument(parser, "-m", help = "Meta data for expressions in tsv")
parser <- add_argument(parser, "-e", help = "circRNA expression in tsv format")
parser <- add_argument(parser, "-s", help = "Sample directory")
parser <- add_argument(parser, "-g", help = "Path to GTF file")
parser <- add_argument(parser, "-t", help = "Path to linear transcript expression")
parser <- add_argument(parser, "-tpm", help = "Path to TPM map")
parser <- add_argument(parser, "-b", help = "biomart_db")
parser <- add_argument(parser, "-o", help = "Output directory", default = "./")
# parameter settings
parser <- add_argument(parser, "--palette", help = "Palette to use for MetBrewer", default = "Renoir")

argv <- parse_args(parser, argv = args)

meta <- read.csv(argv$m, sep = "\t", check.names = F)
circExpression <- read.csv(argv$e, sep = "\t", check.names = F)
circAnnotation <- circExpression[,-which(colnames(circExpression)%in%meta$sample)]
mart <- useMart("ENSEMBL_MART_ENSEMBL", argv$b)
# read circExplorer2 annotate outputs
circExplorerFiles <- list.files(path = argv$s,
                                pattern = "known.txt", recursive = T, full.names = T)
circExons <- bind_rows(sapply(circExplorerFiles, function(annotateOut) {
  c <- read.csv(annotateOut, sep = "\t", check.names = F, header = F)
  c[,ncol(c)] <- as.character(c[,ncol(c)])
  ranges <- apply(c, 1, function(r) {
    exons <- strsplit(r[[length(r)]], "\\|")[[1]]
    exons <- exons[!grepl("None", exons)]
    data.frame(paste0(exons, ":", r[[6]]), 1:length(exons), 
               gsub(" ", "", paste0(r[[1]], ":", r[[2]], "-", r[[3]], "_", r[[6]])),
               r[[15]])
  })
}))
colnames(circExons) <- c("exon", "exon_number", "transcript_id", "gene_symbol")

circExons <- circExons[!startsWith(circExons$exon, ":"),]
circExons$exon_id <- paste0(circExons$exon, ".", circExons$exon_number)
conv <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
              filters = "external_gene_name",
              mart = mart,
              values = unique(circExons$gene_symbol))
circExons <- merge(circExons, conv, by.x = "gene_symbol", by.y = "external_gene_name")
circExonsRanges <- GRanges(circExons$exon)
circExonsRanges$gene_id <- circExons$ensembl_gene_id
circExonsRanges$transcript_id <- circExons$transcript_id
circExonsRanges$exon_number <- circExons$exon_number
circExonsRanges$exon_id <- circExons$exon_id
circExonsRanges$type <- "exon"
circExonsRanges$source <- "CIRCexplorer2"
circExonsRanges$phase <- NA
circExonsRanges <- unique(circExonsRanges)

gtf <- import(argv$g, format = "gtf")
# convert gene names to ENSGs

conv <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
              filters = "external_gene_name",
              mart = mart,
              values = unique(circAnnotation$gene_symbol))
circAnnotation$transcript_id <- rownames(circAnnotation)
circAnnotation <- merge(circAnnotation, conv, by.x = "gene_symbol", by.y = "external_gene_name")
# add gtf columns and rename existing attributes
colnames(circAnnotation)[grep("gene_id", colnames(circAnnotation))] <- "gene_id"
colnames(circAnnotation)[grep("gene_symbol", colnames(circAnnotation))] <- "gene_name"
circAnnotation$source <- "CIRCexplorer2"
circAnnotation$type <- "transcript"
circAnnotation$phase <- NA
circAnnotation$exon_number <- NA
circAnnotation$exon_id <- NA
circAnnotation$key <- NULL
circAnnotation$circBaseID <- NULL
# convert to ranges
circRanges <- makeGRangesFromDataFrame(circAnnotation, keep.extra.columns = T)
# include exons from circExplorer2 annotate output
# combine gtfs
gtfRanges <- sort(c(gtf, circRanges, circExonsRanges))
# export as gtf
export(gtfRanges, file.path(argv$o, "inclCircRNAs.gtf"))
# import circRNA tpms
transcriptTPMs <- read.csv(argv$t, sep = "\t", check.names = F)
circTPMs <- read.csv(argv$tpm, sep = "\t", check.names = F)
circTPMs <- circTPMs[rownames(circExpression),]
allTPMs <- rbind(transcriptTPMs, circTPMs)
write.table(allTPMs, file.path(argv$o, "allTranscriptsTPMs.tsv"), sep = "\t", quote = F, row.names = T)
