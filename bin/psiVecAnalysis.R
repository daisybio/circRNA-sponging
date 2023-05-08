library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(MetBrewer)
library(argparser)


annotatePsi <- function(psivec) {
  psivec[,c("gene_id", "transcript_id")] <- str_split_fixed(rownames(psivec), ";", 2)
  psivec$transcript_type <- ifelse(grepl("c", psivec$transcript_id),
                                   "circRNA", "linearRNA")
  psivec[is.na(psivec)] <- 0
  psivec
}

normalizePerTranscriptType <- function(vec, samples) {
  message("normalizing psi values for mean psi of linear and circular transcripts...")
  means <- vec %>% group_by(transcript_type) %>% summarise_at(vars(samples), mean)
  
  vec[vec$transcript_type=="circRNA",samples] <- vec[vec$transcript_type=="circRNA",samples] /
    as.vector(means[means$transcript_type=="circRNA", samples])
  vec[vec$transcript_type=="linearRNA",samples] <- vec[vec$transcript_type=="linearRNA",samples] / 
    as.vector(means[means$transcript_type=="linearRNA", samples])
  vec %>% group_by(gene_id) %>% mutate(across(samples, function(s) s/sum(s)))
}

sumPerTranscriptType <- function(vec, normalized=T) {
  lab <- ifelse(normalized, "Percent Spliced In (PSI) \n normalized by sample-wise mean of transcript type",
                "Percent Spliced In (PSI)")
  name <- ifelse(normalized, "normalized", "raw")
  # group data by host gene and summarize per transcript type if the gene contains circular transcripts
  rnaTypeRatios <- vec %>%
    group_by(gene_id) %>% filter(length(unique(transcript_type))>1) %>% 
    group_by(gene_id, transcript_type) %>%
    summarise_at(samples, sum)
  # collapse replicates of conditions
  rnaTypeRatiosCollapsed <- as.data.frame(t(rowsum(t(rnaTypeRatios[,samples]), group = meta$condition)/c(table(meta$condition))))
  rnaTypeRatiosCollapsed[,c("gene_id", "transcript_type")] <- rnaTypeRatios[,c("gene_id", "transcript_type")]
  # create boxplots between transcript types
  rnaTypeRatiosCollapsedMelted <- melt(rnaTypeRatiosCollapsed, id.vars = c("gene_id", "transcript_type"))
  rnaTypeRatiosCollapsedMelted$value <- rnaTypeRatiosCollapsedMelted$value*100
  boxP <- ggviolin(rnaTypeRatiosCollapsedMelted, x = "transcript_type", y = "value",
                   fill = "variable", add = "boxplot", add.params = list(fill = "white", width = 0.1)) +
    scale_fill_manual(values = colors) +
    stat_compare_means(method = "t.test", label.x = 1, label.y = 125) +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = "linearRNA", vjust = -2.25) +
    scale_y_continuous(breaks = seq(0, 100, 25)) +
    facet_wrap(~variable) +
    xlab("Transcript Type") +
    ylab(lab) +
    theme(text = element_text(size = 18), strip.text.x = element_text(size = 18),
          legend.title = element_blank())
  ggsave(filename = paste0(name, "_", "violins.png"), plot = boxP, width = 12, height = 8, dpi = 300)
  # save circRNAs with higher sums than linear transcripts of same host gene
  # save circRNAs with higher sums than linear transcripts of same host gene
  dPSIs <- rnaTypeRatios[seq(1, nrow(rnaTypeRatios), 2),samples] -
    rnaTypeRatios[seq(2, nrow(rnaTypeRatios)+1, 2),samples]
  rownames(dPSIs) <- unique(rnaTypeRatios$gene_id)
  signifHostGenes <- dPSIs[rowSums(dPSIs > 0) > 0,]
  vec %>% filter(gene_id%in%rownames(signifHostGenes)) %>% 
    filter(transcript_type=="circRNA") %>%
    write.table(file = paste0(name, "_majorityCircRNAs.tsv"), 
                sep = "\t", quote = F, row.names = F)
}

#### CLI ARGS ####
args <- commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for differenial splicing analysis", name = "DE_parser")
parser <- add_argument(parser, "-m", help = "Meta data for expressions in tsv")
parser <- add_argument(parser, "-i", help = "Suppa2 psiPerIsoform .psi file")
# parameter settings
parser <- add_argument(parser, "--palette", help = "Palette to use for MetBrewer", default = "Renoir")

argv <- parse_args(parser, argv = args)

#### READ INPUT DATA #####
message("reading meta data and diffSplice psivec file...")
meta <- read.csv(argv$m, sep = "\t")
samples <- meta$sample
conditions <- unique(meta$condition)
colors <- setNames(met.brewer(argv$palette, length(conditions)), conditions)
# read diffSplice psivec results
psivec <- read.csv(argv$i, sep = "\t", check.names = F)
#### START SAMPLE-WISE ANALYSIS ####
colnames(psivec) <- samples
psivec <- annotatePsi(psivec)
psivec_norm <- normalizePerTranscriptType(psivec, samples)
message("writing normalized psis")
sapply(unique(meta$condition), function(ct) {
  s <- meta[meta$condition==ct,"sample"]
  write.table(psivec_norm[,s], file = paste0(ct, ".psi"),
              sep = "\t", quote = F)
})

sumPerTranscriptType(psivec, normalized = F)
sumPerTranscriptType(psivec_norm, normalized = T)
