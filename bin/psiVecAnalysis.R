library(dplyr)
library(stringr)
library(EnhancedVolcano)
library(biomaRt)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(MetBrewer)


readdPsi <- function(path) {
  dPsi <- read.csv(path, sep = "\t", check.names = F)
  dPsi$conditions <- sub("_dPSI", "", colnames(dPsi)[2])
  colnames(dPsi)[2:3] <- c("dPSI", "p.val")
  dPsi
}

annotatePsi <- function(psivec) {
  psivec[,c("gene_id", "transcript_id")] <- str_split_fixed(rownames(psivec), ";", 2)
  psivec$transcript_type <- ifelse(grepl("c", psivec$transcript_id),
                                   "circRNA", "linearRNA")
  psivec[is.na(psivec)] <- 0
  psivec
}

normalizePerTranscriptType <- function(psivec, samples) {
  message("normalize psi values for mean of psi of linear and circular transcripts")
  means <- psivec %>% group_by(transcript_type) %>% summarise_at(vars(samples), mean)
  
  psivec[psivec$transcript_type=="circRNA",samples] <- psivec[psivec$transcript_type=="circRNA",samples] /
    as.vector(means[means$transcript_type=="circRNA", samples])
  psivec[psivec$transcript_type=="linearRNA",samples] <- psivec[psivec$transcript_type=="linearRNA",samples] / 
    as.vector(means[means$transcript_type=="linearRNA", samples])
  psivec %>% group_by(gene_id) %>% mutate(across(samples, function(s) s/sum(s)))
}
#### CLI ARGS ####
args = commandArgs(trailingOnly = TRUE)

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
psivec <- normalizePerTranscriptType(psivec, samples)
message("writing normalized psis")
sapply(unique(meta$condition), function(ct) {
  s <- meta[meta$condition==ct,"sample"]
  write.table(psivec[,s], file = paste0(ct, ".psi"),
              sep = "\t", quote = F)
})
# group data by host gene and summarize per transcript type if the gene contains circular transcripts
rnaTypeRatios <- psivec %>%
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
  scale_y_continuous(breaks = seq(0, 100, 25)) +
  facet_wrap(~variable) +
  xlab("Transcript Type") +
  ylab("Percent Spliced In (PSI) value, normalized by sample-wise mean of transcript type") +
  theme(text = element_text(size = 12), strip.text.x = element_text(size = 12))

ggsave(filename = "violins.png", plot = boxP, width = 8, height = 12, dpi = 1624)
