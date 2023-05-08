library(dplyr)
library(argparser)


readdPsi <- function(path) {
  dPsi <- read.csv(path, sep = "\t", check.names = F)
  dPsi$conditions <- sub("_dPSI", "", colnames(dPsi)[2])
  colnames(dPsi)[2:3] <- c("dPSI", "p.val")
  dPsi
}

#### CLI ARGS ####
args = commandArgs(trailingOnly = TRUE)

parser <- arg_parser("Argument parser for differenial splicing analysis", name = "DE_parser")
parser <- add_argument(parser, "-m", help = "Meta file")
parser <- add_argument(parser, "-i", help = "Suppa2 diffSplice output directory", default = "./")
# parameter settings
parser <- add_argument(parser, "-p", help = "P-value cutoff", default = 0.05)
parser <- add_argument(parser, "-d", help = "Minimum difference in psi values between conditions",
                       default = 0.25)

argv <- parse_args(parser, argv = args)

pValue = argv$p
mindPsi = argv$d
# read samplesheet
meta <- read.csv(argv$m, sep = "\t")
conditions <- unique(meta$condition)
# read diffSplice dpsi results
dpsiFiles <- list.files(path = argv$i, pattern = ".dpsi", full.names = T)
dpsiPerConditions <- do.call(rbind, lapply(dpsiFiles, readdPsi))

# split event into gene and transcript ids
dpsiPerConditions[,c("gene_id", "transcript_id")] <- str_split_fixed(dpsiPerConditions$Event_id, ";", 2)
dpsiPerConditions$Event_id <- NULL

# split conditions
cp <- paste0("(", paste(conditions, collapse = "|"), ")")
expr <- paste0("(?<=", cp, ")", "(-)", "(?=", cp, ")")
dpsiPerConditions[, c("Cond1", "Cond2")] <- str_split_fixed(dpsiPerConditions$conditions, regex(expr), 2)
dpsiPerConditions$transcript_type <- ifelse(grepl("c", dpsiPerConditions$transcript_id),
                                                  "circRNA", "linearRNA")

dpsiPerConditions$Spsi <- ifelse(abs(dpsiPerConditions$dPSI) >= mindPsi, "p", "")
dpsiPerConditions$Spval <- ifelse(dpsiPerConditions$p.val <= pValue, "v", "")
dpsiPerConditions$level <- paste0(dpsiPerConditions$Spsi,"_", dpsiPerConditions$Spval)

dpsiPerConditions$p.val <- -log10(dpsiPerConditions$p.val)
signifColors <- setNames(c("darkgrey", "blue2", "darkgreen", "red3"), c("_", "p_", "_v", "p_v"))
signifLabels <- c("not significant",
                  paste0("p-value >= ", pValue),
                  paste0("dPSI >= ", mindPsi),
                  "both")
# volcano plot  
p <- ggplot(dpsiPerConditions, aes(x = dPSI, y = p.val)) +
  geom_point(size = 2, aes(color = level)) +
  geom_hline(yintercept = -log10(pValue), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = -mindPsi, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = mindPsi, linetype = "dashed", color = "grey") +
  scale_color_manual(values = signifColors,
                     labels = signifLabels) +
  ylab("-log10(p-value)") + xlab("PSI (cond2) - PSI (cond1)") +
  theme(legend.title = element_blank(), 
        text = element_text(size = 14), strip.text.x = element_text(size = 12), 
        legend.position = "bottom", legend.direction = "horizontal") +
  facet_wrap(transcript_type~conditions)
ggsave("volcanos.png", plot = p, dpi = 300, width = 8, height = 6)

# filter for pvalue of <= 0.05 and ΔPSI >= 50%
dpsiPerConditions %>% filter(transcript_type == "circRNA") %>%
  filter(p.val >= -log10(pValue), abs(dPSI) >= mindPsi) %>% 
  arrange(desc(p.val), abs(dPSI)) %>%
  write.table(file = "signifCircRNAs.tsv", sep = "\t", row.names = F, quote = F)
dpsiPerConditions %>% filter(transcript_type == "linearRNA") %>%
  filter(p.val >= -log10(pValue), abs(dPSI) >= mindPsi) %>% 
  arrange(desc(p.val), abs(dPSI)) %>%
  write.table(file = "signifLinearRNAs.tsv", sep = "\t", row.names = F, quote = F)
