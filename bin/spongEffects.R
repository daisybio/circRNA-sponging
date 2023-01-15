#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(SPONGE, doParallel, foreach, dplyr, argparser, visNetwork, MetBrewer, pheatmap)

########################
## GENERAL FUNCTIONS ###
########################
# plot ceRNA network with gene annotations
plot_network <- function(ceRNA_network, signif_hits=NULL, gtf=NULL, annotation=NULL) {
  # plot network
  ceRNA_plot <- sponge_plot_network(ceRNA_network, genes_miRNA_candidates, ) %>%
    visNetwork::visEdges(arrows = list(to = list(enabled = T, scaleFactor = 1)))
  ceRNA_plot$x$edges$label <- paste("mscor:", round(ceRNA_network$mscor, 2))
  
  # extract nodes and edges for customization
  nodes <- ceRNA_plot$x$nodes
  edges <- ceRNA_plot$x$edges
  
  # mark differentially expressed RNAs
  if (!is.null(signif_hits)){
    nodes[nodes$id %in% hgncs$hgnc_symbol | nodes$id %in% signif_hits$X,"group"] <- "DE"
  }
  if (!is.null(gtf)) {
    gtf <- rtracklayer::readGFF(gtf)
    annotation <- unique(gtf[!is.na(gtf$transcript_id),c("gene_id", "gene_name", "gene_biotype")])
    colnames(annotation) <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
    rownames(annotation) <- annotation$ensembl_gene_id
  }
  if (!is.null(annotation)){
    # break types of RNA down to lncRNA, coding, and rest
    annotation$gene_biotype[annotation$gene_biotype != "protein_coding" & annotation$gene_biotype != "lncRNA"] <- "other_RNA"
    
    nodes <- merge(nodes, annotation, by = 1, all.x = T)
    
    # change to hgnc
    nodes[!is.na(nodes$hgnc_symbol),1:2] <- nodes[!is.na(nodes$hgnc_symbol),"hgnc_symbol"]
    
    # add circRNA as biotype
    nodes[is.na(nodes$gene_biotype) & !grepl("ENSG", nodes$id),"gene_biotype"] <- "circRNA"
    # label unknown RNAs as other
    nodes[is.na(nodes$gene_biotype), "gene_biotype"] <- "other_RNA"
    # remove preset color and shape
    nodes <- nodes[,-c(3,4)]
    # change to group
    colnames(nodes)[8] <- "group"
    
    # convert geneA and geneB
    edges <- merge(edges, annotation, by.x = "from", by.y = 1, all.x = T)
    edges[!is.na(edges$hgnc_symbol),"from"] <- edges$hgnc_symbol[!is.na(edges$hgnc_symbol)]
    edges <- merge(edges, annotation, by.x = "to", by.y = 1, all.x = T)
    edges[!is.na(edges$hgnc_symbol.y),"to"] <- edges$hgnc_symbol.y[!is.na(edges$hgnc_symbol.y)]
  }
  # plot final graph
  graph <- visNetwork(nodes = nodes, edges = edges) %>%
    visIgraphLayout(type = "full", physics = F) %>%
    visGroups(groupname = "circRNA", shape = "rectangle", color = "#33FF99") %>%
    visGroups(groupname = "protein_coding", shape = "triangle", color = "#0066CC") %>%
    visGroups(groupname = "lncRNA", shape = "square", color = "#F8766D") %>%
    visGroups(groupname = "other_RNA", color = "#DC71FA", shape = "diamond") %>%
    visGroups(groupname = "DE", color = "#CC3333") %>%
    visLegend()
  return(graph)
}

# plot model performance (spongEffects)
plot_performance <- function (trained_model, central_genes_model = NULL, random_model, 
                              training_dataset_name = "TCGA", testing_dataset_name = "TCGA", 
                              subtypes)
{
  trained.model <- trained_model
  CentralGenes.model <- central_genes_model
  Random.model <- random_model
  training_string <- paste0(training_dataset_name, " (Training)")
  testing_string <- paste0(testing_dataset_name, " (Testing)")
  set_names <- c(training_string, testing_string)
  SpongingActiivty.model <- trained_model
  if (!is.null(central_genes_model)) {
    Accuracy.df <- data.frame(Run = rep(set_names, 3), Model = c(rep("Modules", 
                                                                     2), rep("Central Genes", 2), rep("Random", 2)), 
                              Accuracy = c(SpongingActiivty.model$ConfusionMatrix_training[["overall"]][["Accuracy"]], 
                                           SpongingActiivty.model$ConfusionMatrix_testing[["overall"]][["Accuracy"]], 
                                           CentralGenes.model$ConfusionMatrix_training[["overall"]][["Accuracy"]], 
                                           CentralGenes.model$ConfusionMatrix_testing[["overall"]][["Accuracy"]], 
                                           Random.model$ConfusionMatrix_training[["overall"]][["Accuracy"]], 
                                           Random.model$ConfusionMatrix_testing[["overall"]][["Accuracy"]]))
    Accuracy.df$Model <- factor(Accuracy.df$Model, levels = c("Modules", 
                                                              "Random", "Central Genes"))
    Accuracy.df$Run <- factor(Accuracy.df$Run, levels = set_names)
  }
  else {
    Accuracy.df <- data.frame(Run = rep(set_names, 2), Model = c(rep("Modules", 
                                                                     2), rep("Random", 2)), Accuracy = c(SpongingActiivty.model$ConfusionMatrix_training[["overall"]][["Accuracy"]], 
                                                                                                         SpongingActiivty.model$ConfusionMatrix_testing[["overall"]][["Accuracy"]], 
                                                                                                         Random.model$ConfusionMatrix_training[["overall"]][["Accuracy"]], 
                                                                                                         Random.model$ConfusionMatrix_testing[["overall"]][["Accuracy"]]))
    Accuracy.df$Model <- factor(Accuracy.df$Model, levels = c("Modules", 
                                                              "Random"))
    Accuracy.df$Run <- factor(Accuracy.df$Run, levels = set_names)
  }
  Accuracy.plot <- Accuracy.df %>% ggplot(aes(x = Accuracy, 
                                              y = Model)) + geom_line(aes(group = Model)) + geom_point(aes(shape = Run)) + 
    theme_light() + xlab("Subset Accuracy") + ylab("") + 
    theme_bw()
  Metrics.SpongeModules.training <- SpongingActiivty.model$ConfusionMatrix_training$byClass[,"Balanced Accuracy"] %>% 
    as.data.frame() %>% mutate(Model = "Modules") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.SpongeModules.training) = c("Class", "Value", 
                                               "Model")
  Metrics.Random.training <- Random.model$ConfusionMatrix_training$byClass[,"Balanced Accuracy"] %>% 
    as.data.frame() %>% mutate(Model = "Random") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.Random.training) = c("Class", "Value", 
                                        "Model")
  if (!is.null(central_genes_model)) {
    Metrics.CentralGenes.training <- CentralGenes.model$ConfusionMatrix_training[["byClass"]][c(1:length(unique(subtypes)))] %>% 
      as.data.frame() %>% mutate(Model = "Central Genes") %>% 
      tibble::rownames_to_column("Class")
    colnames(Metrics.CentralGenes.training) = c("Class", 
                                                "Value", "Model")
    Metrics.training <- rbind(Metrics.SpongeModules.training, 
                              rbind(Metrics.Random.training, Metrics.CentralGenes.training)) %>% 
      mutate(Run = training_string)
  }
  else {
    Metrics.training <- rbind(Metrics.SpongeModules.training, 
                              Metrics.Random.training) %>% mutate(Run = training_string)
  }
  Metrics.SpongeModules.testing <- SpongingActiivty.model$ConfusionMatrix_testing$byClass[,"Balanced Accuracy"] %>% 
    as.data.frame() %>% mutate(Model = "Modules") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.SpongeModules.testing) = c("Class", "Value", 
                                              "Model")
  Metrics.Random.testing <- Random.model$ConfusionMatrix_testing$byClass[,"Balanced Accuracy"] %>% 
    as.data.frame() %>% mutate(Model = "Random") %>% tibble::rownames_to_column("Class")
  colnames(Metrics.Random.testing) = c("Class", "Value", "Model")
  if (!is.null(central_genes_model)) {
    Metrics.CentralGenes.testing <- CentralGenes.model$ConfusionMatrix_testing[["byClass"]][c(1:length(unique(subtypes)))] %>% 
      as.data.frame() %>% mutate(Model = "Central Genes") %>% 
      tibble::rownames_to_column("Class")
    colnames(Metrics.CentralGenes.testing) = c("Class", 
                                               "Value", "Model")
    Metrics.testing <- rbind(Metrics.SpongeModules.testing, 
                             rbind(Metrics.Random.testing, Metrics.CentralGenes.testing)) %>% 
      mutate(Run = testing_string)
  }
  else {
    Metrics.testing <- rbind(Metrics.SpongeModules.testing, 
                             Metrics.Random.testing) %>% mutate(Run = testing_string)
  }
  Metrics <- rbind(Metrics.training, Metrics.testing)
  if (!is.null(central_genes_model)) {
    Metrics$Model <- factor(Metrics$Model, levels = c("Modules", 
                                                      "Random", "Central Genes"))
  }
  else {
    Metrics$Model <- factor(Metrics$Model, levels = c("Modules", 
                                                      "Random"))
  }
  Metrics$Run <- factor(Metrics$Run, levels = set_names)
  Metrics$Class <- gsub("Class: ", "", Metrics$Class)
  Metrics.plot <- Metrics %>% ggplot(aes(x = Class, y = Value, 
                                         fill = Model)) + geom_bar(position = "dodge", stat = "identity", 
                                                                   width = 0.5) + facet_grid(Metrics$Run) + xlab("Accuracy") + 
    ylab("") + theme_bw()
  metric_plots <- ggarrange(Accuracy.plot, Metrics.plot, ncol = 1, 
                            nrow = 2)
  return(metric_plots)
}

# plot lollipop (spongEffects)
plot_modules <- function (trained_model, k_modules = 25, k_modules_red = 10,
          text_size = 16)
{
  final.model <- trained_model$Model$finalModel
  Variable.importance <- importance(final.model) %>% as.data.frame() %>% 
    tibble::rownames_to_column("Module") %>%
    arrange(desc(MeanDecreaseGini))
  grey_modules = k_modules - k_modules_red
  p <- Variable.importance[1:k_modules, ] %>% mutate(Analysed = c(rep("1", 
                                                                      k_modules_red), rep("0", grey_modules))) %>% arrange(desc(MeanDecreaseGini)) %>% 
    ggplot(aes(x = reorder(Module, MeanDecreaseGini), y = MeanDecreaseGini)) + 
    geom_point() + geom_segment(aes(x = Module, xend = Module, 
                                    y = 0, yend = MeanDecreaseGini, color = Analysed)) + 
    scale_colour_manual(values = c("red", "black"), breaks = c("1", 
                                                               "0")) + coord_flip() + xlab("Module") + ylab("Mean decrease in Gini index") + 
    theme_light() + theme(panel.grid.major.x = element_blank(), 
                          panel.grid.minor.x = element_blank(), axis.ticks.y = element_blank(), 
                          legend.title = element_blank(), legend.position = "none", 
                          legend.background = element_blank(), legend.direction = "horizontal", 
                          panel.border = element_rect(colour = "black", fill = NA, 
                                                      size = 1), text = element_text(size = 16))
  return(p)
}

# Training heat map (spongEffects)
plot_hmap <- function (trained_model, spongEffects, meta_data, label, sampleIDs, 
          Modules_to_Plot = 5, show.rownames = F, show.colnames = F) 
{
  if (label %in% colnames(meta_data) & sampleIDs %in% colnames(meta_data)) {
    final.model <- trained_model$Model$finalModel
    Variable.importance <- importance(final.model) %>% as.data.frame() %>% 
      tibble::rownames_to_column("Module") %>% arrange(desc(MeanDecreaseGini))
    Variable.importance$Module <- gsub("`", "", Variable.importance$Module)
    Annotation.meta <- meta_data[match(colnames(spongEffects), 
                                       meta_data[, sampleIDs]), ]
    unique_subtypes <- unique(Annotation.meta$label)
    number_groups <- length(unique_subtypes)
    col.heatmap <- met.brewer("Renoir", n = number_groups, 
                              type = "continuous")
    col.heatmap <- as.vector(col.heatmap)
    col.heatmap <- setNames(col.heatmap, unique_subtypes)
    col.heatmap <- cell.colors
    Column.Annotation <- HeatmapAnnotation(Group = Annotation.meta[, 
                                                                   label], col = list(Group = col.heatmap))
    spongeEffects.toPlot <- spongEffects[match(Variable.importance$Module, 
                                               rownames(spongEffects)), ]
    print(spongeEffects.toPlot)
    spongeEffects.toPlot <- spongeEffects.toPlot[rowSums(is.na(spongeEffects.toPlot)) == 
                                                   0, ]
    if (Modules_to_Plot > length(rownames(spongeEffects.toPlot))) 
      Modules_to_Plot = length(rownames(spongeEffects.toPlot))
    spongeEffects.toPlot <- spongeEffects.toPlot[1:Modules_to_Plot, 
    ]
    Heatmap.p <- spongeEffects.toPlot %>% t() %>% scale() %>% 
      t() %>% Heatmap(show_row_names = show.rownames, 
                      show_column_names = show.colnames, top_annotation = Column.Annotation, 
                      show_heatmap_legend = TRUE)
    return(Heatmap.p)
  }
  else {
    print("label and/or sampleIDs must be columns in metadata")
  }
}

# Target gene heat map (spongEffects)
plot_target_gene_expressions <- function(target, target_genes, gene_expression, meta, 
                                         log_transform = T, pseudocount = 1, gtf_raw = NULL, 
                                         annotation = NULL, split = "condition", unit = "counts",
                                         show_rows = T, annotation_colors = NULL) {
  # get target expression
  target_expression <- gene_expression[,target, drop = F]
  # get target genes expressions
  target_genes_expression <- gene_expression[,target_genes]
  # combine into one
  data <- t(cbind(target_expression, target_genes_expression))
  # save all conditions of samples
  conditions <- unique(meta[,split])
  
  # convert to hgnc if gtf is given
  if (!is.null(gtf_raw)) {
    print("reading GTF file and converting to HGNC symbols...")
    gtf <- rtracklayer::readGFF(gtf)
    gene.ens.all <- unique(gtf[!is.na(gtf$transcript_id),c("gene_id", "gene_name")])
    colnames(gene.ens.all) <- c("ensembl_gene_id", "hgnc_symbol")
    rownames(gene.ens.all) <- gene.ens.all$ensembl_gene_id
    # convert
    annotation = gene.ens.all
  }
  if (!is.null(annotation)){
    ensgs <- rownames(data)
    ensgs <- merge(ensgs, annotation, by = 1, all.x = T)
    ensgs[!is.na(ensgs$hgnc_symbol),"x"] <- ensgs[!is.na(ensgs$hgnc_symbol),"hgnc_symbol"]
    rownames(data) <- ensgs$x
  }
  # heat color scheme
  colors <- c(colorRampPalette(c("blue", "orange"))(100), colorRampPalette(c("orange", "red"))(100))
  # annotation colors are given
  if (!is.null(annotation_colors)) {
    annotation.colors <- annotation_colors
  } else {
    annotation.colors <- hcl.colors(length(conditions), palette = hcl.pals(type = "diverging")[12])
  }
  # name colors
  names(annotation.colors) <- conditions
  a_c <- list(x=annotation.colors)
  names(a_c) <- split

  if (log_transform) {
    data = log2(data + pseudocount)
  }
  colnames(data) <- gsub("\\.", "-", colnames(data))
  label <- paste0("Module ", unit, " expression for: ", target)
  # create annotation for heat map
  df <- data.frame(meta[,c("sample", split)], row.names = 1)
  pheatmap::pheatmap(data, treeheight_row = 0, treeheight_col = 0,
                     show_colnames = F, show_rownames = show_rows, 
                     cluster_rows = T, cluster_cols = T,
                     color = colors, annotation_col = df,
                     annotation_colors = a_c, 
                     main = label, fontsize_row = 10)
}

args.effects <- commandArgs(trailingOnly = TRUE)

effects_parser <- arg_parser("Argument parser for spongEffects module", name = "spongEffects_parser")
effects_parser <- add_argument(effects_parser, "--spongeData", help = "SPONGE Rdata containing all results, e.g. gene expression, miRNA expression, sponge centralities etc.")
effects_parser <- add_argument(effects_parser, "--meta", help = "Meta data for samples in tsv")
effects_parser <- add_argument(effects_parser, "--gtf", help = "Gene annotation file (GTF)")
########################
##  PARAMETER TUNING  ##
########################
effects_parser <- add_argument(effects_parser, "--train", help = "Percentage of samples to use for training, rest will be used for testing", default = 0.8)
effects_parser <- add_argument(effects_parser, "--cpus", help = "Number of cores to use for backend", default = 4)
effects_parser <- add_argument(effects_parser, "--mscor", help = "Mscore threshold", default = 0.1)
effects_parser <- add_argument(effects_parser, "--genome", help = "Genome version, e.g. GRCh38", default = 0.1)
effects_parser <- add_argument(effects_parser, "--fdr", help = "False discovery rate for padj ceRNA interactions", default = 0.05)
effects_parser <- add_argument(effects_parser, "--modules_cutoff", help = "Modules cutoff", default = 750)
effects_parser <- add_argument(effects_parser, "--bin_size", help = "Total bin size for enrichment", default = 100)
effects_parser <- add_argument(effects_parser, "--min_size", help = "Minimum size for enrichment", default = 100)
effects_parser <- add_argument(effects_parser, "--max_size", help = "Maximum size for enrichment", default = 2000)
effects_parser <- add_argument(effects_parser, "--min_expr", help = "Minimum expression for enrichment", default = 10)
effects_parser <- add_argument(effects_parser, "--method", help = "Method", default = "OE")
effects_parser <- add_argument(effects_parser, "--enrichment_cores", help = "Number of cores to use for enrichment", default = 25)
effects_parser <- add_argument(effects_parser, "--cv_folds", help = "Number of cross validation folds", default = 10)
effects_parser <- add_argument(effects_parser, "--out", help = "Output directory", default = "./")
# parse arguments
argv_eff <- parse_args(effects_parser, argv = args.effects)

#---------------------------OUTPUT DIR--------------------------------
ROOT = argv_eff$out
# create plot directory
PLOT_DIR = file.path(ROOT, "plots")
dir.create(PLOT_DIR, recursive = T, showWarnings = F)
#---------------------------SPONGE DATA-------------------------------
print("loading SPONGE data...")
load(argv_eff$spongeData)

#---------------------------GTF DATA----------------------------------
gtf <- rtracklayer::readGFF(argv_eff$gtf)
gene.ens.all <- unique(gtf[!is.na(gtf$transcript_id),c("gene_id", "gene_name", "gene_biotype")])
colnames(gene.ens.all) <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")
rownames(gene.ens.all) <- gene.ens.all$ensembl_gene_id

#---------------------------PARAMETERS--------------------------------
mscor.threshold = argv_eff$mscor
padj.threshold = argv_eff$fdr
modules_cutoff = argv_eff$modules_cutoff
bin.size = argv_eff$bin_size
min.size = argv_eff$min_size
max.size = argv_eff$max_size
min.expr = argv_eff$min_expr
method = argv_eff$method
enrichment_cores = argv_eff$enrichment_cores
split_training = argv_eff$train
folds = argv_eff$cv_folds

#---------------------------PARALLEL BACKGROUND-----------------------
num.of.cores <- argv_eff$cpus
cat("registering back end with", num.of.cores, "cores\n")
cl <- makeCluster(num.of.cores) 
registerDoParallel(cl)

#---------------------------META FILE/ CONDITIONS---------------------
meta <- read.csv(file = argv$meta, sep = "\t")
# add cell type column
meta$cell_type <- sapply(strsplit(meta$sample, "_"), "[", 1)

# ----------------EXPRESSION SPLITTING--------------------------------
n.train <- round(nrow(meta)*split_training)
n.test <- round(nrow(meta)-n.train)
cat("using", n.train, "samples for training\n", "and", n.test, "for testing\n")
# randomize samples
meta <- meta[sample(1:nrow(meta)), ]
# take n samples from each group
cond.split <- split(meta, meta$cell_type)
cond.ratio <- round(sapply(cond.split, function(x) nrow(x)*split_training))
train.meta <- data.table::rbindlist(mapply(function(x,y) x[1:y,], cond.split, cond.ratio, SIMPLIFY = F))
colnames(train.meta)<-c("sampleIDs",c(colnames(train.meta)[3:length(train.meta)-1],"label"))
test.meta <- data.table::rbindlist(mapply(function(x,y) x[(y+1):nrow(x),], cond.split, cond.ratio, SIMPLIFY = F))
colnames(test.meta)<-c("sampleIDs",c(colnames(test.meta)[3:length(test.meta)],"label"))

rownames(gene_expr) <- gsub("\\.", "-", rownames(gene_expr))
# train gene expression; split_training of samples
train_gene_expr <- gene_expr[rownames(gene_expr) %in% train.meta$sampleIDs,]
train_gene_expr<-t(train_gene_expr)
# test gene expression; rest of samples
test_gene_expr <- gene_expr[rownames(gene_expr) %in% test.meta$sampleIDs,]
test_gene_expr<-t(test_gene_expr)

# train miRNA expression
train_mirna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% train.meta$sampleIDs,]
# test miRNA expression
test_mirna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% test.meta$sampleIDs,]

#-------------------------CE RNA SPLITTING------------------------
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < padj.threshold),]
ceRNA_interactions_fdr <- circ.mRNA.subnetwork(ceRNA_interactions_fdr, "c")
n.train.ce <- round(nrow(ceRNA_interactions_fdr)*split_training)
n.test.ce <- nrow(ceRNA_interactions_fdr)-n.train.ce
train_ceRNA_interactions <- head(ceRNA_interactions_fdr, n.train.ce)
test_ceRNA_interactions <- tail(ceRNA_interactions_fdr, n.test.ce)
#-------------------------CENTRALITIES SPLITTING------------------
network_centralities  <- sponge_node_centralities(ceRNA_interactions_fdr)
n.train.c <- round(nrow(network_centralities)*split_training)
n.test.c <- nrow(network_centralities)-n.train.c
train_network_centralities <- head(network_centralities, n.train.c)
test_network_centralities <- tail(network_centralities, n.test.c)
#-------------------------FILTERING-------------------------------
print("filtering centralities...")
filtered_network_centralities <- filter_ceRNA_network(sponge_effects = train_ceRNA_interactions, 
                                                      Node_Centrality = train_network_centralities,
                                                      add_weighted_centrality=T, 
                                                      mscor.threshold = mscor.threshold, 
                                                      padj.threshold = padj.threshold)
#-------------------------RNAS_OF_INTEREST------------------------
RNAs <- c("lncRNA","circRNA")
RNAs.ofInterest <- ensembl.df %>% dplyr::filter(gene_biotype %in% RNAs) %>%
  dplyr::select(ensembl_gene_id)
# add circRNAs of the data set
new_circRNAs <- data.frame(ensembl_gene_id=colnames(gene_expr)[grep("c", colnames(gene_expr))])
RNAs.ofInterest <- rbind(RNAs.ofInterest, new_circRNAs)
# get central modules
central_gene_modules <- get_central_modules(central_nodes = RNAs.ofInterest$ensembl_gene_id,
                                            node_centrality = filtered_network_centralities$Node_Centrality,
                                            ceRNA_class = RNAs, 
                                            centrality_measure = "Weighted_Degree", 
                                            cutoff = modules_cutoff)
#-------------------------SPONGE MODULES---------------------------
Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, 
                                 central.modules = central_gene_modules,
                                 remove.central = F,
                                 set.parallel = F)
# Module size distribution
Size.modules <- sapply(Sponge.modules, length)
#-------------------------SPLIT MODULES----------------------------
# train and test modules
train.modules <- enrichment_modules(Expr.matrix = train_gene_expr, modules = Sponge.modules, 
                                    bin.size = bin.size, min.size = min.size, max.size = max.size, 
                                    min.expr = min.expr, method = method, cores=enrichment_cores)
test.modules <-  enrichment_modules(Expr.matrix = test_gene_expr, modules = Sponge.modules, 
                                    bin.size = bin.size, min.size = min.size, max.size = max.size, 
                                    min.expr = min.expr, method = method, cores=enrichment_cores)
#-------------------------MODEL PERFORMANCE--------------------------
common_modules = intersect(rownames(train.modules), rownames(test.modules))
train.modules = train.modules[common_modules, ]
test.modules = test.modules[common_modules, ]
trained.model = calibrate_model(Input = train.modules, modules_metadata = train.meta, label = "label", 
                                sampleIDs = "sampleIDs", 
                                Metric = "Exact_match", n_folds = folds, repetitions = 3)
trained.model[["ConfusionMatrix_training"]]

Input.test <- t(test.modules) %>% scale(center = T, scale = T)
Prediction.model <- predict(trained.model$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test.meta$label))
trained.model$ConfusionMatrix_testing <- ConfusionMatrix_testing
#-------------------------RANDOM MODULES----------------------------
Random.modules <- Random_spongEffects(sponge_modules = Sponge.modules,
                                      gene_expr = train_gene_expr, min.size = min.size,bin.size = bin.size, max.size = max.size,
                                      min.expression=min.expr, replace = F,method = method,cores = enrichment_cores)
# We can now use the randomly defined modules to calculate their enrichment in the test set
Random.modules.test <- enrichment_modules(Expr.matrix = test_gene_expr, modules = Random.modules$Random_Modules, 
                                          bin.size = bin.size, min.size = min.size, max.size = max.size, min.expr = min.expr, 
                                          method = method, cores=enrichment_cores)
# We find random modules that were identified both in the train and test and use those as input features for the model
common_modules_random = intersect(rownames(Random.modules$Enrichment_Random_Modules), rownames(Random.modules.test))
Random.modules.train = Random.modules$Enrichment_Random_Modules[common_modules_random, ]
Random.modules.test = Random.modules.test[common_modules_random, ]
Random.model = calibrate_model(Input = Random.modules.train, modules_metadata = train.meta, 
                               label = "label", sampleIDs = "sampleIDs",
                               Metric = "Exact_match", n_folds = folds, repetitions = 1)
Random.model[["ConfusionMatrix_training"]]

# validate on test set
Input.test <- t(Random.modules.test) %>% scale(center = T, scale = T)
Input.test<-Input.test[ , apply(Input.test, 2, function(x) !any(is.na(x)))]
Prediction.model <- predict(Random.model$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing_random <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test.meta$label))
Random.model$ConfusionMatrix_testing_random<-ConfusionMatrix_testing_random
ConfusionMatrix_testing_random

#-------------------------PLOT PERFOMANCES-------------------------
train_name <- paste0(split_training * 100, "%")
test_name <- paste0((1-split_training) * 100, "%")
metrics_plot <- plot_performance(trained_model =  trained.model,
                                random_model = Random.model, 
                                central_genes_model = NULL,
                                training_dataset_name = train_name,
                                testing_dataset_name = test_name,
                                subtypes=as.factor(train.meta$label))
# save metrics plot
ggsave(file.path(PLOT_DIR, "metrics.png"), metrics_plot,
       width = 7.25, height = 5.25, dpi = 1200)

#-------------------------PLOT CLASS DISTRIBUTION------------------
density_plot_train <- plot_density_scores(trained_model = trained.model, spongEffects = train.modules, 
                                          meta_data = train.meta, label = "label", sampleIDs = "sampleIDs")
# save class plot
ggsave(file.path(PLOT_DIR, "classification.png"), density_plot_train,
       width = 7.25, height = 5.25, dpi = 1200)

#-------------------------PLOT TOP RESULTS-------------------------
lollipop_plot <- plot_modules(trained_model = trained.model, k_modules_red = 2,
                             k_modules = 7, text_size = 20)
# save lollipop plot
ggsave(file.path(PLOT_DIR, "lollipop.png"), lollipop_plot,
       width = 7.25, height = 5.25, dpi = 1200)

#-------------------------PLOT NETWORK FOR CENTRALITIES------------
n = 25
central_players <- central_gene_modules %>% 
  arrange(desc(Weighted_Degree)) %>% 
  dplyr::select(gene) %>% slice_head(n = n) %>% unlist
# filter network for central players, filter for circRNA-mRNA connections
network <- circ.mRNA.subnetwork(filtered_network_centralities$Sponge.filtered, "c")
subnetwork <- network %>%
  dplyr::filter(geneA %in% central_players | geneB %in% central_players)
network_plot <- plot_network(subnetwork,
                             annotation = gene.ens.all,
                             node_label_size = 40, edge_label_size = 0)
# save network
visNetwork::visSave(network_plot, file = file.path(PLOT_DIR, paste0("top_", n, "_ceRNA_centralities.html")))

#-------------------------EXTRACT TOP circRNAs---------------------
Variable.importance <- importance(trained.model$Model$finalModel) %>% as.data.frame() %>% 
  tibble::rownames_to_column("Module") %>% arrange(desc(MeanDecreaseGini))
k = 6
candidates <- Variable.importance[1:k, "Module"]
candidates <- gsub("`", "", candidates)
# candidate modules
candidate.modules <- Sponge.modules[candidates]
# continue with top 2 (highest gini indices)
sponged.genes <- candidate.modules[candidates]

#-------------------------PLOT TOP circRNA modules-----------------
top_network_plot <- plot_network(ceRNA_network = network %>% 
                                dplyr::filter(geneA %in% candidates | geneB %in% candidates),
                                annotation = gene.ens.all, node_label_size = 40, edge_label_size = 0)
# plot top k centralities network
visNetwork::visSave(top_network_plot, file = file.path(PLOT_DIR, paste0("top_", k, "_ceRNA_centralities.html")))

#-------------------------MODULE EXPRESSION------------------------
# generate module expression for each candidate
candidate_plots <- list()
for (candidate in candidates) {
  module <- candidate.modules[candidate][[1]]
  module <- module[!grepl("c", module)]
  module_gene_plot <- plot_target_gene_expressions(candidate, module,
                                                  gene_expression = gene_expr, meta = meta,
                                                  annotation = gene.ens.all,
                                                  unit = "counts", log_transform = T, show_rows = F)
  module_miRNA_plot <- plot_miRNAs_per_gene(candidate, genes_miRNA_interactions_manual,
                                            mi_rna_expr = mi_rna_expr, meta = meta,
                                            log_transform = T)
  candidate_plots[[candidate]] <- list(module_gene_plot, module_miRNA_plot)
}
#-------------------------spongEffects HMAP------------------------
train.hmap <- plot_hmap(trained_model = trained.model, spongEffects = train.modules,
                        meta_data = train.meta, label = "label", sampleIDs = "sampleIDs", Modules_to_Plot = 2,
                        show.rownames = T, show.colnames = F)
test.hmap <- plot_hmap(trained_model = trained.model, spongEffects = test.modules,
                        meta_data = test.meta, label = "label", sampleIDs = "sampleIDs", Modules_to_Plot = 2,
                        show.rownames = T, show.colnames = F)
# save train and test heat maps
ggsave(file.path(PLOT_DIR, "train_hm.png"), train.hmap,
       width = 7.25, height = 5.25, dpi = 1200)
ggsave(file.path(PLOT_DIR, "test_hm.png"), test.hmap,
       width = 7.25, height = 5.25, dpi = 1200)
