#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(SPONGE, doParallel, foreach, dplyr, argparser)

args.effects = commandArgs(trailingOnly = TRUE)

parser.e <- arg_parser("Argument parser for spongEffects module", name = "spongEffects_parser")
parser.e <- add_argument(parser.e, "--spongeData", help = "SPONGE Rdata containing all results, e.g. gene expression, miRNA expression, sponge centralities etc.")
parser.e <- add_argument(parser.e, "--meta", help = "Meta data for samples in tsv")
########################
##  PARAMETER TUNING  ##
########################
parser.e <- add_argument(parser.e, "--train", help = "Percentage of samples to use for training, rest will be used for testing", default = 0.8)
parser.e <- add_argument(parser.e, "--cpus", help = "Number of cores to use for backend", default = 4)
parser.e <- add_argument(parser.e, "--mscor", help = "Mscore threshold", default = 0.1)
parser.e <- add_argument(parser.e, "--genome", help = "Genome version, e.g. GRCh38", default = 0.1)
parser.e <- add_argument(parser.e, "--fdr", help = "False discovery rate for padj ceRNA interactions", default = 0.05)
parser.e <- add_argument(parser.e, "--modules_cutoff", help = "Modules cutoff", default = 750)
parser.e <- add_argument(parser.e, "--bin_size", help = "Total bin size for enrichment", default = 100)
parser.e <- add_argument(parser.e, "--min_size", help = "Minimum size for enrichment", default = 100)
parser.e <- add_argument(parser.e, "--max_size", help = "Maximum size for enrichment", default = 2000)
parser.e <- add_argument(parser.e, "--min_expr", help = "Minimum expression for enrichment", default = 10)
parser.e <- add_argument(parser.e, "--method", help = "Method", default = "OE")
parser.e <- add_argument(parser.e, "--enrichment_cores", help = "Number of cores to use for enrichment", default = 25)
# parse arguments
argv.e <- parse_args(parser.e, argv = args.effects)

# parameters
mscor.threshold = argv.e$mscor
padj.threshold = argv.e$fdr
modules_cutoff = argv.e$modules_cutoff
bin.size = argv.e$bin_size
min.size = argv.e$min_size
max.size = argv.e$max_size
min.expr = argv.e$min_expr
method = argv.e$method
enrichment_cores = argv.e$enrichment_cores
split_training = argv.e$train

# load SPONGE data
print("loading SPONGE data...")
load(argv.e$spongeData)

# backend
num.of.cores <- argv.e$cpus
cat("registering back end with", num.of.cores, "cores\n")
cl <- makeCluster(num.of.cores) 
registerDoParallel(cl)

# meta data
meta <- read.csv(file = argv.e$meta, sep = "\t")

# split between train and test
n.train <- round(nrow(meta)*argv.e$train)
n.test <- round(nrow(meta)-n.train)
cat("using", n.train, "samples for training\n", "and", n.test, "for testing\n")
# take n samples from each group TODO: randomize samples
cond.split <- split(meta, meta$condition)
cond.ratio <- round(sapply(cond.split, function(x) nrow(x)*argv.e$train))
train.meta <- data.table::rbindlist(mapply(function(x,y) x[1:y,], cond.split, cond.ratio, SIMPLIFY = F))
colnames(train.meta)<-c("sampleIDs",c(colnames(train.meta)[3:length(train.meta)-1],"label"))
test.meta <- data.table::rbindlist(mapply(function(x,y) x[(y+1):nrow(x),], cond.split, cond.ratio, SIMPLIFY = F))
colnames(test.meta)<-c("sampleIDs",c(colnames(test.meta)[3:length(test.meta)],"label"))

# train gene expression; argv.e$train of samples
train_gene_expr <- gene_expr[rownames(gene_expr) %in% train.meta$sampleIDs,]
train_gene_expr<-t(train_gene_expr)
# test gene expression; rest of samples
test_gene_expr <- gene_expr[rownames(gene_expr) %in% test.meta$sampleIDs,]
test_gene_expr<-t(test_gene_expr)

# train miRNA expression
train_mirna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% train.meta$sampleIDs,]
# test miRNA expression
test_mirna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% test.meta$sampleIDs,]

# ceRNA interactions
ceRNA_interactions_fdr <- ceRNA_interactions_sign[which(ceRNA_interactions_sign$p.adj < argv.e$fdr),]
n.train.ce <- round(nrow(ceRNA_interactions_fdr)*argv.e$train)
n.test.ce <- nrow(ceRNA_interactions_fdr)-n.train.ce
train_ceRNA_interactions <- head(ceRNA_interactions_fdr, n.train.ce)
test_ceRNA_interactions <- tail(ceRNA_interactions_fdr, n.test.ce)
# network centralities
network_centralities  <- sponge_node_centralities(ceRNA_interactions_fdr)
n.train.c <- round(nrow(network_centralities)*argv.e$train)
n.test.c <- nrow(network_centralities)-n.train.c
train_network_centralities <- head(network_centralities, n.train.c)
test_network_centralities <- tail(network_centralities, n.test.c)


# spongEffects filtering
print("filtering centralities...")
filtered_network_centralities <- filter_ceRNA_network(sponge_effects = train_ceRNA_interactions, 
                                                      Node_Centrality = train_network_centralities,
                                                      add_weighted_centrality=T, 
                                                      mscor.threshold = mscor.threshold, 
                                                      padj.threshold = padj.threshold)

# RNAs of interest
RNAs <- c("lncRNA","circRNA")
RNAs.ofInterest <- ensembl.df %>% dplyr::filter(gene_biotype %in% RNAs) %>%
  dplyr::select(ensembl_gene_id)
# add circRNAs of the data set
new_circRNAs <- data.frame(ensembl_gene_id=colnames(gene_expr)[grep("c", colnames(gene_expr))])
RNAs.ofInterest <- rbind(RNAs.ofInterest, new_circRNAs)
central_gene_modules<-get_central_modules(central_nodes = RNAs.ofInterest$ensembl_gene_id,
                                          node_centrality = filtered_network_centralities$Node_Centrality,
                                          ceRNA_class = RNAs, 
                                          centrality_measure = "Weighted_Degree", 
                                          cutoff = modules_cutoff)

# set modules
Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, 
                                 central.modules = central_gene_modules, 
                                 remove.central = F, 
                                 set.parallel = F)
# Module size distribution
Size.modules <- sapply(Sponge.modules, length)

# train and test modules
train.modules <- enrichment_modules(Expr.matrix = train_gene_expr, modules = Sponge.modules, bin.size = bin.size, min.size = min.size, max.size = max.size, min.expr = min.expr, method = method, cores=enrichment_cores)
test.modules <-  enrichment_modules(Expr.matrix = test_gene_expr, modules = Sponge.modules, bin.size = bin.size, min.size = min.size, max.size = max.size, min.expr = min.expr, method = method, cores=enrichment_cores)

# We find modules that were identified both in the train and test and use those as input features for the model
common_modules = intersect(rownames(train.modules), rownames(test.modules))
train.modules = train.modules[common_modules, ]
test.modules = test.modules[common_modules, ]
trained.model = calibrate_model(Input = train.modules, modules_metadata = train.meta, label = "label", sampleIDs = "sampleIDs",Metric = "Exact_match", n_folds = 10, repetitions = 3)
trained.model[["ConfusionMatrix_training"]]

Input.test <- t(test.modules) %>% scale(center = T, scale = T)
Prediction.model <- predict(trained.model$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test.meta$label))
trained.model$ConfusionMatrix_testing<-ConfusionMatrix_testing

# Define random modules
Random.modules <- Random_spongEffects(sponge_modules = Sponge.modules,
                                      gene_expr = train_gene_expr, min.size = min.size,bin.size = bin.size, max.size = max.size,
                                      min.expression=min.expr, replace = F,method = method,cores = enrichment_cores)
# We can now use the randomly defined modules to calculate their enrichment in the test set
Random.modules.test <- enrichment_modules(Expr.matrix = test_gene_expr, modules = Random.modules$Random_Modules, bin.size = bin.size, min.size = min.size, max.size = max.size, min.expr = min.expr, method = method, cores=enrichment_cores)

# Train on randoms
# We find random modules that were identified both in the train and test and use those as input features for the model
common_modules_random = intersect(rownames(Random.modules$Enrichment_Random_Modules), rownames(Random.modules.test))
Random.modules.train = Random.modules$Enrichment_Random_Modules[common_modules_random, ]
Random.modules.test = Random.modules.test[common_modules_random, ]
Random.model = calibrate_model(Input = Random.modules.train, modules_metadata = train.meta, label = "label", sampleIDs = "sampleIDs",Metric = "Exact_match", n_folds = 2, repetitions = 1)
Random.model[["ConfusionMatrix_training"]]

# validate on test set
Input.test <- t(Random.modules.test) %>% scale(center = T, scale = T)
Input.test<-Input.test[ , apply(Input.test, 2, function(x) !any(is.na(x)))]
Prediction.model <- predict(Random.model$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing_random <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test.meta$label))
Random.model$ConfusionMatrix_testing_random<-ConfusionMatrix_testing_random
ConfusionMatrix_testing_random

# Train on central genes
Input.centralgenes.train <- train_gene_expr[rownames(train_gene_expr) %in% names(Sponge.modules), ]
Input.centralgenes.test <- test_gene_expr[rownames(test_gene_expr) %in% names(Sponge.modules), ]
common_modules = intersect(rownames(Input.centralgenes.train), rownames(Input.centralgenes.test))
Input.centralgenes.train = Input.centralgenes.train[common_modules, ]
Input.centralgenes.test = Input.centralgenes.test[common_modules, ]
# Calibrate model
CentralGenes.model = calibrate_model(Input = Input.centralgenes.train, modules_metadata = train.meta, label = "label", sampleIDs = "sampleIDs",Metric = "Exact_match", n_folds = 4, repetitions = 3)
# Validate on test set
Input.centralgenes.test <- t(Input.centralgenes.test) %>% scale(center = T, scale = T)
CentralGenes.prediction <- predict(CentralGenes.model$Model, Input.centralgenes.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(CentralGenes.prediction), as.factor(test.meta$label))
CentralGenes.model$ConfusionMatrix_testing<-ConfusionMatrix_testing
ConfusionMatrix_testing

common_modules_ALLEXPR = intersect(rownames(train_gene_expr), rownames(test_gene_expr))
Input.allexpression.train = train_gene_expr[common_modules, ]
Input.allexpression.test = test_gene_expr[common_modules, ]
# Calibrate model
AllExpression.model = calibrate_model(Input = Input.allexpression.train, modules_metadata = train.meta, label = "label", sampleIDs = "sampleIDs",Metric = "Exact_match", n_folds = 10, repetitions = 3)
# Validate on test set
Input.allexpression.test <- t(Input.allexpression.test) %>% scale(center = T, scale = T)
AllExpression.prediction <- predict(AllExpression.model$Model, Input.centralgenes.test)
# We compute the confusion metrix on the test set
AllExpression.ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(AllExpression.prediction), as.factor(test.meta$label))
AllExpression.model$ConfusionMatrix_testing<-AllExpression.ConfusionMatrix_testing
AllExpression.ConfusionMatrix_testing

# plot perfomance
metrics_plot<-plot_accuracy_sensitivity_specificity(trained_model =  trained.model,random_model = Random.model,central_genes_model = CentralGenes.model, all_expression_model = AllExpression.model,
                                                    training_dataset_name="RPs/MPs",testing_dataset_name="RPs/MPs",
                                                    subtypes=as.factor(train.meta$label))
metrics_plot
s
# plot results
lollipop_plot=plot_top_modules(trained_model=trained.model, k_modules_red = 2, k_modules = 10)
lollipop_plot

# distribution of spongEffect scores
density_plot_train=plot_density_scores(trained_model=trained.model,spongEffects = train.modules, meta_data = train.meta, label = "label", sampleIDs = "sampleIDs")
density_plot_train

# heatmap train
heatmap.train = plot_heatmaps(trained_model = trained.model,spongEffects = train.modules,
                              meta_data = train.meta, label = "label", sampleIDs = "sampleIDs",Modules_to_Plot = 2,
                              show.rownames = F, show.colnames = F)
heatmap.train

# heatmap test
heatmap.test = plot_heatmaps(trained_model = trained.model,spongEffects = test.modules,
                             meta_data = test.meta, label = "label", sampleIDs = "sampleIDs",Modules_to_Plot = 2,
                             show.rownames = F, show.colnames = F)
heatmap.test

# plot modules
plot_involved_miRNAs_to_modules(sponge_modules=Sponge.modules,
                                trained_model=trained.model,
                                gene_mirna_candidates= train_genes_miRNA_candidates,
                                k_modules = 2,
                                filter_miRNAs = 0.0,
                                bioMart_gene_symbol_columns = "hgnc_symbol",
                                bioMart_gene_ensembl = "hsapiens_gene_ensembl")
