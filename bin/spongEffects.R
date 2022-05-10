#!/usr/bin/env Rscript

if(!require(pacman)) install.packages("pacman", repos = "http://cran.us.r-project.org")
pacman::p_load(SPONGE, doParallel, foreach, dplyr, argparser)

args.effects = commandArgs(trailingOnly = TRUE)
# TODO add fine tuning parameters
parser.e <- arg_parser("Argument parser.e for differenial expression analysis", name = "DE_parser.e")
parser.e <- add_argument(parser.e, "--spongeData", help = "SPONGE Rdata containing all results, e.g. gene expression, miRNA expression, sponge centralities etc.")
parser.e <- add_argument(parser.e, "--meta", help = "Meta data for samples in tsv")
parser.e <- add_argument(parser.e, "--train", help = "Percentage of samples to use for training, rest will be used for testing", default = 0.8)
parser.e <- add_argument(parser.e, "--fdr", help = "False discovery rate to use for filtering SPONGE results", default = 0.05)
parser.e <- add_argument(parser.e, "--cpus", help = "Number of cores to use for backend", default = 25)

argv.e <- parse_args(parser.e, argv = args.effects)

argv.e$spongeData <- "/nfs/proj/Bongiovanni-KrdIsar-platelets/circRNA_SPONGE_spongEffects/platelets/output/results/sponging/SPONGE_manual/sponge_0.05fdr.RData"
argv.e$meta <- "/nfs/proj/Bongiovanni-KrdIsar-platelets/circRNA_SPONGE_spongEffects/platelets/resources/samplesheet.tsv"

# load SPONGE data
print("loading SPONGE data...")
load(argv.e$spongeData)

# backend
cat("registering back end with", argv.e$cpus, "cores\n")
num.of.cores <- argv.e$cpus
cl <- makeCluster(num.of.cores) 
registerDoParallel(cl)

# meta data
meta <- read.csv(file = argv.e$meta, sep = "\t")
# update colnames
colnames(meta)[which("sample" == colnames(meta))] <- "sampleID"
colnames(meta)[which("condition" == colnames(meta))] <- "SUBTYPE"
meta$sample <- gsub("-", ".", meta$sample)
# split between train and test
n.train <- round(nrow(meta)*argv.e$train)
n.test <- round(nrow(meta)-n.train)
cat("using", n.train, "samples for training\n", "and", n.test, "for testing\n")
# take n samples from each group
cond.split <- split(meta, meta$condition)
cond.ratio <- round(sapply(cond.split, function(x) nrow(x)*argv.e$train))
train.meta <- data.table::rbindlist(mapply(function(x,y) x[1:y,], cond.split, cond.ratio, SIMPLIFY = F))
colnames(train.meta)<-c("sampleID",c(colnames(train.meta)[3:length(train.meta)-1],"SUBTYPE"))
test.meta <- data.table::rbindlist(mapply(function(x,y) x[(y+1):nrow(x),], cond.split, cond.ratio, SIMPLIFY = F))
colnames(test.meta)<-c("sampleID",c(colnames(test.meta)[3:length(test.meta)],"SUBTYPE"))

# train gene expression; argv.e$train of samples
gene_expr <- t(gene_expr)
train_gene_expr <- gene_expr[rownames(gene_expr) %in% train.meta$sampleID,]
# test gene expression; rest of samples
test_gene_expr <- gene_expr[rownames(gene_expr) %in% test.meta$sampleID,]

# train miRNA expression
mi_rna_expr <- t(mi_rna_expr)
train_mirna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% train.meta$sampleID,]
# test miRNA expression
test_mirna_expr <- mi_rna_expr[rownames(mi_rna_expr) %in% test.meta$sampleID,]

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
                                                      mscor.threshold = 0.1, 
                                                      padj.threshold = 0.05)

# RNAs of interest
RNAs <- c("lncRNA","circRNA")
RNAs.ofInterest <- ensembl.df %>% dplyr::filter(gene_biotype %in% RNAs) %>%
  dplyr::select(ensembl_gene_id)
central_gene_modules<-get_central_modules(central_nodes = RNAs.ofInterest$ensembl_gene_id,
                                          node_centrality = filtered_network_centralities$Node_Centrality,
                                          ceRNA_class = RNAs, 
                                          centrality_measure = "Weighted_Degree", 
                                          cutoff = 750)

# set modules
Sponge.modules <- define_modules(network = filtered_network_centralities$Sponge.filtered, 
                                 central.modules = central_gene_modules, 
                                 remove.central = F, 
                                 set.parallel = F)
# Module size distribution
Size.modules <- sapply(Sponge.modules, length)

# train and test modules
# TODO: correct expression matrices?
train.modules <- enrichment_modules(Expr.matrix = train_gene_expr, modules = Sponge.modules, bin.size = 100, min.size = 10, max.size = 2000, min.expr = 10, method = "OE", cores=25)
test.modules <-  enrichment_modules(Expr.matrix = test_gene_expr, modules = Sponge.modules, bin.size = 10, min.size = 1, max.size = 2000, min.expr = 1, method = "OE", cores=25)

# We find modules that were identified both in the train and test and use those as input features for the model
common_modules = intersect(rownames(train.modules), rownames(test.modules))
train.modules = train.modules[common_modules, ]
test.modules = test.modules[common_modules, ]
trained.model = calibrate_model(Input = train.modules, modules_metadata = train.meta, label = "SUBTYPE", sampleIDs = "sampleID",Metric = "Exact_match", n_folds = 10, repetitions = 3)
trained.model[["ConfusionMatrix_training"]]

Input.test <- t(test.modules) %>% scale(center = T, scale = T)
Prediction.model <- predict(trained.model$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test_cancer_metadata$SUBTYPE))
trained.model$ConfusionMatrix_testing<-ConfusionMatrix_testing

# Define random modules
Random.modules <- Random_spongEffects(sponge_modules = Sponge.modules,
                                      gene_expr = train_cancer_gene_expr, min.size = 1,bin.size = 10, max.size = 200,
                                      min.expression=1, replace = F,method = "OE",cores = 1)
# We can now use the randomly defined modules to calculate their enrichment in the test set
Random.modules.test <- enrichment_modules(Expr.matrix = test_cancer_gene_expr, modules = Random.modules$Random_Modules, bin.size = 10, min.size = 1, max.size = 2000, min.expr = 1, method = "OE", cores=1)

# Train on randoms
# We find random modules that were identified both in the train and test and use those as input features for the model
common_modules_random = intersect(rownames(Random.modules$Enrichment_Random_Modules), rownames(Random.modules.test))
Random.modules.train = Random.modules$Enrichment_Random_Modules[common_modules_random, ]
Random.modules.test = Random.modules.test[common_modules_random, ]
Random.model = calibrate_model(Input = Random.modules.train, modules_metadata = train_cancer_metadata, label = "SUBTYPE", sampleIDs = "sampleID",Metric = "Exact_match", n_folds = 2, repetitions = 1)
Random.model[["ConfusionMatrix_training"]]

# validate on test set
Input.test <- t(Random.modules.test) %>% scale(center = T, scale = T)
Input.test<-Input.test[ , apply(Input.test, 2, function(x) !any(is.na(x)))]
Prediction.model <- predict(Random.model$Model, Input.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing_random <- caret::confusionMatrix(as.factor(Prediction.model), as.factor(test_cancer_metadata$SUBTYPE))
Random.model$ConfusionMatrix_testing_random<-ConfusionMatrix_testing_random
ConfusionMatrix_testing_random

# Train on central genes
Input.centralgenes.train <- train_cancer_gene_expr[rownames(train_cancer_gene_expr) %in% names(Sponge.modules), ]
Input.centralgenes.test <- test_cancer_gene_expr[rownames(test_cancer_gene_expr) %in% names(Sponge.modules), ]
common_modules = intersect(rownames(Input.centralgenes.train), rownames(Input.centralgenes.test))
Input.centralgenes.train = Input.centralgenes.train[common_modules, ]
Input.centralgenes.test = Input.centralgenes.test[common_modules, ]
# Calibrate model
CentralGenes.model = calibrate_model(Input = Input.centralgenes.train, modules_metadata = train_cancer_metadata, label = "SUBTYPE", sampleIDs = "sampleID",Metric = "Exact_match", n_folds = 1, repetitions = 1)
# Validate on test set
Input.centralgenes.test <- t(Input.centralgenes.test) %>% scale(center = T, scale = T)
CentralGenes.prediction <- predict(CentralGenes.model$Model, Input.centralgenes.test)
# We compute the confusion metrix on the test set
ConfusionMatrix_testing <- caret::confusionMatrix(as.factor(CentralGenes.prediction), as.factor(test_cancer_metadata$SUBTYPE))
CentralGenes.model$ConfusionMatrix_testing<-ConfusionMatrix_testing
ConfusionMatrix_testing

# plot perfomance
plot_accuracy_sensitivity_specificity(trained.model,CentralGenes.model,Random.model,
                                      training_dataset_name="TCGA",testing_dataset_name="TCGA",
                                      as.factor(test_cancer_metadata$SUBTYPE))

# plot results
lollipop_plot=plot_top_modules(trained_model=trained.model, k_modules_red = 2, k_modules = 4)
lollipop_plot

# distribution of spongEffect scores
density_plot_train=plot_density_scores(trained_model=trained.model,spongEffects = train.modules, meta_data = train_cancer_metadata, label = "SUBTYPE", sampleIDs = "sampleID")
density_plot_train

# heatmap train
heatmap.train = plot_heatmaps(trained_model = trained.model,spongEffects = train.modules,
                              meta_data = train_cancer_metadata, label = "SUBTYPE", sampleIDs = "sampleID",Modules_to_Plot = 2,
                              show.rownames = F, show.colnames = F)
heatmap.train

# heatmap test
heatmap.test = plot_heatmaps(trained_model = trained.model,spongEffects = test.modules,
                             meta_data = test_cancer_metadata, label = "SUBTYPE", sampleIDs = "sampleID",Modules_to_Plot = 2,
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
