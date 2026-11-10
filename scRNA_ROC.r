#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
rm(list = ls())
options(timeout = 1000)

# 1. Parameter settings
# ---- Define argument parser ----
suppressPackageStartupMessages({
  if (!require(argparse, quietly = TRUE)) {
    install.packages("argparse")
    library(argparse)
  }
})

parser <- ArgumentParser(description = paste0(
  "
A script for ROC visualization for specific genes of single-cell data.
Example: scRNA_ROC -I SeuratObject.rds -SC orig.ident -G group.csv -F genelist.txt -C group1,group2 -O result
"
))
# ---- Required arguments ----
parser$add_argument("-I", "--input_rds",
                    help = "Path to an RDS file containing a Seurat object for analysis.",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-SC", "--sample_column",
                    help = "Name of the column in Seurat object metadata that indicates sample identity (e.g., orig.ident)",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-G", "--group_info",
                    help = "Path to a CSV file containing group information. The file should have two columns: 'sample' and 'group'.",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-F", "--gene_file",
                    help = "Path to a file containing feature (gene) information, with one feature per line.",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-C", "--compare",
                    help = "Names of the groups to compare, separated by commas or spaces (e.g., HC,UC)"
                    , required = TRUE
)
# ---- Optional arguments ----
parser$add_argument("-GS", "--geneset",
                    help = "Whether to use gene set scoring instead of individual genes (optional, default: FALSE)",
                    action = "store_true",
                    default = FALSE
)
parser$add_argument("-W", "--plot_width",
                    help = "Width of the ROC plot (optional, default: 12)",
                    type = "numeric",
                    default = 12)
parser$add_argument("-H", "--plot_height",
                    help = "Height of the ROC plot (optional, default: 12)",
                    type = "numeric",
                    default = 12)
parser$add_argument("-P", "--prefix",
                    help = "Prefix for output files (optional, default: none)",
                    type = "character",
                    default = "")
parser$add_argument("-O", "--outpath",
                    help = "Output directory (optional, default: current directory)",
                    type = "character",
                    default = getwd())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_rds <- opts$input_rds
sample_column <- opts$sample_column
group_info <- opts$group_info
compare <- unlist(strsplit(opts$compare, "[, ]+"))
geneset <- opts$geneset
gene_file <- opts$gene_file
plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath
file_n <- 0
plot_n <- 0

message("\n========================scRNA ROC anaysis======================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(pROC,quietly = TRUE)){install.packages("pROC")}
})

# 3. Load data
# 3.1 Load seurat object
cat("\nLoading seurat object...\n")
seo <- readRDS(input_rds)
cat("Seurat object loaded. Dimensions:", dim(seo), "\n")
# 3.2 Genes check
raw_genes <- readLines(gene_file, warn = FALSE)
cat("\nTotal genes in input list:", length(raw_genes), "\n")
dup_genes <- raw_genes[duplicated(raw_genes)]
cat("Number of duplicated genes:", length(dup_genes), "\n")
if (length(dup_genes) > 0) {
  cat("Duplicated genes:\n")
  print(unique(dup_genes))
}
unique_genes <- unique(raw_genes)
genes_in_seo <- unique_genes[unique_genes %in% rownames(seo)]
genes_not_in_seo <- unique_genes[!(unique_genes %in% rownames(seo))]
cat("Number of genes found in Seurat object:", length(genes_in_seo), "\n")
cat("Number of genes NOT found in Seurat object:", length(genes_not_in_seo), "\n")
if (length(genes_not_in_seo) > 0) {
  cat("Genes NOT found in Seurat object:\n")
  print(genes_not_in_seo)
}
gene_list <- genes_in_seo
cat("Final gene list length:", length(gene_list), "\n")
# 3.3 group information
group_df <- read.csv(group_info, header = TRUE, stringsAsFactors = FALSE)
if (length(compare) != 2) {
  stop("Error: compare must contain exactly 2 group names")
}
missing_groups <- compare[!(compare %in% unique(group_df$group))]
if (length(missing_groups) > 0) {
  cat("\nError: The following groups in 'compare' are NOT found in group file:\n")
  print(missing_groups)
  stop("\nStopping because compare groups are not valid.")
}
group_df_sub <- group_df %>% filter(group %in% compare)
cell_sample <- seo[[sample_column, drop = TRUE]]
samples <- unique(cell_sample)
cat("\nNumber of samples in Seurat object:", length(samples), "\n")
cat("Number of samples in selected groups:", nrow(group_df_sub), "\n")
missing_samples <- group_df_sub$sample[!(group_df_sub$sample %in% samples)]
if (length(missing_samples) > 0) {
  cat("\nERROR: The following samples from compare groups are NOT found in the Seurat object:\n")
  print(missing_samples)
  stop("\nStopping because some compare group samples are missing in Seurat.")
} else {
  cat("All compare-group samples are found in the Seurat object.\n")
}

# 4. Calculating featuring value
samples_to_use <- group_df_sub$sample
cat("\nCalculating feature value only for samples in compare groups ......\n")
seo <- NormalizeData(seo, normalization.method = "LogNormalize", scale.factor = 10000)
if (!geneset) {
  mat <- LayerData(seo, layer = "data")
  mean_mat <- matrix(0,
                    nrow = length(samples_to_use),
                    ncol = length(gene_list))
  rownames(mean_mat) <- samples_to_use
  colnames(mean_mat) <- gene_list
  for (smp in samples_to_use) {
    cells <- names(cell_sample[cell_sample == smp])
    if (length(cells) == 0) {
      stop("No cells found in Seurat object for sample: ", smp)
    }
    mean_mat[smp, ] <- Matrix::rowMeans(mat[gene_list, cells, drop = FALSE])
  }
  feature_df <- as.data.frame(mean_mat)
  feature_df$sample <- rownames(feature_df)
  feature_df <- feature_df %>%
    left_join(group_df_sub, by = "sample")
  feature_df$group <- factor(feature_df$group)
} else {
  seo <- AddModuleScore(
    seo,
    features = list(gene_list),
    name = "ModuleScore"
  )
  score_col <- "ModuleScore1"
  mean_vec <- numeric(length(samples_to_use))
  names(mean_vec) <- samples_to_use
  for (smp in samples_to_use) {
    cells <- names(cell_sample[cell_sample == smp])
    if (length(cells) == 0) stop("No cells found in Seurat for sample:", smp)
    tmp <- FetchData(seo, vars = score_col, cells = cells)[, score_col]
    mean_vec[smp] <- mean(tmp)
  }
  feature_df <- data.frame(
    sample = names(mean_vec),
    GeneSet = mean_vec
  )
  feature_df <- feature_df %>%
    left_join(group_df_sub, by = "sample")
  feature_df$group <- factor(feature_df$group)
}

# 5. Plot ROC curve
auc_list <- list()
roc_list <- list() 
direction_list <- list()
genes_plot <- colnames(feature_df)[ ! colnames(feature_df) %in% c("sample", "group")]
for (gene in genes_plot){
  cat("Plotting ROC curve for gene ", gene, " ......\n")
  roc_obj <- roc(response = feature_df$group,
                 predictor = feature_df[[gene]],
                 levels = compare)
  
  roc_df <- coords(
    roc_obj,
    x = "all",
    ret = c("specificity", "sensitivity", "threshold"),
    transpose = FALSE
  )
  roc_df$FPR <- 1 - roc_df$specificity
  roc_df$TPR <- roc_df$sensitivity
  roc_df$gene <- gene
  roc_list[[gene]] <- roc_df

  auc_value <- round(auc(roc_obj), 3)
  auc_list[[gene]] <- auc_value
  if (roc_obj$direction == "<") {
    direction_list[[gene]] <- paste0(compare[2], "_high")
  } else if (roc_obj$direction == ">") {
    direction_list[[gene]] <- paste0(compare[1], "_high")
  } else {
    direction_list[[gene]] <- NA
  }
  # saving plot
  file_name <- file.path(outpath, paste0(prefix, "ROC_", gsub("[^A-Za-z0-9]", "_", gene), ".png"))
  png(file_name, width = plot_width, height = plot_height, res = 300, units = "in")
  plot(roc_obj, main = paste0(gene, " , AUC = ", auc_value) , col="red")
  dev.off()
  file_name <- file.path(outpath, paste0(prefix, "ROC_", gsub("[^A-Za-z0-9]", "_", gene), ".pdf"))
  pdf(file_name, width = plot_width, height = plot_height)
  plot(roc_obj, main = paste0(gene, " , AUC = ", auc_value) , col="red")
  dev.off()
  plot_n <- plot_n + 1
  cat("Plot ", plot_n, " : ", prefix, "ROC_", gsub("[^A-Za-z0-9]", "_", gene), ".png (pdf) has been generated successfully !\n", sep = "")
}

# 5. AUC recording
auc_df <- data.frame(
  AUC = unlist(auc_list),
  direction = unlist(direction_list)
) %>% arrange(desc(AUC))
high_auc_genes <- rownames(auc_df)[auc_df$AUC > 0.7]
cat("\nNumber of genes with AUC > 0.7:", length(high_auc_genes), "\n")
print(high_auc_genes)
write.csv(auc_df, file.path(outpath, paste0(prefix, "AUC_result.csv")), row.names = TRUE, quote = FALSE)
roc_all_df <- dplyr::bind_rows(roc_list)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "AUC_result.csv saved successfully !\n", sep = "")
write.csv(roc_all_df,
          file.path(outpath, paste0(prefix, "ROC_coordinates.csv")),
          row.names = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "ROC_coordinates.csv saved successfully !\n", sep = "")
message("\n=========================All Done !============================\n")

