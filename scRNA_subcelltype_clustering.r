#!/usr/bin/env Rscript
rm(list = ls())
options(timeout = 1000)

# 1. Parameter settings
suppressPackageStartupMessages({
  if (!require(argparse, quietly = TRUE)) {
    install.packages("argparse")
  }
})

parser <- ArgumentParser(
  description = paste0(
    '
    A script for single-cell data sub-celltype clustering analysis.
  Example: scRNA_subcelltype_clustering -I seo_annotated.rds -AC celltype -CA Leiden -Res resolution_recording.csv -T 12 -P test_ -O test
  '
  )
)

#Required
parser$add_argument(
  "-I",
  "--input",
  dest = "input_file",
  type = "character",
  required = TRUE,
  help = 'Path to an rds file for a seurat object used for scRNA data sub-celltype clustering. (Required)'
)
parser$add_argument(
  "-AC",
  "--annotation_column",
  dest = "annotation_column",
  type = "character",
  required = TRUE,
  help = 'Column name in metadata of the Seurat object that contains cell type annotations. (Required)'
)

# Optional
parser$add_argument(
  '-Res',
  dest = 'resolution_file',
  type = "character",
  help = 'Path to a CSV file containing cell types and their corresponding resolutions for sub-clustering. The file should have two columns: "celltype" and "resolution".'
)
parser$add_argument(
  '-PC',
  dest = 'PC_num',
  help = 'Number of PCs used for dimension reduction and clustering. (Optional, default: 30)',
  type = "integer",
  default = 30
)
parser$add_argument(
  '-CA',
  dest = 'cluster_algorithm',
  help = 'Algorithm used in "FindClusters" function, one of "Louvain" or "Leiden". If specified as "Leiden", you must configure the Python environment and use the -P_env parameter to specify the path to that Python environment.(Optional, default: Louvain)',
  type = "character",
  default = "Louvain"
)
parser$add_argument(
  '-A',
  dest = 'assay',
  help = 'The assay name that contain expression data. (Optional, default : RNA)',
  type = "character",
  default = "RNA"
)
parser$add_argument(
  '-P_env',
  dest = 'python_env',
  help = 'The path to the Python environment required to use the Leiden algorithm. When the -CA parameter is "Leiden"", it must be specified.',
  type = "character"
)
parser$add_argument(
  '-T',
  dest = 'threads',
  help = 'Number of threads used for clustering analysis. (Optional, default: 10)',
  type = "integer",
  default = 10
)
parser$add_argument(
  '-Red',
  dest = 'reduction',
  help = 'Reduction method for clustering. One of "pca" or "harmony". (Optional, default: harmony)',
  type = "character",
  default = "harmony"
)
parser$add_argument(
  '-P',
  '--prefix',
  dest = 'prefix',
  help = 'prefix for output files (Optional, default: none)',
  type = "character",
  default = ""
)
parser$add_argument(
  '-O',
  '--outpath',
  dest = 'outpath',
  help = 'output directory. (Optional, default: current directory)',
  default = getwd()
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_file <- opts$input_file
annotation_column <- opts$annotation_column
PC_num <- opts$PC_num
cluster_algorithm <- opts$cluster_algorithm
assay <- opts$assay
python_env <- opts$python_env
threads <- opts$threads
resolution_file <- opts$resolution_file
reduction <- opts$reduction
prefix <- opts$prefix
outpath <- opts$outpath

if (!is.null(python_env) && !file.exists(python_env)) {
  stop("The specified Python environment path does not exist.")
}

message("\n=================scRNA sub-celltype Clustering================\n")
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

# 2. Load required libraries
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(clustree,quietly = TRUE)){install.packages("clustree")}
  if(!require(reticulate,quietly = TRUE)){install.packages("reticulate")}
  if(!require(harmony,quietly = TRUE)){install.packages("harmony")}
})

if (cluster_algorithm == "Louvain") {
  algorithm_num <- 1
} else if (cluster_algorithm == "Leiden") {
  algorithm_num <- 4
  if (!is.null(python_env)) {
    use_python(python_env)
    reticulate::py_config()
    leiden_env <- reticulate::py_module_available("leidenalg")
    if (!leiden_env) {
      stop(
        "Error: Python module 'leidenalg' is not available.\n",
        "Please ensure the correct Python environment is set.\n",
        "Install 'leidenalg' using: pip install leidenalg\n"
      )
    }
  } else{
    reticulate::py_config()
    reticulate::py_require("leidenalg")
    leiden_env <- reticulate::py_module_available("leidenalg")
    if (!leiden_env) {
      stop("Error: can not get the leiden module !")
    }
  }
} else {
  stop("Invalid cluster_algorithm. Choose either 'Louvain' or 'Leiden'.")
}

# 3. Load the data
options(future.plan = NULL)
future::plan("multicore", workers = threads)
options(future.globals.maxSize = threads * 1024^3, future.seed = TRUE)
seo <- readRDS(input_file)
cat("\nLoaded Seurat object from:", input_file, "\n")
if (algorithm_num == 1) {
  resolutions <- seq(0, 1, by = 0.1)   # Louvain
} else {
  resolutions <- seq(0.1, 1, by = 0.1) # Leiden
}
cts <- unique(seo@meta.data[[annotation_column]])
plot_n <- 0
file_n <- 0
for (ct in cts){
  cat("\nProcessing cell type:", ct, "\n")
  if (file.exists(file.path(outpath, paste0(prefix, ct, "_clustering.rds")))){
    cat("\nClustering results for cell type ", ct, " already exist. Skipping...\n", sep = "")
    next
  }
  seo_subset <- subset(seo, celltype == ct)
  seo_subset <- NormalizeData(seo_subset)
  seo_subset <- FindVariableFeatures(seo_subset, selection.method = "vst", nfeatures = 2000)
  seo_subset <- ScaleData(seo_subset)
  seo_subset <- RunPCA(seo_subset, npcs = PC_num, verbose = FALSE)
  seo_subset <- RunHarmony(seo_subset, group.by.vars = "orig.ident", max_iter = 20, ncores = threads, early_stop = TRUE, verbose = TRUE)
  seo_subset <- FindNeighbors(seo_subset, dims = 1:PC_num, reduction = reduction)
  seo_subset <- FindClusters(seo_subset, resolution = resolutions, algorithm = algorithm_num, random.seed = 123, verbose = FALSE)
  p <- clustree(seo_subset, prefix = paste0(assay, "_snn_res.")) +
    labs(title = "Cluster trees for different resolution clustering") +
    coord_flip()
  ct <- gsub(" ", "_", ct)
  ggsave(file.path(outpath, paste0(prefix, ct, "_clustree.png")), p, width = 24, height = 16, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, ct, "_clustree.pdf")), p, width = 24, height = 16)
  plot_n <- plot_n + 1
  cat("Plot ", plot_n, " : ", prefix, ct, "_clustree.png(pdf) has been successfully generated !\n", sep = "")
  saveRDS(seo_subset, file.path(outpath, paste0(prefix,ct, "_clustering.rds")))
  file_n <- file_n + 1
  cat("File", file_n, ": ", prefix, ct, "_clustering.rds saved successfully !\n", sep = "")
}

# 4. Merge all cell types with sub-clustering results
if (is.null(resolution_file) || !file.exists(resolution_file)) {
  stop("Please provide a valid resolution file using the -Res parameter.")
}
celltype_res <- read.csv(resolution_file, header = TRUE)
cat("\nCell type resolutions for sub-clustering loaded from:", resolution_file, "\n")
cat("Cell types and their resolutions:\n")
print(celltype_res)
seo_ct_list <- list()
for (ct in cts){
  cat("\nProcessing cell type:", ct, "\n")
  celltype <- gsub(" ", "_", ct)
  seo_ct <- readRDS(file.path(outpath, paste0(prefix, celltype, "_clustering.rds")))
  res <- celltype_res[celltype_res$celltype == celltype, "resolution"]
  seo_ct <- FindClusters(seo_ct, resolution = res, algorithm = algorithm_num, random.seed = 123, verbose = FALSE)
  seo_ct$sub_celltype <- paste0(celltype, "_", seo_ct$seurat_clusters)
  Idents(seo_ct) <- "sub_celltype"
  p <- DimPlot(seo_ct, reduction = "umap", label = F, label.size = 4)
  ggsave(file.path(outpath, paste0(prefix, celltype, "_res", res,"_umap.png")), p, width = 10, height = 8, dpi = 300 , bg = "white")
  ggsave(file.path(outpath, paste0(prefix, celltype, "_res", res,"_umap.pdf")), p, width = 10, height = 8)
  plot_n <- plot_n + 1
  cat("Plot ", plot_n, " : ", prefix, celltype, "_res", res,"_umap.png(pdf) has been successfully generated !\n", sep = "")
  seo_ct_list[[celltype]] <- seo_ct
  file_n <- file_n + 1
  saveRDS(seo_ct, file.path(outpath, paste0(prefix, celltype, "_clustering_res", res, "_results.rds")))
  cat("File", file_n, ": ", prefix, celltype, "_clustering_res", res, "_results.rds saved successfully !\n", sep = "")
  
  #DE analysis for sub-celltypes
  cat("Start differential expression analysis for the sub-celltypes of ", celltype, " ...\n", sep = "")
  ct_DE_dir <- file.path(outpath, paste0(celltype,"_res", res, "_DE"))
  if (!dir.exists(ct_DE_dir)) {
    dir.create(ct_DE_dir)
  }
  Idents(seo_ct) <- "sub_celltype"
  ct_markers_df <- FindAllMarkers(seo_ct, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = TRUE)
  write.csv(ct_markers_df, file = file.path(ct_DE_dir, paste0(prefix, celltype, "_markers.csv")), row.names = FALSE, quote = FALSE)
  file_n <- file_n + 1
  cat("File", file_n, ": ", prefix, celltype, "_markers.csv saved successfully !\n", sep = "")
  # Top 10 markers for each cluster
  ct_top10_markers <- ct_markers_df %>%
    filter(p_val_adj < 0.05, pct.1 > 0.3) %>% 
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
    slice_head(n = 10)
  write.csv(ct_top10_markers, file = file.path(ct_DE_dir, paste0(prefix,celltype, "_top10_markers.csv")), row.names = FALSE, quote = FALSE)
  file_n <- file_n + 1
  cat("File", file_n, ": ", prefix, celltype, "_top10_markers.csv saved successfully !\n", sep = "")
  for (cluster in sort(unique(seo_ct@meta.data$sub_celltype))){
    cat("Generating dotplot for cluster:", cluster, "...\n")
    cluster_markers <- ct_top10_markers %>% 
      filter(cluster == !!cluster)
    cat("Number of markers for cluster", cluster, ":", nrow(cluster_markers), "\n")
    if (nrow(cluster_markers) == 0){ next }
    p <- DotPlot(seo_ct, features = cluster_markers$gene, cols = c("white", "orangered2"), dot.scale = 8, cluster.idents = T) + 
      RotatedAxis() + 
      scale_x_discrete("Markers") + 
      scale_y_discrete("sub_celltype") +
      ggtitle(paste("Cluster", cluster))
    ggsave(file.path(ct_DE_dir, paste0(prefix, cluster, "_dotplot.png")), p, width = 10, height = 10, dpi = 300, bg = "white")
    ggsave(file.path(ct_DE_dir, paste0(prefix, cluster, "_dotplot.pdf")), p, width = 10, height = 10)
    plot_n <- plot_n + 1
    cat("Plot ", plot_n, " : ", prefix, cluster, "_dotplot.png(pdf) has been successfully generated !\n", sep = "")
  }
}

# 5. Adding the sub-celltype annotation to the original Seurat object
sub_labels <- lapply(seo_ct_list, function(obj){
  data.frame(
    cell = Cells(obj),
    sub_celltype = obj$sub_celltype
  )
})
sub_labels_df <- do.call(rbind, sub_labels)
subtype_vector <- setNames(sub_labels_df$sub_celltype, sub_labels_df$cell)
seo$sub_celltype <- subtype_vector[Cells(seo)]
Idents(seo) <- "sub_celltype"
p1 <- DimPlot(seo, reduction = "umap", label = T, label.size = 4)
ggsave(file.path(outpath, paste0(prefix, "subcelltype_umap.png")), p1, width = 16, height = 12, dpi = 300 , bg = "white")
ggsave(file.path(outpath, paste0(prefix, "subcelltype_umap.pdf")), p1, width = 16, height = 12)
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "subcelltype_umap.png(pdf) has been successfully generated !\n", sep = "")
saveRDS(seo, file.path(outpath, paste0(prefix, "subcelltype.rds")))
file_n <- file_n + 1
cat("File", file_n, ": ", prefix, "subcelltype.rds saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")







