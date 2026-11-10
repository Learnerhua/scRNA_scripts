#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
rm(list = ls())
options(timeout = 1000)
options(future.globals.maxSize = 200 * 1024^3, future.seed = TRUE) 

# 1. Parameter settings
suppressPackageStartupMessages({
  if (!require(argparse, quietly = TRUE)) {
    install.packages("argparse")
  }
})

parser <- ArgumentParser(
  description = paste0(
    '
    A script for decontamination and doublets removal in single-cell RNA-seq data.
  Example: scRNA_dedoublets -I seo_afterQC.rds -D -PC 30 -R 0.1 -CA Leiden -P_env /bin/python3.12 -DR 0.075 -T 20 -P 03. -O 03_deContamination
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
  help = 'Path to an rds file for a seurat object used for contamination ans doublets removal. (Required)'
)

# Optional
parser$add_argument(
  '-D',
  '--decontx',
  dest = 'decontx',
  action = 'store_true',
  help = 'Apply DecontX to remove ambient RNA contamination. Specify this flag to enable; omit it to disable. (Default: disabled)',
  default = FALSE
)
parser$add_argument(
  '-PC',
  dest = 'PC_num',
  help = 'Number of PCs used for dimension reduction and clustering. (Optional, default: 30)',
  type = "integer",
  default = 30
)
parser$add_argument(
  '-R',
  dest = 'resolution',
  help = 'Resolution parameter for clustering. (Optional, default: 0.1)',
  type = "numeric",
  default = 0.1
)
parser$add_argument(
  '-CA',
  dest = 'cluster_algorithm',
  help = 'Algorithm used in "FindClusters" function, one of "Louvain" or "Leiden". If specified as "Leiden", you must configure the Python environment and use the -P_env parameter to specify the path to that Python environment.(Optional, default: Louvain)',
  type = "character",
  default = "Louvain"
)
parser$add_argument(
  '-P_env',
  dest = 'python_env',
  help = 'The path to the Python environment required to use the Leiden algorithm.',
  type = "character"
)
parser$add_argument(
  '-DR',
  dest = 'double_rate',
  help = 'Expected doublet rate (proportion of droplets containing two cells). (Optional, default: 0.075)',
  type = "numeric",
  default = 0.075
)
parser$add_argument(
  '-T',
  dest = 'threads',
  help = 'Number of threads used for dimension reduction and clustering. (Optional, default: 10)',
  type = "integer",
  default = 10
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
decontx <- opts$decontx
PC_num <- opts$PC_num
resolution <- opts$resolution
cluster_algorithm <- opts$cluster_algorithm
python_env <- opts$python_env
double_rate <- opts$double_rate
threads <- opts$threads
prefix <- opts$prefix
outpath <- opts$outpath

if (!is.null(python_env) && !file.exists(python_env)) {
  stop("The specified Python environment path does not exist.")
}

message("\n===============scRNA Decontamination & Dedoublets==============\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

# 设置输出目录
if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(celda,quietly = TRUE)){install.packages("celda")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(reticulate,quietly = TRUE)){install.packages("reticulate")}
  if(!require(DoubletFinder,quietly = TRUE)){install.packages("DoubletFinder")}
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

# 3. Load data
seo <- readRDS(input_file)
cat("\nLoaded Seurat object from:", input_file, "\n")
seo <- seo %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData() %>% 
  RunPCA(ndims = PC_num, features = VariableFeatures(object = seo))
cat("\nSeurat object normalized and PCA performed.\n")

options(future.plan = NULL)
future::plan("multicore", workers = threads)

seo <- seo %>%
  RunUMAP(dims = 1:PC_num) %>%
  FindNeighbors(dims = 1:PC_num) %>%
  FindClusters(resolution = resolution, algorithm = algorithm_num, random.seed = 123, verbose = FALSE)
cat("\nUMAP, t-SNE, neighbors, and clustering completed.\n")

# 4. Decontamination
# Based on 'celda' R package
if (decontx){
  cat("\nApplying DecontX for decontamination...\n")
  sce <- as.SingleCellExperiment(seo)
  sce.decont <- decontX(x = sce)
  seo@meta.data$decontX_contamination <- colData(sce.decont)$decontX_contamination
  seo$decontX_class <- "L-prob"
  seo$decontX_class[which(seo$decontX_contamination > quantile(seo$decontX_contamination, 0.9))] <- "H-prob"
  seo$decontX_class <- factor(seo$decontX_class, levels = c("L-prob", "H-prob"))
  saveRDS(seo, file = file.path(outpath, paste0(prefix, "seo_contx.rds")))
  cat("File : ", prefix, "seo_contx.rds saved successfully !\n", sep = "")
  
  
  p1 <- FeaturePlot(seo,features="decontX_contamination", pt.size = 1,)+
    scale_color_viridis_c()+
    ggtitle("DecontX Contamination Levels")
  p2 <- DimPlot(seo,group.by="decontX_class")+
    scale_color_manual(values = c("L-prob" = "grey40", "H-prob" = "red"))+
    ggtitle("DecontX Class (Low vs High Probability)")
  plot1 <- wrap_plots(p1, p2, ncol =2)
  ggsave(file.path(outpath, paste0(prefix, "decontX_umap.png")), plot1 ,dpi=300, width = 16, height = 8, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "decontX_umap.pdf")), plot1, width = 16, height = 8)
  cat("\nPlot : ", prefix, "decontX_umap.png(pdf) has been successfully generated !\n", sep = "")
  
  p3 <- VlnPlot(seo, group.by = "decontX_class", 
                features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
                pt.size = 0, ncol = 3,
                cols = c("L-prob" = "lightblue", "H-prob" = "orange"))
  ggsave(file.path(outpath, paste0(prefix, "decontX_violin.png")), p3, dpi=300, width = 24, height = 8, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "decontX_violin.pdf")), p3, width = 24, height = 8)
  cat("\nPlot : ", prefix, "decontX_violin.png(pdf) has been successfully generated !\n", sep = "")
}

# 5. Doublet detection and removal
sweep.res.list <- paramSweep(seo, PCs = 1:PC_num, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
png("pk_plot.png", width = 800, height = 800)
bcmvn <- find.pK(sweep.stats)
dev.off()
file.remove("pk_plot.png")
pk_best <- bcmvn %>%
  slice_max(BCmetric) %>%
  pull(pK) %>% as.character() %>% as.numeric()
homotypic.pro <- modelHomotypic(seo$seurat_clusters)
nExp_poi <- round(double_rate * nrow(seo@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.pro))
seo <- doubletFinder(seo, PCs = 1:PC_num, pN = 0.25, pK = pk_best, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
col_names <- colnames(seo@meta.data)
colnames(seo@meta.data)[length(col_names)-1] <- "Double_score"
colnames(seo@meta.data)[length(col_names)] <-  "Is_Double"
seo$Is_Double <- factor(seo$Is_Double, levels = c("Singlet", "Doublet"))
saveRDS(seo, file = file.path(outpath, paste0(prefix, "seo_doublets.rds")))
cat("File : ", prefix, "seo_doublets.rds saved successfully !\n", sep = "")

# Visulization
p4 <- FeaturePlot(seo,features="Double_score", reduction = "umap", pt.size = 1)+
  scale_color_viridis_c()+
  ggtitle("Doublet score per cell")
p5 <- DimPlot(seo, group.by="Is_Double", reduction = "umap")+
  scale_color_manual(values = c("Singlet" = "grey40", "Doublet" = "red"))+
  ggtitle("Predicted singlets and doublets")
plot2 <- wrap_plots(p4, p5, ncol =2)
ggsave(file.path(outpath, paste0(prefix, "doublet_umap.png")), plot2, dpi=300, width = 16, height = 8, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "doublet_umap.pdf")), plot2, width = 16, height = 8)
cat("\nPlot : ", prefix, "doublet_umap.png(pdf) has been successfully generated !\n", sep = "")
p6 <- VlnPlot(seo, group.by = "Is_Double", 
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), 
              pt.size = 0, ncol = 3,
              cols = c("Singlet" = "lightblue", "Doublet" = "orange"))
ggsave(file.path(outpath, paste0(prefix, "doublet_violin.png")), p6, dpi=300, width = 24, height = 8, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "doublet_violin.pdf")), p6, width = 24, height = 8)
cat("\nPlot : ", prefix, "doublet_violin.png(pdf) has been successfully generated !\n", sep = "")

# Filter out doublets
cat("\nNumber of cells per sample before contamination and doublets removal:\n")
table(seo$orig.ident)
if ("decontX_class" %in% colnames(seo@meta.data)){
  seo_filtered <- subset(seo, subset = Is_Double == "Singlet" & decontX_class == "L-prob")
  cat("\nNumber of cells per sample after contamination and doublets removal:\n")
  table(seo_filtered$orig.ident)
} else {
  seo_filtered <- subset(seo, subset = Is_Double == "Singlet")
  cat("\nNumber of cells per sample after doublets removal:\n")
  table(seo_filtered$orig.ident)
}
saveRDS(seo_filtered, file = file.path(outpath, paste0(prefix, "seo_deCont.rds")))
cat("File : ", prefix, "seo_deCont.rds saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")
