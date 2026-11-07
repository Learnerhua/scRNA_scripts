#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
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
    A script for automatically annotating single-cell data using the SingleR package.
  Example: scRNA_SingleR -I test_query.rds -RDS test_ref.rds -P test_rds. -O result -AC celltype
  scRNA_SingleR -I test_query.rds -S mouse -P test_species. -O result
  
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
  help = 'Path to an rds file for a seurat object used for scRNA data annotation. The Seurat object must contain clustering information in the "seurat_clusters" metadata. (Required)'
)

# Optional
parser$add_argument(
  "-S",
  "--species",
  dest = "species",
  type = "character",
  help = 'Specify the species for SingleR annotation using its built-in reference datasets.Options: "human" or "mouse". (Optional; default: none)'
)
parser$add_argument(
  "-RD",
  "--rdata_file",
  dest = "rdata_file",
  type = "character",
  help = 'Reference file for automatic annotation in RDATA format (SCE object). (Optional; default: none)'
)
parser$add_argument(
  "-RDS",
  "--rds_file",
  dest = "rds_file",
  type = "character",
  help = 'Reference file for automatic annotation in RDS format (Seurat object). (Optional; default: none)'
)
parser$add_argument(
  "-AC",
  "--annotation_column",
  dest = "annotation_column",
  type = "character",
  help = "The column name in the reference data that contains the annotation information. This parameter is required when using -RD or -RDS to provide a reference file. (Optional; default: none)"
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
species <- opts$species
rdata_file <- opts$rdata_file
rds_file <- opts$rds_file
annotation_column <- opts$annotation_column
prefix <- opts$prefix
outpath <- opts$outpath

num_specified <- sum(!is.null(species), !is.null(rdata_file), !is.null(rds_file))
if (num_specified == 0) {
  stop("Please specify either the species (-S) or a reference file (-RD or -RDS) for annotation.")
} else if (num_specified > 1) {
  stop("Please specify only ONE option: either species (-S) OR a reference file (-RD or -RDS), not multiple.")
}

message("\n=====================scRNA SingleR Annotation==================\n")
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
  if (!require(SingleR, quietly = TRUE)){BiocManager::install("SingleR")}
  if(!require(SingleCellExperiment,quietly = TRUE)){install.packages("SingleCellExperiment")}
  if (!require(scater, quietly = TRUE)){BiocManager::install("scater")}
})

# 3. Load input data
if (!is.null(species)){
  suppressPackageStartupMessages({
    if(!require(celldex,quietly = TRUE)){BiocManager::install("celldex")}
  })
  if (species == "human"){
    cat("\nLoading built-in human reference data from celldex package...\n")
    ref_data <- celldex::HumanPrimaryCellAtlasData()
  } else if (species == "mouse"){
    cat("\nLoading built-in mouse reference data from celldex package...\n")
    ref_data <- celldex::MouseRNAseqData()
  } else {
    stop('Invalid species specified. Please use "human" or "mouse".')
  }
}

if (!is.null(rdata_file)){
  cat("\nLoading reference data from RDATA file...\n")
  ref_data <- get(load(rdata_file))
  if (is.null(annotation_column)){
    cat("\nThe metadata information of the provided reference data:\n")
    print(colData(ref_data))
    stop("You must specify -AC/--annotation_column when using -RD for annotation.")
  }
  colData(ref_data)$label.main <- colData(ref_data)[, annotation_column]
}

if(!is.null(rds_file)){
  cat("\nLoading reference data from RDS file...\n")
  ref_seo <- readRDS(rds_file)
  if (is.null(annotation_column)){
    cat("\nThe metadata information of the provided reference data:\n")
    print(head(ref_seo@meta.data))
    stop("You must specify -AC/--annotation_column when using -RDS for annotation.")
  }
  cat("\nConverting Seurat object to SingleCellExperiment object...\n")
  counts_mat <- LayerData(ref_seo, assay = "RNA", layer = "counts")
  data_mat <- LayerData(ref_seo, assay = "RNA", layer = "data")
  if (!all(dim(counts_mat) == 0)){
    ref_data <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts_mat))
    ref_data <- scater::logNormCounts(ref_data)
  } else if (!all(dim(data_mat) == 0)){
    ref_data <- SingleCellExperiment::SingleCellExperiment(assays=list(logcounts=data_mat))
  } else {
    stop("No counts or data slot found in the Seurat object.")
  }
  colData(ref_data)$label.main <- ref_seo@meta.data[, annotation_column]
  save(ref_data, file = file.path(outpath, paste0(prefix, "SingleR_ref.rdata")))
  cat("\nThe provided reference RDS file has been successfully saved as an RData file: ", prefix, "SingleR_ref.rdata. You can directly use this RData file in subsequent analyses without conversion.\n", sep = "")
}

cat("\nLoading query data from RDS file...\n")
seo <- readRDS(file = input_file)
query_counts <- GetAssayData(object = seo@assays$RNA, layer = "counts")

#SingleR
cat("\nRunning SingleR for cell type annotation...\n")
cellpred <- SingleR(test = query_counts, ref = ref_data, clusters = seo$seurat_clusters, labels = ref_data@colData$label.main)
celltype_df <- data.frame(Cluster = rownames(cellpred), celltype = cellpred$labels, stringsAsFactors = F)
write.csv(celltype_df, file.path(outpath, paste0(prefix, "SingleR_result.csv")), row.names = F, quote = F)
cat("\nFile 1 : ", prefix, "SingleR_result.csv", " saved successfully ! \n", sep = "")
seo@meta.data$SingleR = "NA"
for (i in 1:nrow(celltype_df)){
  seo@meta.data[which(seo$seurat_clusters == celltype_df$Cluster[i]), "SingleR"] <- celltype_df$celltype[i]
}
#Visulization
p1 <- plotScoreHeatmap(cellpred, show_colnames=TRUE)
ggsave(file.path(outpath, paste0(prefix, "singleR_heatmap.png")), p1, height = 8, width = 12, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "singleR_heatmap.pdf")), p1, height = 8, width = 12)
cat("\nPlot 1 : ",  prefix, "singleR_heatmap.png(pdf) has been generated successfully !\n", sep = "")
p2 <- DimPlot(seo, reduction = "umap", group.by = "SingleR", label = TRUE, label.size = 6) + 
  labs(title = paste0("UMAP plot group by SingleR result")) +
  theme(legend.position = "right")
ggsave(file.path(outpath, paste0(prefix, "singleR_umap.png")), p2, height = 8, width = 12, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "singleR_umap.pdf")), p2, height = 8, width = 12)
cat("\nPlot 2 : ",  prefix, "singleR_umap.png(pdf) has been generated successfully !\n", sep = "")
#save RDS
saveRDS(seo, file.path(outpath, paste0(prefix, "seo_SingleR.rds")))
cat("\nFile 2 : ", prefix, "seo_SingleR.rds", " saved successfully ! \n", sep = "")

