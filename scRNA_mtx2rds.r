#!/usr/bin/env Rscript
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
A script for creating seurat object file based on matrix file.
Example: scRNA_mtx2rds -M test.mtx -C test_cell.tsv -G test_gene.tsv -m test_meta.txt -CC NAME -P test_ -O test; scRNA_mtx2rds -E expressionmatrix.tsv.gz -m metadata.csv
"
))
# 1.1 Required parameter
parser$add_argument(
  "-E",
  "--expr_file",
  help = "Path to an expression matrix file (tsv/tsv.gz, CSV, or similar). The file must contain gene identifiers in the FIRST column and cell names in the header row.",
  type = "character"
)
parser$add_argument(
  "-M",
  "--mtx_file",
  help = "Path to a matrix file containing non-normalized count data or normalized data (e.g. XXX.mtx).",
  type = "character"
)
parser$add_argument(
  "-C",
  "--cell_file",
  help = "Path to a file containing a cell name in each line. Number of lines must equal number of columns in the matrix.",
  type = "character"
)
parser$add_argument(
  "-G",
  "--gene_file",
  help = "Path to a file containing a gene name in each line. Number of lines must equal number of rows in the matrix.",
  type = "character"
)
parser$add_argument(
  "-m",
  "--meta_file",
  help = "Path to a metadata file containing cell information. The file may be TSV or CSV; the delimiter will be detected automatically. Must include a column with cell names.",
  type = "character",
  required = TRUE
)
parser$add_argument(
  "-CC",
  "--cell_column",
  help = "Column in meta file containing cell names. (default: NAME)",
  type = "character",
  default = "NAME"
)
# 1.2 optional parameter
parser$add_argument("-P", "--prefix",
                    help = "Prefix for output files (optional, default: none)",
                    type = "character",
                    default = "")
parser$add_argument("-O", "--outpath",
                    help = "Output directory (optional, default: current directory)",
                    type = "character",
                    default = getwd())


# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
expr_file <- opts$expr_file
mtx_file <- opts$mtx_file
cell_file <- opts$cell_file
gene_file <- opts$gene_file
meta_file <- opts$meta_file
cell_column <- opts$cell_column
prefix <- opts$prefix
outpath <- opts$outpath

valid_expr_only <- !is.null(expr_file) &&
  is.null(mtx_file)  &&
  is.null(cell_file) &&
  is.null(gene_file)

valid_mtx_triplet <- is.null(expr_file) &&
  !is.null(mtx_file)  &&
  !is.null(cell_file) &&
  !is.null(gene_file)

if (!(valid_expr_only || valid_mtx_triplet)) {
  stop(
    "Invalid argument combination.\n",
    "You must provide EITHER:\n",
    "  (1) --expr_file ONLY, OR\n",
    "  (2) --mtx_file + --cell_file + --gene_file together.\n\n",
    "Examples:\n",
    "  Rscript script.R --expr_file matrix.tsv.gz\n",
    "  Rscript script.R --mtx_file matrix.mtx --cell_file barcodes.tsv --gene_file genes.tsv\n"
  )
}

message("\n==========================Matrix to Seurat======================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

if (!dir.exists(outpath)){
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(Matrix,quietly = TRUE)){install.packages("Matrix")}
  if(!require(data.table,quietly = TRUE)){install.packages("data.table")}
})

# 3. Load data
cat("\nReading expression matrix...\n")
if (valid_expr_only){
  df <- tryCatch({
    fread(expr_file, check.names = FALSE)
  }, error = function(e) {
    cat("fread failed, trying read.csv...\n")
    read.csv(expr_file, check.names = FALSE, stringsAsFactors = FALSE)
  })
  df <- as.data.frame(df)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  expr <- df
  rownames(expr) <- make.unique(gsub("_", "-", rownames(expr)))
} else {
  expr <- readMM(mtx_file)
  cat("Reading cell and gene names...\n")
  cells <- read.csv(cell_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  genes <- read.csv(gene_file, header = FALSE, stringsAsFactors = FALSE)[, 1]
  if (ncol(expr) != length(cells)) stop("Number of columns in matrix != number of cells")
  if (nrow(expr) != length(genes)) stop("Number of rows in matrix != number of genes")
  colnames(expr) <- cells
  rownames(expr) <- genes
  rownames(expr) <- make.unique(gsub("_", "-", rownames(expr)))
}
cat("Expression matrix dimensions: ", nrow(expr), " rows x ", ncol(expr), " columns\n")

cat("\nReading metadata...\n")
meta <- fread(meta_file, check.names = FALSE)
meta <- as.data.frame(meta)
cat("Metadata dimensions: ", nrow(meta), " rows x ", ncol(meta), " columns\n")
if (grepl("TYPE", meta[1, 1], ignore.case = TRUE)){ meta <- meta[-1, ] }
rownames(meta) <- meta[[cell_column]]
meta <- meta[, !(colnames(meta) %in% cell_column), drop = FALSE]
common_cells <- intersect(colnames(expr), rownames(meta))
meta <- meta[common_cells, , drop = FALSE]
meta[is.na(meta)] <- "Unknown"

# 4. building expression matrix
cat("\nChecking expression values...\n")
if (inherits(expr, "dgCMatrix") || inherits(expr, "dgTMatrix")) {
  is_integer_like <- !any(expr@x %% 1 != 0)
} else {
  is_integer_like <- !any(as.matrix(expr) %% 1 != 0, na.rm = TRUE)
}
if (!is_integer_like) {
  cat("Detected non-integer expression values (likely normalized data).\n")
  fake_counts <- Matrix::Matrix(0, nrow = nrow(expr), ncol = ncol(expr), sparse = TRUE)
  rownames(fake_counts) <- rownames(expr)
  colnames(fake_counts) <- colnames(expr)
  cat("Fake counts matrix created.\n")
  counts_mat <- fake_counts
  data_mat <- expr
} else {
  cat("Detected integer counts (raw data).\n")
  counts_mat <- expr
  data_mat <- log1p(expr)
}

# 5.Creating Seurat object
cat("\nCreating Seurat object...\n")
seo <- CreateSeuratObject(counts = counts_mat[, common_cells], meta.data = meta)
seo <- SetAssayData(seo, assay = "RNA", layer = "data", new.data = data_mat[, common_cells])
if (!is_integer_like) {
  seo@assays$RNA@layers$counts <- NULL #only for V5
  cat("Removed counts slot since original data was non-integer.\n")
}
cat("\nBasic information of the Seurat object:\n")
print(seo)
cat("Metadata columns:\n")
print(head(seo))
cat("Feature name:\n")
print(head(rownames(seo)))
saveRDS(seo, file = file.path(outpath, paste0(prefix, "seurat.rds")))
cat("File ", prefix, "seurat.rds saved successfully !\n", sep = "")
cat("\n==========================All Done!===============================\n")
