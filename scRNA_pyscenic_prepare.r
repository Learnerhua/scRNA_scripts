#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
# =============================================================================
# Build the loom input file for pySCENIC from a Seurat object.
#
# Pull the expression matrix from one assay/layer, drop genes detected in too
# few cells, optionally subsample cells, then write a loom file whose Gene and
# CellID attributes are what the pyscenic CLI expects.
#
# Example:
#   scRNA_pyscenic_prepare.r -I seu_annotated.rds -A RNA -L data \
#       -MG 50 -NC 2000 -P 01. -O 01_pyscenic
# =============================================================================

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
    Build the loom input file for pySCENIC from a Seurat object: extract one
    assay/layer, filter genes detected in too few cells, optionally subsample
    cells, and write a loom file for the pyscenic CLI.
  Example: scRNA_pyscenic_prepare.r -I seu_annotated.rds -A RNA -L data -MG 50 -NC 2000 -P 01. -O 01_pyscenic
  '
  )
)

# Required
parser$add_argument(
  "-I",
  "--input_rds",
  dest = "input_rds",
  type = "character",
  required = TRUE,
  help = 'Path to a Seurat rds. (Required)'
)

# Optional
parser$add_argument(
  "-A",
  "--assay",
  dest = "assay",
  type = "character",
  default = "RNA",
  help = 'Assay holding the expression matrix. (Optional, default: RNA)'
)
parser$add_argument(
  "-L",
  "--layer",
  dest = "layer",
  type = "character",
  default = "data",
  help = 'Layer to use; "data" (log-normalized) is recommended, "counts" also works. (Optional, default: data)'
)
parser$add_argument(
  "-MG",
  "--min_cells",
  dest = "min_cells",
  type = "integer",
  default = 50,
  help = 'Keep genes detected (expression > 0) in at least this many cells. (Optional, default: 50)'
)
parser$add_argument(
  "-NC",
  "--n_cells",
  dest = "n_cells",
  type = "integer",
  default = 0,
  help = 'Number of cells to subsample; 0 keeps all cells. (Optional, default: 0 = all cells)'
)
parser$add_argument(
  "-SD",
  "--seed",
  dest = "seed",
  type = "integer",
  default = 123,
  help = 'Random seed for cell subsampling. (Optional, default: 123)'
)
parser$add_argument(
  "--title",
  dest = "title",
  type = "character",
  default = "pySCENIC input",
  help = 'Title written into the loom file. (Optional, default: pySCENIC input)'
)
parser$add_argument(
  "-P",
  "--prefix",
  dest = "prefix",
  type = "character",
  default = "",
  help = 'Prefix for output files. (Optional, default: none)'
)
parser$add_argument(
  "-O",
  "--outpath",
  dest = "outpath",
  type = "character",
  default = getwd(),
  help = 'Output directory. (Optional, default: current directory)'
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_rds <- opts$input_rds
assay <- opts$assay
layer <- opts$layer
min_cells <- opts$min_cells
n_cells <- opts$n_cells
seed <- opts$seed
title <- opts$title
prefix <- opts$prefix
outpath <- opts$outpath

if (!file.exists(input_rds)) {
  stop("Input file does not exist: ", input_rds)
}
if (min_cells < 1) {
  stop("-MG/--min_cells must be >= 1.")
}
if (n_cells < 0) {
  stop("-NC/--n_cells must be >= 0 (0 means keep all cells).")
}
if (!layer %in% c("data", "counts")) {
  stop('-L/--layer must be "data" or "counts"; scale.data is dense and not supported.')
}

message("\n==================scRNA pySCENIC input prepare==================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2026 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if (!require(Seurat, quietly = TRUE)) install.packages("Seurat")
  if (!require(SCopeLoomR, quietly = TRUE)) install.packages("SCopeLoomR")
  if (!require(Matrix, quietly = TRUE)) install.packages("Matrix")
})

file_n <- 0

# 3. Load data
cat(">>> Loading Seurat object\n")
seu <- readRDS(input_rds)
cat(">>> Object: ", ncol(seu), " cells x ", nrow(seu), " genes\n", sep = "")

if (!assay %in% Assays(seu)) {
  stop("Assay not found: ", assay, "; available: ", paste(Assays(seu), collapse = ", "))
}
if (!layer %in% Layers(seu[[assay]])) {
  split_layers <- Layers(seu[[assay]])
  if (any(grepl(paste0("^", layer, "\\."), split_layers))) {
    cat(">>> Assay has split ", layer, " layers (e.g. minor\'s .1 / .2); joining them\n", sep = "")
    seu <- JoinLayers(seu, assay = assay, layers = layer)
  }
}
if (!layer %in% Layers(seu[[assay]])) {
  stop("Layer not found in assay ", assay, ": ", layer,
       "; available: ", paste(Layers(seu[[assay]]), collapse = ", "))
}

# 4. Extract expression matrix
cat(">>> Extracting ", layer, " layer from assay ", assay, "\n", sep = "")
expr <- GetAssayData(seu, assay = assay, layer = layer)
if (!inherits(expr, "dgCMatrix")) {
  expr <- as(as(expr, "CsparseMatrix"), "dgCMatrix")
}

# 5. Filter genes by detection rate
n_gene_before <- nrow(expr)
expr <- expr[Matrix::rowSums(expr > 0) >= min_cells, , drop = FALSE]
n_gene_after <- nrow(expr)
cat(">>> Genes kept: ", n_gene_after, " / ", n_gene_before,
    " (detected in >= ", min_cells, " cells)\n", sep = "")
if (n_gene_after == 0) {
  stop("No gene passes the -MG threshold; lower -MG and retry.")
}

# 6. Subsample cells
n_cell_before <- ncol(expr)
if (n_cells > 0 && n_cells < n_cell_before) {
  set.seed(seed)
  cells_use <- sample(colnames(expr), n_cells)
  expr <- expr[, cells_use, drop = FALSE]
  cat(">>> Cells subsampled: ", n_cells, " / ", n_cell_before,
      " (seed ", seed, ")\n", sep = "")
} else {
  cat(">>> Cells kept: ", n_cell_before, " (no subsampling)\n", sep = "")
  if (n_cell_before > 10000) {
    cat(">>> Note: pySCENIC runtime grows with cell count; consider -NC 2000 if speed matters.\n")
  }
}

# 7. Write loom
loom_file <- paste0(prefix, "expr_norm.loom")
cat(">>> Writing loom file\n")
invisible(capture.output(
  suppressWarnings(suppressMessages(
    build_loom(
      file.name = file.path(outpath, loom_file),
      dgem = expr,
      title = title,
      loom.spec.version = 3,
      display.progress = FALSE
    )
  )),
  type = "output"
))
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", loom_file, " saved successfully !\n", sep = "")

cat("\n>>> Final matrix: ", nrow(expr), " genes x ", ncol(expr), " cells\n", sep = "")
message("\n=========================All Done !============================\n")
