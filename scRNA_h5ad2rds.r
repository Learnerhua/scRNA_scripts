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
    A script for converting single-cell data files from h5ad format to rds format.
  Example: scRNA_h5ad2rds -I h5ad_files -O rds_files; 
  scRNA_h5ad2rds -I h5ad_files -M sceasy -P sceasy_ -O rds_files -P_env /bin/python
  
  '
  )
)

#Required
parser$add_argument(
  "-I",
  "--input",
  dest = "input_dir",
  type = "character",
  required = FALSE,
  help = "Path to a directory containing h5ad files (suffix '.h5ad', case-insensitive). Note: the X layer of the h5ad files must contain raw count data. (Mutually exclusive with -F)"
)

parser$add_argument(
  "-F",
  "--file",
  dest = "input_file",
  type = "character",
  required = FALSE,
  help = "Path to a single h5ad file (mutually exclusive with -I)"
)
# Optional
parser$add_argument(
  '-M',
  '--method',
  dest = 'method',
  help = 'R package to use for conversion. Options: "anndataR" (default) or "sceasy".',
  type = "character",
  default = "anndataR"
)

parser$add_argument(
  '-T',
  '--threads',
  dest = 'threads',
  help = 'Number of threads used for conversion. (when using "sceasy")',
  type = "integer",
  default = 12
)

parser$add_argument(
  '-P_env',
  dest = 'python_env',
  help = 'The path to the Python environment required to use sceasy for conversion. This environment must have the "numpy" and "anndata" packages installed in advance. When the -M parameter is "sceasy", this option must be specified.',
  type = "character"
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
input_dir <- opts$input_dir
input_file <- opts$input_file
method <- opts$method
threads <- opts$threads
python_env <- opts$python_env
prefix <- opts$prefix
outpath <- opts$outpath

# -I and -F are mutually exclusive; exactly one must be provided
if (is.null(input_dir) && is.null(input_file)) {
  stop("Please provide either -I (input_dir) or -F (single file).")
}
if (!is.null(input_dir) && !is.null(input_file)) {
  stop("Please provide either -I (input_dir) or -F (single file), not both.")
}
if (!is.null(input_file) && !file.exists(input_file)) {
  stop(paste0("Input file not found: ", input_file))
}

if (method == "sceasy" & is.null(python_env)){
  stop("When using 'sceasy' method, you must specify a valid python path through '-P_env' parameter.")
}

message("\n=====================scRNA Format Conversion===================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

# 设置输出目录
# showWarnings=FALSE + 无条件调用 => 避免多进程并发时 "file already exists" 错误
dir.create(outpath, recursive = TRUE, showWarnings = FALSE)

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
})

if (method == "anndataR") {
  suppressPackageStartupMessages({
    if(!require(anndataR, quietly = TRUE)){remotes::install_github("scverse/anndataR")}
  })
} else if (method == "sceasy") {
  suppressPackageStartupMessages({
    if(!require(sceasy, quietly = TRUE)){devtools::install_github("cellgeni/sceasy")}
    if(!require(reticulate,quietly = TRUE)){install.packages("reticulate")}
    if(!require(BiocParallel,quietly = TRUE)){install.packages("BiocParallel")}
  })
  use_python(python_env)
  if ((!py_module_available("numpy")) | (!py_module_available("anndata"))){
    stop("Package numpy and anndata must be installed in advance !")
  }
  os <- Sys.info()[["sysname"]]
  if (os == "Windows") {
    register(SnowParam(workers = threads, progressbar = TRUE))
  } else {
    register(MulticoreParam(workers = threads, progressbar = TRUE))
  }
}else {
  stop("Invalid method specified. Please choose either 'anndataR' or 'sceasy'.")
}

# 3. Convert h5ad files to rds format
# Build file list: either a single file (-F) or all h5ad in directory (-I)
if (!is.null(input_file)) {
  h5ad_files <- basename(input_file)
  input_dir_eff <- dirname(input_file)
} else {
  h5ad_files <- list.files(input_dir, pattern = "\\.h5ad$", full.names = FALSE, ignore.case = TRUE)
  input_dir_eff <- input_dir
  if (length(h5ad_files) == 0) {
    stop("No h5ad files found in the specified directory.")
  }
}
if (method == "anndataR") {
  cat("\nUsing anndataR for conversion.\n")
  n <- 0
  for (h5ad_file in h5ad_files) {
    n <- n + 1
    cat("\nProcessing file: ", h5ad_file, "...\n", sep = "")
    input_path <- file.path(input_dir_eff, h5ad_file)
    rds_file <- sub("\\.h5ad$", ".rds", h5ad_file, ignore.case = TRUE)
    seurat_obj <- anndataR::read_h5ad( input_path, as = 'Seurat')
    # seurat_obj[["RNA"]] <- SetAssayData(
    #   seurat_obj[["RNA"]],
    #   layer = "counts",
    #   new.data = seurat_obj[["RNA"]]@layers$X
    # )
    seurat_obj[["RNA"]]@layers$X <- NULL
    saveRDS(seurat_obj, file = file.path(outpath, paste0(prefix, rds_file)))
    cat("\nFile ", n, ": ", prefix, rds_file, " saved successfully ! \n", sep = "")
  }
} else if (method == "sceasy") {
  cat("\nUsing sceasy for conversion.\n")
  n <- 0
  for (h5ad_file in h5ad_files) {
    n <- n + 1
    cat("\nProcessing file: ", h5ad_file, "...\n", sep = "")
    input_path <- file.path(input_dir_eff, h5ad_file)
    rds_file <- sub("\\.h5ad$", ".rds", h5ad_file, ignore.case = TRUE)
    sceasy::convertFormat(
      obj = input_path,
      from = "anndata",
      to = "seurat",
      outFile = file.path(outpath, paste0(prefix, rds_file))
    )
    cat("\nFile ", n, ": ", prefix, rds_file, " saved successfully ! \n", sep = "")
  }
}
