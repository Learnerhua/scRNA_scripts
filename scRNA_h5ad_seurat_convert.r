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
    A script for converting between Seurat objects (.rds) and AnnData objects (.h5ad).
    Example: h5ad_seurat_convert.r -I input.rds -O output_dir -P output_prefix
    Example: h5ad_seurat_convert.r -I input.h5ad -O output_dir -P output_prefix
    '
  )
)

# Required
parser$add_argument(
  "-I",
  "--input",
  dest = "input_file",
  type = "character",
  required = TRUE,
  help = 'Path to an rds file (Seurat object) or h5ad file (AnnData object). (Required)'
)

# Optional
parser$add_argument(
  '-O',
  '--outpath',
  dest = 'outpath',
  help = 'Output directory. (Optional, default: current directory)',
  default = getwd()
)

parser$add_argument(
  '-P',
  '--prefix',
  dest = 'prefix',
  help = 'Prefix for output files. (Optional, default: input file name without extension)',
  type = "character",
  default = NULL
)

parser$add_argument(
  '-A',
  '--assay_name',
  dest = 'assay_name',
  help = 'Name of the assay to be converted for Seurat objects. (Optional, default: RNA)',
  type = "character",
  default = "RNA"
)

parser$add_argument(
  '-X',
  '--x_mapping',
  dest = 'x_mapping',
  help = 'Data to map to the X slot in AnnData. (Optional, default: counts)',
  type = "character",
  default = "counts"
)

parser$add_argument(
  '-C',
  '--compression',
  dest = 'compression',
  help = 'Compression type for h5ad file: "none", "gzip", or "lzf". (Optional, default: gzip)',
  type = "character",
  default = "gzip"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_file <- opts$input_file
outpath <- opts$outpath
prefix <- opts$prefix
assay_name <- opts$assay_name
x_mapping <- opts$x_mapping
compression <- opts$compression

# Start timing
start_time <- Sys.time()

message("\n======================Format Conversion====================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2026 OYJH. All Rights Reserved.")
message("\n===========================================================\n")

# Set output directory
if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(Seurat, quietly = TRUE)){install.packages("Seurat")}
  if(!require(anndataR, quietly = TRUE)){
    if (!requireNamespace("BiocManager", quietly = TRUE)){
      install.packages("BiocManager")
    }
    BiocManager::install("anndataR")
  }
})

# 3. Determine file type and perform conversion
# Get file extension
file_ext <- tolower(tools::file_ext(input_file))

# Get input file basename without extension
input_basename <- tools::file_path_sans_ext(basename(input_file))

# Set output filename prefix
if (is.null(prefix)) {
  output_prefix <- input_basename
} else {
  output_prefix <- paste0(prefix, input_basename)
}

# Check if input file exists
if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# Perform conversion based on file type
if (file_ext == "h5ad") {
  # Convert h5ad to rds
  message("\nInput file type: h5ad (AnnData object)")
  message("Converting h5ad to rds (Seurat object)...\n")

  # Read h5ad file and convert directly to Seurat object
  seurat_obj <- anndataR::read_h5ad(input_file, as = 'Seurat')

  # Save as rds
  output_file <- file.path(outpath, paste0(output_prefix, ".rds"))
  saveRDS(seurat_obj, output_file)

  message("Conversion completed successfully!")
  message("Output file: ", output_file)
  message("Seurat object info:")
  print(seurat_obj)

} else if (file_ext == "rds") {
  # Convert rds to h5ad
  message("\nInput file type: rds (Seurat object)")
  message("Converting rds to h5ad (AnnData object)...\n")

  # Read rds file
  seu <- readRDS(input_file)

  # Convert to AnnData with specified assay and x_mapping
  message(paste0("Using assay: ", assay_name, ", x_mapping: ", x_mapping))
  chr_cols <- sapply(seu@meta.data, is.factor)
  seu@meta.data[, chr_cols] <- lapply(seu@meta.data[, chr_cols], as.character)
  adata <- as_AnnData(seu,
                      assay_name = assay_name,
                      x_mapping = x_mapping)

  # Save as h5ad with compression
  output_file <- file.path(outpath, paste0(output_prefix, ".h5ad"))
  message(paste0("Writing h5ad file with compression: ", compression))
  write_h5ad(adata, output_file, compression = compression, mode = "w")

  message("Conversion completed successfully!")
  message("Output file: ", output_file)
  message("AnnData object info:")
  print(adata)

} else {
  stop("Unsupported file format: ", file_ext, ". Only .h5ad and .rds files are supported.")
}

# End timing
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "auto")

message("\n===========================================================")
message("Total elapsed time: ", format(elapsed_time, digits = 3))
message("Analysis completed at: ", end_time)
message("===========================================================\n")
