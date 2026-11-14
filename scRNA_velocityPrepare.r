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
A script for building a velocity-ready Seurat object by aligning a full scRNA-seq Seurat object
containing spliced and unspliced assays to a reference (subset) Seurat object.
Example: scRNA_velocityPrepare -R ref_seurat.rds -F full_seurat.rds -P prefix -O output
"
))

# ---- Required arguments ----
parser$add_argument(
  "-R",
  "--ref_seurat",
  type = "character",
  required = TRUE,
  help = "Path to the reference Seurat object (RDS file) defining the target cell and gene set. Metadata and dimensional reductions will be taken from this object. (Required)" 
)

parser$add_argument(
  "-F",
  "--full_seurat",
  type = "character",
  required = TRUE,
  help = "Path to the full Seurat object with spliced and unspliced assays (RDS file). (Required)"
)

# ---- Optional arguments ----
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
ref_seurat <- opts$ref_seurat
full_seurat <- opts$full_seurat
prefix <- opts$prefix
outpath <- opts$outpath

file_n <- 0

message("\n=================scRNA Cell velocity data prepare===============\n")
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

# 2. Load libraries
suppressPackageStartupMessages({
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(anndataR,quietly = TRUE)){install.packages("anndataR")}
})

# 3. Load data
cat("Loading reference seurat objects...\n")
seu1 <- readRDS(ref_seurat)    # analysis subset (reference)
cat("Loading full seurat objects...\n")
seu2 <- readRDS(full_seurat)    # full object with spliced/unspliced
cells <- colnames(seu1)
genes <- rownames(seu1)
cell_idx <- match(cells, colnames(seu2))
gene_idx <- match(genes, rownames(seu2))
stopifnot(
  all(!is.na(cell_idx)),
  all(!is.na(gene_idx))
)

rna_counts <- GetAssayData(
  seu2,
  assay = "RNA",
  layer = "counts"
)[gene_idx, cell_idx]
spliced_counts <- GetAssayData(
  seu2,
  assay = "spliced",
  layer = "counts"
)[gene_idx, cell_idx]
unspliced_counts <- GetAssayData(
  seu2,
  assay = "unspliced",
  layer = "counts"
)[gene_idx, cell_idx]

seu_velocity <- CreateSeuratObject(
  counts = rna_counts,
  meta.data = seu1@meta.data
)

stopifnot(
  identical(colnames(seu_velocity), colnames(seu1)),
  identical(rownames(seu_velocity), rownames(seu1))
)

seu_velocity[["spliced"]]   <- CreateAssayObject(counts = spliced_counts)
seu_velocity[["unspliced"]] <- CreateAssayObject(counts = unspliced_counts)

for (rd in Reductions(seu1)) {
  seu_velocity[[rd]] <- seu1[[rd]]
}

cat("\nBasic information of rebuilt seu_velocity: \n")
print(seu_velocity)


cat("\nConverting to AnnData format ......\n")
X <- t(as.matrix(GetAssayData(seu_velocity, assay = "RNA", layer = "counts")))
spliced_t <- t(as.matrix(GetAssayData(seu_velocity, assay = "spliced", layer = "counts")))
unspliced_t <- t(as.matrix(GetAssayData(seu_velocity, assay = "unspliced", layer = "counts")))
adata <- AnnData(
  X = X,
  obs = data.frame(seu_velocity@meta.data, row.names = colnames(seu_velocity)),
  var = data.frame(gene_symbols = rownames(seu_velocity), row.names = rownames(seu_velocity))
)
adata$layers[["spliced"]] <- spliced_t
adata$layers[["unspliced"]] <- unspliced_t

get_embedding <- function(seu, reduction) {
  if (!reduction %in% Reductions(seu)) {
    return(NULL)
  }
  emb <- Embeddings(seu, reduction = reduction)
  emb
}

# UMAP
umap <- get_embedding(seu1, "umap")
if (!is.null(umap)) {
  cell_order <- rownames(adata$obs)
  stopifnot(all(cell_order %in% rownames(umap)))
  adata$obsm[["X_umap"]] <- umap[cell_order, , drop = FALSE]
}
# PCA
pca <- get_embedding(seu1, "pca")
if (!is.null(pca)) {
  adata$obsm[["X_pca"]] <- pca[rownames(adata$obs), , drop = FALSE]
}
# Harmony
harmony <- get_embedding(seu1, "harmony")
if (!is.null(harmony)) {
  adata$obsm[["X_harmony"]] <- harmony[rownames(adata$obs), , drop = FALSE]
}
stopifnot(
  identical(rownames(adata$obs), colnames(seu_velocity)),
  identical(rownames(adata$obs), colnames(seu1))
)
cat("\nBasic information of rebuilt adata_velocity: \n")
print(adata)

# Saving results
cat("Saving result...\n")
saveRDS(
  seu_velocity,
  file = file.path(
    outpath,
    paste0(prefix, "scRNA_velocity_input.rds")
  )
)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, "scRNA_velocity_input.rds"), " saved successfully !\n", sep = "")
write_h5ad(adata, file.path(outpath, paste0(prefix, "scRNA_velocity_input.h5ad")), mode = "w")
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, "scRNA_velocity_input.h5ad"), " saved successfully !\n", sep = "")
message("\n=========================All Done !============================\n")