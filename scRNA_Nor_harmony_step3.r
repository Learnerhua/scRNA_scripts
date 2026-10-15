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
    A script for normalization, HVG selection, scaling, PCA and batch-effect
    correction with Harmony (based on Seurat v5). Input is the counts-only
    Seurat object from scRNA_dedoublets_step2 (seu_clean.rds).
    The batch column defaults to "sample" (7 donors) - never set it to "group",
    otherwise the biological comparison (treatment vs control) would be removed.
  Example: scRNA_Nor_harmony_step3 -I 02_deContamination/02.seu_clean.rds -PC 30 -HT 2 -T 20 -P 03. -O 03_Nor_Harmony
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
  help = 'Path to the counts-only Seurat object rds (output of scRNA_dedoublets_step2). (Required)'
)
# Optional
parser$add_argument(
  '-BC',
  '--batch_col',
  dest = 'batch_col',
  help = 'Metadata column used as the Harmony batch variable. Use the per-sample column ("sample"), NOT the disease group column. (Optional, default: sample)',
  type = "character",
  default = "sample"
)
parser$add_argument(
  '-NF',
  '--n_features',
  dest = 'n_features',
  help = 'Number of highly variable genes selected by vst. (Optional, default: 2000)',
  type = "integer",
  default = 2000
)
parser$add_argument(
  '-HT',
  '--harmony_theta',
  dest = 'Harmony_theta',
  help = 'Diversity clustering penalty for Harmony. Larger values force stronger mixing; 0 disables the penalty. (Optional, default: 2)',
  type = "numeric",
  default = 2
)
parser$add_argument(
  '-PC',
  dest = 'PC_num',
  help = 'Number of PCs used for dimension reduction. (Optional, default: 30)',
  type = "integer",
  default = 30
)
parser$add_argument(
  '-T',
  dest = 'threads',
  help = 'Number of threads for the harmony algorithm. (Optional, default: 10)',
  type = "integer",
  default = 10
)
parser$add_argument(
  "-J", "--doJackStraw",
  action = "store_true",
  default = FALSE,
  help = "Perform JackStraw analysis for PCA (default: FALSE)"
)
parser$add_argument(
  '-W',
  '--width',
  dest = 'plot_width',
  help = 'width of plots in inches (Optional, default: 12)',
  type = "numeric",
  default = 12
)
parser$add_argument(
  '-H',
  '--height',
  dest = 'plot_height',
  help = 'height of plots in inches (Optional, default: 8)',
  type = "numeric",
  default = 8
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
batch_col <- opts$batch_col
n_features <- opts$n_features
Harmony_theta <- opts$Harmony_theta
PC_num <- opts$PC_num
threads <- opts$threads
doJackStraw <- opts$doJackStraw
plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath

if (batch_col == "group") {
  stop("Error: batch_col must NOT be 'group'. Harmony batch correction would remove the disease-group signal. Use the per-sample column ('sample').")
}

message("==================scRNA Normalization & Harmony===================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2026 OYJH. All Rights Reserved.")
message("\n==================================================================\n")

# 设置输出目录
if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# Summary 文件（终端输出的结构化记录）
summary_file <- file.path(outpath, paste0(prefix, "Nor_harmony_summary.txt"))
summary_con <- file(summary_file, open = "wt")
write_summary <- function(...) {
  # 入参任一为 length>1 向量时会被 paste0 collapse 压行 → 报错提示改用 cat(sep="\n")
  args_list <- list(...)
  multi <- vapply(args_list, length, integer(1)) > 1
  if (any(multi)) {
    which_arg <- which(multi)[1]
    stop("write_summary arg #", which_arg,
         " is a length-", length(args_list[[which_arg]]),
         " vector which would be collapsed into a single line by paste0. ",
         "Use cat(..., sep=\"\\n\", file = summary_con) instead for multi-line content.")
  }
  txt <- paste0(..., collapse = "")
  cat(txt, file = summary_con)
}
write_summary("==============scRNA Normalization & Harmony Summary===============\n")
write_summary("Current time: ", as.character(Sys.time()), "\n")
write_summary("Input RDS   : ", input_file, "\n")
write_summary("Batch col   : ", batch_col, " | n_features: ", n_features,
              " | theta: ", Harmony_theta, " | PC: ", PC_num, " | threads: ", threads, "\n")

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(harmony,quietly = TRUE)){install.packages("harmony")}
})

# 3. Load data
seu <- readRDS(input_file)
cat("\nLoaded Seurat object from:", input_file, "\n")
cat("Cells:", ncol(seu), " | Genes:", nrow(seu), " | Samples:",
    length(unique(seu[[batch_col, drop = TRUE]])), "\n")
if (!(batch_col %in% colnames(seu@meta.data))) {
  stop("Error: batch column '", batch_col, "' not found in metadata. Available: ",
       paste(colnames(seu@meta.data), collapse = ", "))
}
cat("Batch variable distribution (", batch_col, "):\n", sep = "")
print(table(seu[[batch_col, drop = TRUE]]))
write_summary("\n---------- Cells per batch ----------\n")
tb_batch <- as.data.frame(table(seu[[batch_col, drop = TRUE]]), stringsAsFactors = FALSE)
colnames(tb_batch) <- c(batch_col, "cells")
cat(capture.output(write.table(tb_batch, sep = "\t", row.names = FALSE, quote = FALSE)),
    file = summary_con, sep = "\n")
cat("\n", file = summary_con)

harmony_rds <- file.path(outpath, paste0(prefix, "seu_harmony.rds"))

if (file.exists(harmony_rds)) {
  # Skip 重计算：读缓存 rds 直接重画图（改 -W/-H 时用）
  cat("\n>>> Found existing seu_harmony.rds, skipping computation and re-plotting only.\n")
  cat("    Path:", harmony_rds, "\n")
  seu <- readRDS(harmony_rds)
  if (!("harmony" %in% names(seu@reductions))) {
    stop("Existing seu_harmony.rds lacks the 'harmony' reduction - looks stale. ",
         "Delete the rds to re-run the full pipeline.")
  }
  cat("    Loaded cells:", ncol(seu), " | harmony reduction present.\n")
} else {
  # 4. Normalize / HVG / Scale / PCA
  cat("\n>>> NormalizeData (LogNormalize, scale.factor=10000) ...\n")
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

  cat(">>> FindVariableFeatures (vst, nfeatures=", n_features, ") ...\n", sep = "")
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = n_features)
  cat("    HVGs selected:", length(VariableFeatures(seu)), "\n")

  cat(">>> ScaleData ...\n")
  seu <- ScaleData(seu)

  cat(">>> RunPCA (npcs=", PC_num, ") ...\n", sep = "")
  seu <- RunPCA(seu, npcs = PC_num, verbose = FALSE)
  cat("PCA completed.\n")

  # 5. Harmony batch correction
  cat(">>> RunHarmony (group.by.vars=", batch_col, ", theta=", Harmony_theta,
      ", ncores=", threads, ") ...\n", sep = "")
  seu <- RunHarmony(seu, group.by.vars = batch_col, max_iter = 20,
                    theta = Harmony_theta, ncores = threads,
                    early_stop = TRUE, verbose = TRUE)
  cat("Harmony completed.\n")

  saveRDS(seu, file = harmony_rds)
  cat("\nFile 1 : ", prefix, "seu_harmony.rds saved successfully !\n", sep = "")
}

# 6. Visualization (always runs, both full pipeline and cache re-plot)
cat("\n>>> Generating plots ...\n")

p_pca <- DimPlot(seu, reduction = "pca", group.by = batch_col)
ggsave(file.path(outpath, paste0(prefix, "pca_dimplot.png")), p_pca,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "pca_dimplot.pdf")), p_pca,
       width = plot_width, height = plot_height)
cat("Plot : ", prefix, "pca_dimplot.png(pdf) has been successfully generated !\n", sep = "")

p_elbow <- ElbowPlot(seu, ndims = PC_num)
ggsave(file.path(outpath, paste0(prefix, "pca_elbowplot.png")), p_elbow,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "pca_elbowplot.pdf")), p_elbow,
       width = plot_width, height = plot_height)
cat("Plot : ", prefix, "pca_elbowplot.png(pdf) has been successfully generated !\n", sep = "")

# JackStraw plot for PCA (optional)
if (doJackStraw) {
  cat("\nPerforming JackStraw analysis for PCA...\n")
  seu <- seu %>%
    JackStraw(dims = PC_num, num.replicate = 100, verbose = TRUE) %>%
    ScoreJackStraw(dims = 1:PC_num)
  p_js <- JackStrawPlot(seu, dims = 1:PC_num) +
    ggtitle(paste0("JackStraw Plot (Top ", PC_num, " PCs)")) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(outpath, paste0(prefix, "pca_jackstrawplot.png")), p_js,
         width = plot_width, height = plot_height, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "pca_jackstrawplot.pdf")), p_js,
         width = plot_width, height = plot_height)
  cat("Plot : ", prefix, "pca_jackstrawplot.png(pdf) has been successfully generated !\n", sep = "")
}

top10 <- head(VariableFeatures(seu), 10)
p_hvg <- VariableFeaturePlot(seu) + NoLegend()
p_hvg_lab <- LabelPoints(p_hvg, points = top10, repel = TRUE, size = 4)
ggsave(file.path(outpath, paste0(prefix, "HVGs.png")), p_hvg_lab,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "HVGs.pdf")), p_hvg_lab,
       width = plot_width, height = plot_height)
cat("Plot : ", prefix, "HVGs.png(pdf) has been successfully generated !\n", sep = "")

p_har <- DimPlot(seu, reduction = "harmony", group.by = batch_col)
ggsave(file.path(outpath, paste0(prefix, "harmony_dimplot.png")), p_har,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "harmony_dimplot.pdf")), p_har,
       width = plot_width, height = plot_height)
cat("Plot : ", prefix, "harmony_dimplot.png(pdf) has been successfully generated !\n", sep = "")

close(summary_con)
cat("\nFile 2 : ", prefix, "Nor_harmony_summary.txt saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")
