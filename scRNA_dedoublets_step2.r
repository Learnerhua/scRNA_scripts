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
    A script for doublets removal (scDblFinder) and optional ambient-RNA
    decontamination (DecontX, per-sample, off by default) in single-cell data.
    Input is the merged QC-passed Seurat object from scRNA_QC_step1 (counts only).
    A lightweight Normalize/HVG/Scale/PCA/cluster pass is run internally to aid
    doublet modelling; these normalized values are NOT saved - the output object
    keeps raw counts only for the downstream normalization step.
  Example: scRNA_dedoublets_step2 -I 01_QC/01.seu_afterQC.rds -PC 30 -R 0.1 -T 20 -P 02. -O 02_deContamination
  Example (with DecontX): scRNA_dedoublets_step2 -I 01_QC/01.seu_afterQC.rds -D -PC 30 -R 0.1 -T 20 -P 02. -O 02_deContamination
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
  help = 'Path to an rds file for the merged QC-passed Seurat object (output of scRNA_QC_step1). (Required)'
)
# Optional
parser$add_argument(
  '-D',
  '--decontx',
  dest = 'decontx',
  action = 'store_true',
  help = 'Apply DecontX per sample to remove ambient RNA contamination (off by default; enable with -D).'
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
  help = 'Algorithm used in "FindClusters" function, one of "Louvain" or "Leiden". If "Leiden", you must provide a Python environment with leidenalg via -P_env. (Optional, default: Louvain)',
  type = "character",
  default = "Louvain"
)
parser$add_argument(
  '-P_env',
  dest = 'python_env',
  help = 'Path to the Python environment required to use the Leiden algorithm.',
  type = "character"
)
parser$add_argument(
  '-DR',
  dest = 'double_rate',
  help = 'Expected doublet rate passed to scDblFinder. If omitted, scDblFinder estimates it automatically from cell count. (Optional)',
  type = "numeric",
  default = NULL
)
parser$add_argument(
  '-SC',
  dest = 'sample_column',
  help = 'Column name in the Seurat metadata used as the sample identifier for per-sample doublet modelling. (Optional, default: sample)',
  type = "character",
  default = "sample"
)
parser$add_argument(
  '-T',
  dest = 'threads',
  help = 'Number of threads for scDblFinder (MulticoreParam). (Optional, default: 10)',
  type = "integer",
  default = 10
)
parser$add_argument(
  '-W',
  '--width',
  dest = 'plot_width',
  help = 'width of plots in inches (Optional, default: 16)',
  type = "numeric",
  default = 16
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
decontx <- opts$decontx
PC_num <- opts$PC_num
resolution <- opts$resolution
cluster_algorithm <- opts$cluster_algorithm
python_env <- opts$python_env
double_rate <- opts$double_rate
sample_column <- opts$sample_column
threads <- opts$threads
plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath

if (!is.null(python_env) && !file.exists(python_env)) {
  stop("The specified Python environment path does not exist.")
}
if (!(cluster_algorithm %in% c("Louvain", "Leiden"))) {
  stop("Invalid cluster_algorithm. Choose either 'Louvain' or 'Leiden'.")
}

message("\n=========================scRNA Dedoublets=========================\n")
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
summary_file <- file.path(outpath, paste0(prefix, "Dedoublets_summary.txt"))
summary_con <- file(summary_file, open = "wt")
write_summary <- function(...) {
  # 检测：任一入参是长度 >1 的向量时会被 paste0 collapse 压成一行。
  # 多行表格请改用 cat(..., sep="\n") 直接写 summary_con。
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
write_summary("=====================scRNA Dedoublets Summary=====================\n")
write_summary("Current time: ", as.character(Sys.time()), "\n")
write_summary("Input RDS  : ", input_file, "\n")
write_summary("DecontX    : ", if (decontx) "enabled (per-sample)" else "disabled", "\n")
write_summary("Sample col : ", sample_column, " | PC: ", PC_num,
              " | resolution: ", resolution, " | algorithm: ", cluster_algorithm, "\n")
write_summary("scDblFinder dbr: ", if (is.null(double_rate)) "auto" else as.character(double_rate),
              " | threads: ", threads, "\n")

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(scDblFinder,quietly = TRUE)){BiocManager::install("scDblFinder")}
  if(!require(BiocParallel,quietly = TRUE)){BiocManager::install("BiocParallel")}
  if(decontx && !require(celda,quietly = TRUE)){BiocManager::install("celda")}
  if(decontx && !require(reticulate,quietly = TRUE)){install.packages("reticulate")}
})

# Leiden 需要 python 环境
if (cluster_algorithm == "Leiden") {
  algorithm_num <- 4
  if (!is.null(python_env)) {
    use_python(python_env)
    reticulate::py_config()
    leiden_env <- reticulate::py_module_available("leidenalg")
    if (!leiden_env) {
      stop("Error: Python module 'leidenalg' is not available in the specified environment.\n",
           "Install it with: pip install leidenalg")
    }
  } else {
    reticulate::py_config()
    reticulate::py_require("leidenalg")
    leiden_env <- reticulate::py_module_available("leidenalg")
    if (!leiden_env) {
      stop("Error: cannot get the leiden module. Provide -P_env pointing to a Python env with leidenalg.")
    }
  }
} else {
  algorithm_num <- 1
}

# 3. Load data
seu <- readRDS(input_file)
cat("\nLoaded Seurat object from:", input_file, "\n")
cat("Cells:", ncol(seu), " | Genes:", nrow(seu), " | Samples:", 
    length(unique(seu[[sample_column, drop = TRUE]])), "\n")
if (!(sample_column %in% colnames(seu@meta.data))) {
  stop("Error: sample column '", sample_column, "' not found in metadata. Available: ",
       paste(colnames(seu@meta.data), collapse = ", "))
}

doublets_rds <- file.path(outpath, paste0(prefix, "seu_doublets.rds"))

if (file.exists(doublets_rds)) {
  # Skip 重计算，直接读 rds 重画图（与- W / / -H 参数对齐）
  cat("\n>>> Found existing seu_doublets.rds, skipping analysis and re-plotting only.\n")
  cat("    Path:", doublets_rds, "\n")
  seu <- readRDS(doublets_rds)
  cat("    Loaded cells:", ncol(seu), " | Genes:", nrow(seu), "\n")
  # 当 rds 来自旧版本脚本（没有 Double_score 列），提示用户先全跑一次
  if (!("Is_Double" %in% colnames(seu@meta.data))) {
    stop("Existing seu_doublets.rds lacks 'Is_Double' column — looks stale. ",
         "Delete the rds to re-run the full pipeline.")
  }
  cat("\nLoaded cached Seurat object. Skipping Normalize/PCA/UMAP/DecontX/scDblFinder.\n")
} else {
  # 4. Lightweight preprocessing (for clustering only; normalized values are NOT saved)
  cat("\n>>> Lightweight preprocessing for doublet modelling ...\n")
  seu <- seu %>%
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData() %>%
    RunPCA(npcs = PC_num, verbose = FALSE)
  cat("PCA completed.\n")

  seu <- seu %>%
    RunUMAP(dims = 1:PC_num, verbose = FALSE) %>%
    FindNeighbors(dims = 1:PC_num, verbose = FALSE) %>%
    FindClusters(resolution = resolution, algorithm = algorithm_num,
                 random.seed = 123, verbose = FALSE)
  cat("UMAP / clustering completed (clusters:", length(unique(seu$seurat_clusters)), ").\n")

  # 5. DecontX (per-sample, optional)
  if (decontx) {
    cat("\n>>> Running DecontX per sample ...\n")
    contam_vec <- rep(NA_real_, ncol(seu))
    names(contam_vec) <- colnames(seu)
    hprob_cells <- c()
    for (s in unique(seu[[sample_column, drop = TRUE]])) {
      cells <- colnames(seu)[seu[[sample_column, drop = TRUE]] == s]
      cat("  DecontX on sample:", s, "(cells:", length(cells), ") ...\n")
      sce_s <- as.SingleCellExperiment(seu[, cells])
      sce.decont <- decontX(x = sce_s)
      contam_vec[cells] <- colData(sce.decont)$decontX_contamination
    }
    seu$decontX_contamination <- contam_vec
    # Per-sample 90% quantile threshold → L-prob / H-prob class
    seu$decontX_class <- "L-prob"
    for (s in unique(seu[[sample_column, drop = TRUE]])) {
      idx <- seu[[sample_column, drop = TRUE]] == s
      thr <- quantile(seu$decontX_contamination[idx], 0.9, na.rm = TRUE)
      seu$decontX_class[idx & seu$decontX_contamination > thr] <- "H-prob"
      cat("  Sample", s, "| contamination 90% threshold:", round(thr, 4),
          "| H-prob cells:", sum(idx & seu$decontX_contamination > thr), "\n")
    }
    seu$decontX_class <- factor(seu$decontX_class, levels = c("L-prob", "H-prob"))

    p1 <- FeaturePlot(seu, features = "decontX_contamination", pt.size = 0.5) +
      scale_color_viridis_c() +
      ggtitle("DecontX Contamination Levels")
    p2 <- DimPlot(seu, group.by = "decontX_class") +
      scale_color_manual(values = c("L-prob" = "grey40", "H-prob" = "red")) +
      ggtitle("DecontX Class (Low vs High Probability)")
    plot1 <- wrap_plots(p1, p2, ncol = 2)
    ggsave(file.path(outpath, paste0(prefix, "decontX_umap.png")), plot1,
           width = plot_width, height = plot_height, dpi = 300, bg = "white")
    ggsave(file.path(outpath, paste0(prefix, "decontX_umap.pdf")), plot1,
           width = plot_width, height = plot_height)
    cat("\nPlot : ", prefix, "decontX_umap.png(pdf) has been successfully generated !\n", sep = "")

    p3 <- VlnPlot(seu, group.by = "decontX_class",
                  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                  pt.size = 0, ncol = 3,
                  cols = c("L-prob" = "lightblue", "H-prob" = "orange"))
    ggsave(file.path(outpath, paste0(prefix, "decontX_violin.png")), p3,
           width = plot_width, height = plot_height, dpi = 300, bg = "white")
    ggsave(file.path(outpath, paste0(prefix, "decontX_violin.pdf")), p3,
           width = plot_width, height = plot_height)
    cat("Plot : ", prefix, "decontX_violin.png(pdf) has been successfully generated !\n", sep = "")
  } else {
    cat("\nDecontX disabled (skip).\n")
  }

  # 6. Doublet detection (scDblFinder, per-sample via samples argument)
  cat("\n>>> Running scDblFinder ... Current time:", format(Sys.time()), "\n")
  set.seed(123)
  sce_dbl <- as.SingleCellExperiment(seu)
  sce_dbl <- scDblFinder(sce_dbl,
                         dbr = double_rate,
                         samples = sample_column,
                         BPPARAM = MulticoreParam(threads))
  cat("scDblFinder completed. Current time:", format(Sys.time()), "\n")

  seu$Double_score <- colData(sce_dbl)$scDblFinder.score
  seu$Is_Double    <- colData(sce_dbl)$scDblFinder.class
  seu$Is_Double    <- ifelse(seu$Is_Double == "doublet", "Doublet", "Singlet")
  seu$Is_Double    <- factor(seu$Is_Double, levels = c("Singlet", "Doublet"))
  cat("Doublets detected:", sum(seu$Is_Double == "Doublet"), "/", ncol(seu), "cells\n")
  write_summary("\n---------- Doublet summary (per sample) ----------\n")
  tb_dbl <- as.data.frame.matrix(table(seu[[sample_column, drop = TRUE]], seu$Is_Double))
  tb_dbl$sample <- rownames(tb_dbl)
  tb_dbl <- tb_dbl[, c("sample", "Singlet", "Doublet")]
  cat(capture.output(write.table(tb_dbl, sep = "\t", row.names = FALSE, quote = FALSE)),
      file = summary_con, sep = "\n")
  cat("\n", file = summary_con)
  write_summary("\n")
  # 保存带标记的对象（未过滤）
  saveRDS(seu, file = doublets_rds)
  cat("\nFile 1 : ", prefix, "seu_doublets.rds saved successfully !\n", sep = "")
}

# 7. Visualization (always, regardless of skip vs full pipeline)
p4 <- FeaturePlot(seu, features = "Double_score", reduction = "umap", pt.size = 0.5) +
  scale_color_viridis_c() +
  ggtitle("Doublet score per cell")
p5 <- DimPlot(seu, group.by = "Is_Double", reduction = "umap") +
  scale_color_manual(values = c("Singlet" = "grey40", "Doublet" = "red")) +
  ggtitle("Predicted singlets and doublets")
plot2 <- wrap_plots(p4, p5, ncol = 2)
ggsave(file.path(outpath, paste0(prefix, "doublet_umap.png")), plot2,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "doublet_umap.pdf")), plot2,
       width = plot_width, height = plot_height)
cat("\nPlot : ", prefix, "doublet_umap.png(pdf) has been successfully generated !\n", sep = "")

p6 <- VlnPlot(seu, group.by = "Is_Double",
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              pt.size = 0, ncol = 3,
              cols = c("Singlet" = "lightblue", "Doublet" = "orange"))
ggsave(file.path(outpath, paste0(prefix, "doublet_violin.png")), p6,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "doublet_violin.pdf")), p6,
       width = plot_width, height = plot_height)
cat("Plot : ", prefix, "doublet_violin.png(pdf) has been successfully generated !\n", sep = "")

# 8. Filter out doublets (+ H-prob contamination when DecontX enabled)
#    统计信息（before/after 细胞数）总是计算并写入 summary；
#    仅当 clean rds 缺失时才重新过滤并保存 rds（缓存与 clean rds 均已存在时
#    不重复保存，画完图即退出）。
clean_rds <- file.path(outpath, paste0(prefix, "seu_clean.rds"))

cat("\nNumber of cells per sample before doublets removal:\n")
print(table(seu[[sample_column, drop = TRUE]]))
write_summary("\n---------- Cells per sample: before / after ----------\n")
if (decontx) {
  keep <- seu$Is_Double == "Singlet" & seu$decontX_class == "L-prob"
} else {
  keep <- seu$Is_Double == "Singlet"
}
seu_clean <- subset(seu, cells = colnames(seu)[keep])
cat("Number of cells per sample after doublets removal:\n")
print(table(seu_clean[[sample_column, drop = TRUE]]))
for (s in unique(seu[[sample_column, drop = TRUE]])) {
  n_in <- sum(seu[[sample_column, drop = TRUE]] == s)
  n_out <- sum(seu_clean[[sample_column, drop = TRUE]] == s)
  write_summary(s, "\tbefore=", n_in, "\tafter=", n_out, "\n")
}
write_summary("\n")

if (!file.exists(clean_rds)) {
  # 仅保留 counts 层（丢弃轻量预处理产生的 data/scale.data/PCA/UMAP/聚类）
  seu_clean <- DietSeurat(seu_clean, assays = "RNA", layers = "counts",
                          dimreducs = NULL, graphs = NULL)
  saveRDS(seu_clean, file = clean_rds)
  cat("\nFile 2 : ", prefix, "seu_clean.rds saved successfully !\n", sep = "")
  cat("Note: output keeps raw counts only; normalization is done in the next step.\n")
} else {
  cat("\nseu_clean.rds already exists, skipping rds re-save (statistics above are from the cached object).\n")
}

close(summary_con)
cat("\nFile 3 : ", prefix, "Dedoublets_summary.txt saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")
