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
    A script for single-cell data clustering analysis (Seurat v5).
    Input is the harmony-corrected object from scRNA_Nor_harmony_step3.
    Workflow: FindNeighbors -> multi-resolution clustree exploration ->
    (user picks -Res) single-resolution FindClusters + UMAP -> FindAllMarkers DE.
    Two output products are cached independently:
      1) clustree.png        : exploration result (resolution sweep)
      2) seu_res<Res>_clustered.rds : final clustered object (incl. UMAP)
    If the clustered rds already exists and is complete, the script skips all
    computation and only re-plots from the cached object (use for -W/-H changes).
  Example: scRNA_clustering_step4 -I 03_Nor_Harmony/03.seu_harmony.rds -Red harmony -PC 30 -Res 0.7 -T 20 -P 04. -O 04_clustering
  Example (exploration only, no -Res): scRNA_clustering_step4 -I 03_Nor_Harmony/03.seu_harmony.rds -Red harmony -PC 30 -T 20 -P 04. -O 04_clustering
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
  help = 'Path to an rds file for a seurat object used for scRNA data clustering. (Required)'
)

# Optional
parser$add_argument(
  '-PC',
  dest = 'PC_num',
  help = 'Number of PCs used for dimension reduction and clustering. (Optional, default: 30)',
  type = "integer",
  default = 30
)
parser$add_argument(
  '-Red',
  dest = 'reduction',
  help = 'Reduction method for clustering. One of "pca" or "harmony". (Optional, default: harmony)',
  type = "character",
  default = "harmony"
)
parser$add_argument(
  '-Res',
  dest = 'resolution',
  help = 'Resolution value for the final clustering. Omit to run clustree exploration only.',
  type = "double"
)
parser$add_argument(
  '-CA',
  dest = 'cluster_algorithm',
  help = 'Algorithm used in "FindClusters" function, one of "Louvain" or "Leiden". If specified as "Leiden", you must configure the Python environment and use the -P_env parameter to specify the path to that Python environment. (Optional, default: Louvain)',
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
  help = 'The path to the Python environment required to use the Leiden algorithm. When the -CA parameter is "Leiden", it must be specified.',
  type = "character"
)
parser$add_argument(
  '-skipDE',
  action = "store_true",
  default = FALSE,
  help = "Skip the FindAllMarkers differential expression analysis and per-cluster dotplots. (Optional; default: FALSE)"
)
parser$add_argument(
  '-T',
  dest = 'threads',
  help = 'Number of threads used for clustering analysis. (Optional, default: 10)',
  type = "integer",
  default = 10
)
parser$add_argument(
  '-W',
  '--width',
  dest = 'plot_width',
  help = 'width of UMAP plots in inches (Optional, default: 10)',
  type = "numeric",
  default = 10
)
parser$add_argument(
  '-H',
  '--height',
  dest = 'plot_height',
  help = 'height of UMAP plots in inches (Optional, default: 8)',
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
PC_num <- opts$PC_num
reduction <- opts$reduction
resolution <- opts$resolution
cluster_algorithm <- opts$cluster_algorithm
assay <- opts$assay
python_env <- opts$python_env
skipDE <- opts$skipDE
threads <- opts$threads
plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath

if (!is.null(python_env) && !file.exists(python_env)) {
  stop("The specified Python environment path does not exist.")
}
if (!(reduction %in% c("pca", "harmony"))) {
  stop("Invalid reduction. Choose either 'pca' or 'harmony'.")
}

message("=========================scRNA Clustering=========================\n")
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
summary_file <- file.path(outpath, paste0(prefix, "Clustering_summary.txt"))
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
write_summary("=====================scRNA Clustering Summary=====================\n")
write_summary("Current time: ", as.character(Sys.time()), "\n")
write_summary("Input RDS : ", input_file, "\n")
write_summary("Reduction : ", reduction, " | PC: ", PC_num,
              " | resolution: ", if (is.null(resolution)) "NULL (exploration)" else resolution,
              " | algorithm: ", cluster_algorithm, " | skipDE: ", skipDE, "\n")

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(clustree,quietly = TRUE)){install.packages("clustree")}
  if(!require(reticulate,quietly = TRUE)){install.packages("reticulate")}
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
  } else {
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

# ============================================================================
# 3. 产物路径定义 + rds 完整性检查
# ============================================================================
clustree_png <- file.path(outpath, paste0(prefix, "clustree.png"))
res_str <- if (is.null(resolution)) "" else paste0("res", resolution)
clustered_rds <- file.path(outpath, paste0(prefix, "seu_", res_str, "_clustered.rds"))

is_complete_rds <- function(path) {
  # 快速预检：文件存在且大小合理（rds 完整对象通常 > 1 MB，半成品多为 0 字节/极小）
  if (!file.exists(path)) return(FALSE)
  if (file.info(path)$size < 1e6) {
    cat("[Warning] Found", path, "but it is too small (", file.info(path)$size,
        " bytes) - treating as incomplete, will re-run.\n")
    return(FALSE)
  }
  # 真实读取验证（消耗时间，但确保不是半成品）
  ok <- tryCatch({
    obj <- readRDS(path)
    !is.null(obj)
  }, error = function(e) {
    cat("[Warning] Found", path, "but readRDS failed - file is incomplete/corrupted:",
        conditionMessage(e), "\n")
    cat("          Will re-run. If this persists, delete the file manually.\n")
    FALSE
  })
  ok
}

if (!is.null(resolution) && file.exists(clustered_rds)) {
  # ---- 分支 A：最终 rds 已存在且完整 → 只重画图，跳过一切计算 ----
  cat("\n>>> Found existing clustered rds: ", clustered_rds, "\n", sep = "")
  if (is_complete_rds(clustered_rds)) {
    cat("    File is complete. Skipping FindNeighbors/clustree/FindClusters/UMAP/DE ...\n")
    cat("    Loading cached object and re-plotting only.\n")
    seu <- readRDS(clustered_rds)
    cat("    Loaded cells:", ncol(seu), " | clusters:",
        length(unique(seu$seurat_clusters)), "\n")
    if (!("umap" %in% names(seu@reductions))) {
      stop("Cached rds lacks the 'umap' reduction - looks stale. Delete the rds to re-run.")
    }
    Idents(seu) <- "seurat_clusters"
    seu$seurat_clusters <- factor(seu$seurat_clusters,
                                  levels = sort(unique(as.numeric(as.character(seu$seurat_clusters)))))
    cached_res <- resolution
    # 画图（同下游公共画图函数）
    source_plot <- "cached"
  } else {
    # rds 存在但不完整 → 走完整计算（会覆盖半成品）
    cat("\n>>> Clustered rds is incomplete; running the full clustering pipeline.\n")
    seu <- readRDS(input_file)
    cat("\nLoaded Seurat object from:", input_file, "\n")
    cat("Cells:", ncol(seu), " | Genes:", nrow(seu), "\n")
    if (!(reduction %in% names(seu@reductions))) {
      stop("Error: reduction '", reduction, "' not found in the input object. Available: ",
           paste(names(seu@reductions), collapse = ", "))
    }
    # ---- 完整流程 ----
    cat("\n>>> FindNeighbors (reduction = ", reduction, ", dims = 1:", PC_num, ") ...\n", sep = "")
    seu <- FindNeighbors(seu, reduction = reduction, dims = 1:PC_num)
    # 若 clustree.png 不存在 → 跑多分辨率探索
    if (!file.exists(clustree_png)) {
      if (algorithm_num == 1) {
        resolutions <- seq(0, 1, by = 0.1)
      } else {
        resolutions <- seq(0.1, 1, by = 0.1)
      }
      cat(">>> clustree exploration across resolutions:", paste(resolutions, collapse = ", "), "\n")
      seu <- FindClusters(seu, resolution = resolutions, algorithm = algorithm_num,
                          random.seed = 123, verbose = FALSE)
      p1 <- clustree(seu, prefix = paste0(assay, "_snn_res.")) +
        labs(title = "Cluster trees for different resolution clustering") +
        coord_flip()
      ggsave(file.path(outpath, paste0(prefix, "clustree.png")), p1,
             width = 24, height = 16, dpi = 300, bg = "white")
      ggsave(file.path(outpath, paste0(prefix, "clustree.pdf")), p1,
             width = 24, height = 16)
      cat("\nPlot 1 : ", prefix, "clustree.png(pdf) has been successfully generated !\n", sep = "")
    } else {
      cat("\n>>> clustree.png already exists, skipping multi-resolution exploration.\n")
      cat("    (delete ", clustree_png, " to re-run exploration)\n", sep = "")
    }
    cat("\n>>> FindClusters (resolution = ", resolution, ") ...\n", sep = "")
    seu <- FindClusters(seu, resolution = resolution, algorithm = algorithm_num,
                        random.seed = 123, verbose = FALSE)
    cat(">>> RunUMAP (reduction = ", reduction, ") ...\n", sep = "")
    seu <- RunUMAP(seu, reduction = reduction, dims = 1:PC_num, verbose = FALSE)
    cat("UMAP completed.\n")
    seu$seurat_clusters <- factor(seu$seurat_clusters,
                                  levels = sort(unique(as.numeric(as.character(seu$seurat_clusters)))))
    Idents(seu) <- "seurat_clusters"
    # 保存 rds（完整）
    saveRDS(seu, file = clustered_rds)
    cat("\nFile 2 : ", prefix, "seu_", res_str, "_clustered.rds saved successfully !\n", sep = "")
    source_plot <- "computed"
  }
} else {
  # ---- 分支 B：无 -Res 或 rds 不存在 → 需要计算 ----
  seu <- readRDS(input_file)
  cat("\nLoaded Seurat object from:", input_file, "\n")
  cat("Cells:", ncol(seu), " | Genes:", nrow(seu), "\n")
  if (!(reduction %in% names(seu@reductions))) {
    stop("Error: reduction '", reduction, "' not found in the input object. Available: ",
         paste(names(seu@reductions), collapse = ", "))
  }
  cat("\n>>> FindNeighbors (reduction = ", reduction, ", dims = 1:", PC_num, ") ...\n", sep = "")
  seu <- FindNeighbors(seu, reduction = reduction, dims = 1:PC_num)

  if (!file.exists(clustree_png)) {
    if (algorithm_num == 1) {
      resolutions <- seq(0, 1, by = 0.1)
    } else {
      resolutions <- seq(0.1, 1, by = 0.1)
    }
    cat(">>> clustree exploration across resolutions:", paste(resolutions, collapse = ", "), "\n")
    seu <- FindClusters(seu, resolution = resolutions, algorithm = algorithm_num,
                        random.seed = 123, verbose = FALSE)
    p1 <- clustree(seu, prefix = paste0(assay, "_snn_res.")) +
      labs(title = "Cluster trees for different resolution clustering") +
      coord_flip()
    ggsave(file.path(outpath, paste0(prefix, "clustree.png")), p1,
           width = 24, height = 16, dpi = 300, bg = "white")
    ggsave(file.path(outpath, paste0(prefix, "clustree.pdf")), p1,
           width = 24, height = 16)
    cat("\nPlot 1 : ", prefix, "clustree.png(pdf) has been successfully generated !\n", sep = "")
  } else {
    cat("\n>>> clustree.png already exists, skipping multi-resolution exploration.\n")
  }

  if (is.null(resolution)) {
    cat("\n[Info] No -Res provided: exploration mode only.\n")
    cat("       Inspect ", prefix, "clustree.png(pdf), choose a resolution,\n", sep = "")
    cat("       then re-run with -Res <value> to produce the final clustered object.\n")
    write_summary("\n[Info] Exploration mode only (no -Res): clustree produced. Re-run with -Res to cluster.\n")
    close(summary_con)
    cat("\n============================All Done============================\n")
    quit(status = 0)
  }

  cat("\n>>> FindClusters (resolution = ", resolution, ") ...\n", sep = "")
  seu <- FindClusters(seu, resolution = resolution, algorithm = algorithm_num,
                      random.seed = 123, verbose = FALSE)
  cat(">>> RunUMAP (reduction = ", reduction, ") ...\n", sep = "")
  seu <- RunUMAP(seu, reduction = reduction, dims = 1:PC_num, verbose = FALSE)
  cat("UMAP completed.\n")
  seu$seurat_clusters <- factor(seu$seurat_clusters,
                                levels = sort(unique(as.numeric(as.character(seu$seurat_clusters)))))
  Idents(seu) <- "seurat_clusters"
  saveRDS(seu, file = clustered_rds)
  cat("\nFile 2 : ", prefix, "seu_", res_str, "_clustered.rds saved successfully !\n", sep = "")
  source_plot <- "computed"
}

# ============================================================================
# 4. 公共画图（两种来源共用）
# ============================================================================
res_lab <- if (is.null(resolution)) "" else paste0("res", resolution, "_")
p2 <- DimPlot(seu, reduction = "umap", label = TRUE, label.size = 6)
ggsave(file.path(outpath, paste0(prefix, res_lab, "cluster_umap.png")), p2,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, res_lab, "cluster_umap.pdf")), p2,
       width = plot_width, height = plot_height)
cat("\nPlot 2 : ", prefix, res_lab, "cluster_umap.png(pdf) has been successfully generated !\n", sep = "")

p3 <- DimPlot(seu, reduction = "umap", group.by = "sample", label = FALSE)
ggsave(file.path(outpath, paste0(prefix, res_lab, "sample_umap.png")), p3,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, res_lab, "sample_umap.pdf")), p3,
       width = plot_width, height = plot_height)
cat("Plot 3 : ", prefix, res_lab, "sample_umap.png(pdf) has been successfully generated !\n", sep = "")

p3b <- DimPlot(seu, reduction = "umap", group.by = "group", label = FALSE)
ggsave(file.path(outpath, paste0(prefix, res_lab, "group_umap.png")), p3b,
       width = plot_width, height = plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, res_lab, "group_umap.pdf")), p3b,
       width = plot_width, height = plot_height)
cat("Plot 4 : ", prefix, res_lab, "group_umap.png(pdf) has been successfully generated !\n", sep = "")

write_summary("\n---------- Cluster composition ----------\n")
cat("\nCells per cluster:\n")
print(table(seu$seurat_clusters))
tb_clu <- as.data.frame(table(seu$seurat_clusters), stringsAsFactors = FALSE)
colnames(tb_clu) <- c("cluster", "cells")
cat(capture.output(write.table(tb_clu, sep = "\t", row.names = FALSE, quote = FALSE)),
    file = summary_con, sep = "\n")
cat("\n", file = summary_con)
write_summary("\n---------- Cells per sample ----------\n")
tb_smp <- as.data.frame(table(seu$sample), stringsAsFactors = FALSE)
colnames(tb_smp) <- c("sample", "cells")
cat(capture.output(write.table(tb_smp, sep = "\t", row.names = FALSE, quote = FALSE)),
    file = summary_con, sep = "\n")
cat("\n", file = summary_con)

# ============================================================================
# 5. DE 分析（仅当计算来源 && 未 skip）
# ============================================================================
if (is.null(resolution)) {
  cat("\n[Info] Exploration mode - DE analysis skipped.\n")
} else if (skipDE) {
  cat("\n[Info] -skipDE set - skipping differential expression analysis.\n")
} else {
  # 若是从缓存加载，DE csv 已存在则跳过重算
  cluster_marker_dir <- file.path(outpath, paste0(prefix, "clusters_DE"))
  markers_csv <- file.path(cluster_marker_dir, paste0(prefix, res_str, "_cluster_markers.csv"))
  if (exists("source_plot") && source_plot == "cached" && file.exists(markers_csv)) {
    cat("\n[Info] Cached run and cluster markers csv already exists - skipping DE re-computation.\n")
  } else {
    cat("\nStart differential expression analysis for the clusters ...\n")
    if (!dir.exists(cluster_marker_dir)) {
      dir.create(cluster_marker_dir)
    }
    Idents(seu) <- "seurat_clusters"
    cluster_markers_df <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.1,
                                         logfc.threshold = 0.25, verbose = TRUE)
    write.csv(cluster_markers_df,
              file = markers_csv, row.names = FALSE, quote = FALSE)
    cat("\nFile 3 : ", prefix, res_str, "_cluster_markers.csv saved successfully !\n", sep = "")
    # Top 10 markers for each cluster
    cluster_top10_markers <- cluster_markers_df %>%
      filter(p_val_adj < 0.05, pct.1 > 0.3) %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = 10)
    top10_csv <- file.path(cluster_marker_dir, paste0(prefix, res_str, "_cluster_top10_markers.csv"))
    write.csv(cluster_top10_markers, file = top10_csv, row.names = FALSE, quote = FALSE)
    cat("File 4 : ", prefix, res_str, "_cluster_top10_markers.csv saved successfully !\n", sep = "")

    for (clu in sort(unique(seu$seurat_clusters))) {
      cat("\nGenerating dotplot for cluster:", clu, "...\n")
      cluster_markers <- cluster_top10_markers %>%
        filter(cluster == !!clu)
      cat("Number of markers for cluster", clu, ":", nrow(cluster_markers), "\n")
      if (nrow(cluster_markers) == 0) { next }
      p <- DotPlot(seu, features = cluster_markers$gene,
                   cols = c("white", "orangered2"), dot.scale = 8, cluster.idents = TRUE) +
        RotatedAxis() +
        scale_x_discrete("Markers") +
        scale_y_discrete("Seurat_Clusters") +
        ggtitle(paste("Cluster", clu))
      ggsave(file.path(cluster_marker_dir, paste0(prefix, res_str, "_cluster_", clu, "_dotplot.png")),
             p, width = 10, height = 10, dpi = 300, bg = "white")
      ggsave(file.path(cluster_marker_dir, paste0(prefix, res_str, "_cluster_", clu, "_dotplot.pdf")),
             p, width = 10, height = 10)
      rm(p)
    }
    cat("\nDE analysis completed !\n")
  }
}

close(summary_con)
cat("\nFile 5 : ", prefix, "Clustering_summary.txt saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")
