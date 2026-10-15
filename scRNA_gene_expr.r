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
    A script for visualizing target genes across cell types in single-cell data.
    Input is the annotated Seurat object from scRNA_annotate_step5
    (seu_annotated.rds, with UMAP + celltype). The script:
      1. Verifies input genes exist in the object (warns + skips missing)
      2. Computes average expression per cell type (mean of data layer)
      3. Exports a CSV (genes x celltypes) of average expression
      4. Renders 5 plots: Ridge / Violin / Feature / Dot / Heatmap
    All plots use Seurat built-in functions. UMAP-based plots require
    the input rds to have a UMAP reduction.
    Heatmap: cells are downsampled per celltype (-DS, default 100) and the
    top bar shows celltype colour blocks with names at 45 degrees.
  Example: scRNA_gene_expr -I 05_annotation/05.seu_annotated.rds -G CSF3R,S100A9,S100A8 -DS 100 -P 06. -O 06_geneExpr
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
  help = 'Path to an rds file for an annotated Seurat object (output of scRNA_annotate_step5). (Required)'
)
parser$add_argument(
  "-G",
  "--genes",
  dest = "genes",
  type = "character",
  required = TRUE,
  help = 'Comma-separated gene list, OR a path to a TSV/TXT file with one gene per line.'
)
# Optional
parser$add_argument(
  '-CT',
  '--celltype_col',
  dest = 'celltype_col',
  help = 'Metadata column name used for grouping (cell type). (Optional, default: celltype)',
  type = "character",
  default = "celltype"
)
parser$add_argument(
  '-CTO',
  '--celltype_order',
  dest = 'celltype_order',
  help = 'Order of cell types in the plots (comma- or space-separated, e.g. T,B,Macro). Optional; defaults to alphabetical order of the column.',
  type = "character"
)
parser$add_argument(
  '-RI',
  '--ridge_cols',
  dest = 'ridge_cols',
  help = 'Number of columns for RidgePlot facets (=0 uses the default 4 columns). (Optional, default: 0)',
  type = "integer",
  default = 0
)
parser$add_argument(
  '-DS',
  '--downsample',
  dest = 'downsample',
  type = "integer",
  default = 100,
  help = 'Downsample this many cells per celltype before drawing the heatmap (default: 100).'
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
  help = 'height of plots in inches (Optional, default: 10)',
  type = "numeric",
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
celltype_col <- opts$celltype_col
celltype_order_str <- opts$celltype_order
ridge_cols <- opts$ridge_cols
downsample_n <- opts$downsample
plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath

if (!(celltype_col %in% c("celltype"))) {
  # still accepted via -CT but keep flexibility
}

message("\n===============scRNA Gene Expression Visualization================\n")
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

# Summary 文件
summary_file <- file.path(outpath, paste0(prefix, "GeneExpression_summary.txt"))
summary_con <- file(summary_file, open = "wt")
write_summary <- function(...) {
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
write_summary("==================scRNA Gene Expression Summary===================\n")
write_summary("Current time: ", as.character(Sys.time()), "\n")
write_summary("Input RDS   : ", input_file, "\n")
write_summary("Celltype col: ", celltype_col, " | ridge_cols: ", ridge_cols, "\n")
write_summary("Downsample  : ", downsample_n, " cells per celltype (heatmap)\n")

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
})

# 3. Load data
seu <- readRDS(input_file)
cat("\nLoaded Seurat object from:", input_file, "\n")
cat("Cells:", ncol(seu), " | Genes:", nrow(seu), " | Layers:",
    paste(Layers(seu), collapse = ", "), "\n")
if (!(celltype_col %in% colnames(seu@meta.data))) {
  stop("Error: celltype column '", celltype_col, "' not found in metadata. Available: ",
       paste(colnames(seu@meta.data), collapse = ", "))
}
cat("Celltype groups (", celltype_col, "):\n", sep = "")
print(table(seu[[celltype_col, drop = TRUE]]))

# 4. Parse gene list
genes_arg <- opts$genes
if (file.exists(genes_arg)) {
  cat("\nReading gene list from file:", genes_arg, "\n")
  gene_lines <- readLines(genes_arg)
  genes <- unique(trimws(unlist(strsplit(paste(gene_lines, collapse = ","), "[, \t\n]+"))))
  genes <- genes[genes != ""]
} else {
  genes <- unique(trimws(unlist(strsplit(genes_arg, "[, \t]+"))))
  genes <- genes[genes != ""]
}
cat("\nInput gene count:", length(genes), "\n")

# Verify genes in object
genes_present <- intersect(genes, rownames(seu))
genes_missing <- setdiff(genes, rownames(seu))
cat("Genes found in Seurat object:", length(genes_present), "\n")
cat("Genes missing from Seurat object:", length(genes_missing), "\n")
if (length(genes_missing) > 0) {
  cat("Missing genes:\n")
  print(genes_missing)
}
if (length(genes_present) == 0) {
  stop("Error: none of the input genes are in the Seurat object. Check gene names/symbols.")
}

write_summary("Input genes : ", length(genes), " total, ",
              length(genes_present), " present in Seurat object\n")
write_summary("Missing     : ", if (length(genes_missing) > 0) paste(genes_missing, collapse = ", ") else "none", "\n")

# Set identities to cell type (used by VlnPlot / RidgePlot / DotPlot)
Idents(seu) <- celltype_col

# Optional: enforce cell type order
if (!is.null(celltype_order_str)) {
  ct_order <- unlist(strsplit(celltype_order_str, "[, ]+"))
  ct_order <- ct_order[ct_order != ""]
  ct_present <- intersect(ct_order, unique(seu[[celltype_col, drop = TRUE]]))
  ct_others <- setdiff(unique(seu[[celltype_col, drop = TRUE]]), ct_present)
  final_order <- c(ct_present, sort(ct_others))
  seu[[celltype_col, drop = TRUE]] <- factor(seu[[celltype_col, drop = TRUE]], levels = final_order)
  Idents(seu) <- seu[[celltype_col, drop = TRUE]]
  cat("\nCelltype order enforced:\n")
  print(levels(seu))
}

# 5. Compute average expression per cell type (mean of data layer)
cat("\n>>> Computing average expression per cell type (data layer, log-normalized) ...\n")
avg_exp <- AverageExpression(seu, features = genes_present,
                             group.by = celltype_col, layer = "data",
                             verbose = FALSE)$RNA
# avg_exp: genes (rows) x celltypes (cols), dense matrix
write_summary("\n---------- Average expression shape ----------\n")
write_summary("Genes x celltypes: ", nrow(avg_exp), " x ", ncol(avg_exp), "\n")

# Export as CSV
avg_csv <- file.path(outpath, paste0(prefix, "avg_expression.csv"))
df_avg <- as.data.frame(avg_exp)
df_avg <- cbind(gene = rownames(df_avg), df_avg)
write.csv(df_avg, file = avg_csv, row.names = FALSE, quote = FALSE)
cat("File 1 : ", prefix, "avg_expression.csv saved successfully !\n", sep = "")
write_summary("CSV saved    : ", avg_csv, "\n")

# 6. Ridge plot
cat("\n>>> Ridge plot ...\n")
# RidgePlot 用 Seurat 默认参数，仅控制列数。
# 注意：RidgePlot 的列数参数是 ncol（cols 是颜色向量，误传会报
# "Insufficient values in manual scale"）。
n_ct <- length(unique(seu[[celltype_col, drop = TRUE]]))
if (length(genes_present) == 1) {
  # 单基因：单面板（每 celltype 一行山脊）
  p_ridge <- RidgePlot(seu, features = genes_present, group.by = celltype_col) +
    theme(legend.position = "none")
  w_ridge <- plot_width
  h_ridge <- max(plot_height, n_ct * 0.3 + 1.5)
} else {
  # 多基因：ncol 控制列数（默认 4）；图例不需要（x 轴已是细胞类型）
  n_cols_r <- if (ridge_cols > 0) ridge_cols else 4L
  p_ridge <- RidgePlot(seu, features = genes_present, group.by = celltype_col,
                       ncol = n_cols_r) +
    theme(legend.position = "none")
  w_ridge <- plot_width
  n_rows_r <- ceiling(length(genes_present) / n_cols_r)
  # 每行高度 = celltype 数 × 行距系数（不能过大，保持整体宽高比协调）
  h_ridge <- max(plot_height, n_rows_r * (n_ct * 0.22 + 1.2))
}
# 上限控制：高度不超过 -W * 2（避免 30+ 基因变成超长纸带）
h_ridge <- min(h_ridge, plot_width * 2)
cat("RidgePlot size:", w_ridge, "x", h_ridge, "\n")
ggsave(file.path(outpath, paste0(prefix, "ridgeplot.png")), p_ridge,
       width = w_ridge, height = h_ridge, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(outpath, paste0(prefix, "ridgeplot.pdf")), p_ridge,
       width = w_ridge, height = h_ridge, limitsize = FALSE)
cat("Plot 1 : ", prefix, "ridgeplot.png(pdf) has been successfully generated !\n", sep = "")

# 7. Violin plot
cat("\n>>> Violin plot ...\n")
# 按基因数自动选列数（取整以让整体宽高比接近 -W/-H）；默认每列 = 一个基因
n_v <- length(genes_present)
target_ratio <- plot_width / plot_height
best_ncol_v <- 2L
best_err_v <- Inf
for (nc in 1:min(8L, n_v)) {
  nr <- ceiling(n_v / nc)
  r <- nc / nr
  e <- abs(r - target_ratio)
  if (e < best_err_v) { best_err_v <- e; best_ncol_v <- nc }
}
n_cols_v <- best_ncol_v
n_rows_v <- ceiling(n_v / n_cols_v)
unit_w_v <- plot_width / n_cols_v
h_v <- n_rows_v * unit_w_v   # 1:1 子图布局
w_v <- plot_width
# 上限：避免基因过多（>30）时高度超 50 英寸（ggsave 硬限制）
h_v <- min(h_v, plot_width * 2)
cat("ViolinPlot layout:", n_cols_v, "cols x", n_rows_v, "rows  (", n_v, "genes, ",
    sprintf("%.1f x %.1f inches", w_v, h_v), " )\n")
p_v <- VlnPlot(seu, features = genes_present, group.by = celltype_col,
               pt.size = 0, ncol = n_cols_v) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")
ggsave(file.path(outpath, paste0(prefix, "violinplot.png")), p_v,
       width = w_v, height = h_v, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(outpath, paste0(prefix, "violinplot.pdf")), p_v,
       width = w_v, height = h_v, limitsize = FALSE)
cat("Plot 2 : ", prefix, "violinplot.png(pdf) has been successfully generated !\n", sep = "")

# 8. Feature plot (UMAP)
cat("\n>>> Feature plot ...\n")
if ("umap" %in% names(seu@reductions)) {
  # 每个子图 aspect.ratio=1（正方形）；自动选列数使整体宽高比最接近 -W/-H 比例
  # 注意：列数上限不得超过基因数，否则（单基因时）会按多列画布出图、右侧留白
  n_f <- length(genes_present)
  target_ratio <- plot_width / plot_height
  best_ncol <- 1L
  best_err <- Inf
  for (nc in seq_len(min(n_f, 5))) {
    nr <- ceiling(n_f / nc)
    r <- nc / (nr * 1.1)   # 1.1：行间留白
    e <- abs(r - target_ratio)
    if (e < best_err) { best_err <- e; best_ncol <- nc }
  }
  n_cols_f <- best_ncol
  n_rows_f <- ceiling(n_f / n_cols_f)
  unit_w <- plot_width / n_cols_f
  h_f <- n_rows_f * unit_w * 1.1
  w_f <- plot_width
  cat("FeaturePlot layout:", n_cols_f, "cols x", n_rows_f, "rows\n")
  p_f <- FeaturePlot(seu, features = genes_present, reduction = "umap",
                     cols = c("lightgrey", "red"), ncol = n_cols_f) &
    theme(aspect.ratio = 1)
  ggsave(file.path(outpath, paste0(prefix, "featureplot.png")), p_f,
         width = w_f, height = h_f, dpi = 300, bg = "white", limitsize = FALSE)
  ggsave(file.path(outpath, paste0(prefix, "featureplot.pdf")), p_f,
         width = w_f, height = h_f, limitsize = FALSE)
  cat("Plot 3 : ", prefix, "featureplot.png(pdf) has been successfully generated !\n", sep = "")
} else {
  cat("[Info] No 'umap' reduction in the input object - skipping FeaturePlot.\n")
}

# 9. Dot plot
cat("\n>>> Dot plot ...\n")
# DotPlot 横轴是基因（多→长），纵轴是 celltype（少→短）；
# 整体宽度按基因数线性扩展，高度按 celltype 数 × 行高算
n_d <- length(genes_present)
n_ct_d <- length(unique(seu[[celltype_col, drop = TRUE]]))
# 宽度：每个基因占 ~0.55 英寸（RotatedAxis 后垂直标签约 0.5-0.6），最低不小于 -W
w_d <- max(plot_width, n_d * 0.55 + 2)
# 高度：每个 celltype 约 0.32 英寸，最低不小于 -H
h_d <- max(plot_height, n_ct_d * 0.32 + 1.5)
cat("DotPlot size:", w_d, "x", h_d, "(", n_d, "genes,", n_ct_d, "celltypes )\n")
p_d <- DotPlot(seu, features = genes_present, group.by = celltype_col,
               cols = c("lightgrey", "red"), dot.scale = 6) +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "right")
ggsave(file.path(outpath, paste0(prefix, "dotplot.png")), p_d,
       width = w_d, height = h_d, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(outpath, paste0(prefix, "dotplot.pdf")), p_d,
       width = w_d, height = h_d, limitsize = FALSE)
cat("Plot 4 : ", prefix, "dotplot.png(pdf) has been successfully generated !\n", sep = "")

# 10. Heatmap (Seurat DoHeatmap)
cat("\n>>> Heatmap ...\n")
# DoHeatmap 默认用 scale.data 显示；scale.data 通常只含 HVG，
# 非 HVG 的目标基因会被省略 → 先对目标基因单独 ScaleData，再画图。
# group.bar=FALSE 移除顶部 celltype 色条。
# celltype 顺序：沿用对象默认顺序（或前面 -CTO 指定的顺序）；
# 因为下面已按每组 downsample=100 采样，细胞数已拉平，无需再按计数排序。

# 对目标基因做 scale（features 参数确保只 scale 用户基因）
cat(">>> ScaleData on target genes ...\n")
seu <- ScaleData(seu, features = genes_present, verbose = FALSE)

# 按 celltype 下采样（每组最多 -DS 个细胞，默认 100），降低长时渲染
# 注：Seurat 5.5.0 的 subset(downsample=) 经实测是【确定性】的——两个独立
# Rscript 进程抽到的细胞完全相同（1800/1800），因此该步骤天然可复现。
# 仍固定随机种子作为防御：若未来 Seurat 版本改为真随机抽样，可保证结果不漂移。
set.seed(42)
seu_ds <- subset(seu, downsample = downsample_n)
cat("Downsampled cells:", ncol(seu_ds), "across",
    length(unique(Idents(seu_ds))), "celltypes (max", downsample_n, "per celltype)\n")
cat("Heatmap celltype order:\n")
print(levels(Idents(seu_ds)))

# 顶部保留 celltype 颜色块 + 名字（Seurat 默认 angle = 45 斜排）
# Heatmap 宽高随基因数 + celltype 数自适应：
#   宽度：每个基因占 ~0.45 英寸（45° 标签宽度），最低不小于 -W
#   高度：每个 celltype 约 0.45 英寸，最低不小于 -H
# 同时加上下限避免极端情况（图过长或过矮）
n_g_h <- length(genes_present)
n_ct_h <- length(levels(Idents(seu_ds)))
w_h <- max(plot_width, n_g_h * 0.45 + 1.5)
h_h <- max(plot_height, n_ct_h * 0.45 + 2)
# 上限：宽度不超过 -W*1.5，高度不超过 -W*1.5（避免极端长图）
w_h <- min(w_h, plot_width * 1.5)
h_h <- min(h_h, plot_width * 1.5)
cat("Heatmap size:", w_h, "x", h_h, "(", n_g_h, "genes,", n_ct_h, "celltypes )\n")
ph <- DoHeatmap(seu_ds, features = genes_present, group.by = celltype_col,
                group.bar = TRUE,             # 顶部 celltype 颜色块
                label = TRUE,                 # 显示细胞类型名字
                angle = 45, size = 5,
                draw.lines = TRUE, raster = FALSE) +
  theme(axis.text.y = element_text(size = 8)) +
  labs(title = NULL)
ggsave(file.path(outpath, paste0(prefix, "heatmap.png")), ph,
       width = w_h, height = h_h, dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(outpath, paste0(prefix, "heatmap.pdf")), ph,
       width = w_h, height = h_h, limitsize = FALSE)
cat("Plot 5 : ", prefix, "heatmap.png(pdf) has been successfully generated !\n", sep = "")

close(summary_con)
cat("\nFile 2 : ", prefix, "GeneExpression_summary.txt saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")