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
    A script for cell-type annotation of single-cell data (Seurat v5).
    Input is the clustered object from scRNA_clustering_step4
    (seu_res<Res>_clustered.rds, with UMAP + cluster + DE in clusters_DE/).

    Two modes:
      -MB <markers.tsv>   : exploration (cluster x marker DotPlot). Does NOT write rds.
                            TSV columns: V1=celltype, V2=markers (comma-separated).
      -MA <annotation.tsv> : final annotation (cluster -> celltype map). Writes
                            seu_annotated.rds and produces celltype UMAP + DotPlot.
                            TSV columns: Cluster / Markers / CellType.

    Exactly one of -MB/-MA must be supplied.
  Example: scRNA_annotate_step5 -I 04_clustering_louvain/04.seu_res0.1_clustered.rds -MB Markers.txt -P 05. -O 05_annotation
  Example: scRNA_annotate_step5 -I 04_clustering_louvain/04.seu_res0.1_clustered.rds -MA 05_annotation/annotation.txt -P 05. -O 05_annotation
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
  help = 'Path to an rds file for a seurat object used for scRNA data annotation. (Required)'
)
parser$add_argument(
  "-MB",
  "--markers_before",
  dest = "markers_before_annotation",
  type = "character",
  help = "Path to a tsv file used for annotation that contains 2 columns: the first column is the cell type, and the second column is the marker gene(s) (separated by commas)."
)
parser$add_argument(
  "-MA",
  "--markers_after",
  dest = "markers_after_annotation",
  type = "character",
  help = "Path to a tsv file used for annotation that contains 3 columns: (1) Cluster: cluster number, (2)Markers: marker gene(s) (separated by commas), and (3)CellType: annotated cell type."
)

# Optional
parser$add_argument(
  '-S',
  "--split",
  dest = 'split',
  action = 'store_true',
  help = 'Whether to split the markers by cell type according to the first column in the markers file provided by the -MB parameter. This should be used when the number of markers for each cell type is too large. (Optional; default: FALSE)'
)
parser$add_argument("-CO", "--celltype_order",
                    help = 'Order of cell types in the dot plot (comma- or space-separated, e.g. B_cell, T_cell, NK).'
)
parser$add_argument(
  '-W',
  '--width',
  dest = 'plot_width',
  help = 'width of plots in inches (Optional, default: 24)',
  type = "numeric",
  default = 24
)
parser$add_argument(
  '-H',
  '--height',
  dest = 'plot_height',
  help = 'height of plots in inches (Optional, default: 14)',
  type = "numeric",
  default = 14
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
markers_before <- opts$markers_before_annotation
markers_after <- opts$markers_after_annotation
split <- opts$split


plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath
if (!is.null(opts$celltype_order)) {
  celltype_order <- unlist(strsplit(opts$celltype_order, "[, ]+"))
} else {
  celltype_order <- NULL
}

if (!is.null(split) & is.null(markers_before) & is.null(markers_after)) {
  stop("If you want to split the markers by cell type, you must provide a file with cell types and marker genes using the -MB parameter.")
}
if (is.null(markers_before) == is.null(markers_after)) {
  stop("Either the -MB parameter or the -MA parameter must be specified, but not both.")
}

message("=========================scRNA Annotation=========================\n")
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
summary_file <- file.path(outpath, paste0(prefix, "Annotation_summary.txt"))
summary_con <- file(summary_file, open = "wt")
write_summary <- function(...) {
  # 入参任一为 length>1 向量时会被 paste0 collapse 压行 → 报错
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
write_summary("=====================scRNA Annotation Summary=====================\n")
write_summary("Current time: ", as.character(Sys.time()), "\n")
write_summary("Input RDS   : ", input_file, "\n")
write_summary("Mode        : ", if (!is.null(markers_before)) "-MB (exploration)" else "-MA (final annotation)",
              " | split: ", split, "\n")

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
})

# 3. 加载数据 + rds 完整性检查（与 step4 同款逻辑）
annotated_rds <- file.path(outpath, paste0(prefix, "seu_annotated.rds"))

is_complete_rds <- function(path) {
  if (!file.exists(path)) return(FALSE)
  if (file.info(path)$size < 1e6) {
    cat("[Warning] Found", path, "but it is too small (", file.info(path)$size,
        " bytes) - treating as incomplete, will re-run.\n")
    return(FALSE)
  }
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

use_cached <- FALSE
if (file.exists(annotated_rds) && is_complete_rds(annotated_rds)) {
  cat("\n>>>Found existing seu_annotated.rds (complete). Loading cached object.\n")
  cat("    Path:", annotated_rds, "\n")
  cached_seu <- readRDS(annotated_rds)
  cat("    Loaded cells:", ncol(cached_seu), " | celltype levels:",
      paste(levels(cached_seu$celltype), collapse = ", "), "\n")
  use_cached <- TRUE
} else {
  if (file.exists(annotated_rds)) {
    cat("\n>>>seu_annotated.rds exists but is incomplete; will re-run and overwrite.\n")
  }
  cached_seu <- NULL
}

# 即使缓存可用，-MA 也要重新读 TSV（选项 B：rds 完整时仍读 -MA，便于改映射后重写）
seu <- readRDS(input_file)
genes_in_seu <- rownames(seu)
cat("\nLoaded clustered Seurat object from:", input_file, "\n")
cat("Cells:", ncol(seu), " | Genes:", nrow(seu), "\n")

# 4. Before annotate (exploration mode, -MB)
if (!is.null(markers_before) && !split) {
  if (use_cached) cat("\n[Info] Cached seu_annotated.rds loaded, but -MB was given. Re-painting exploration DotPlot only.\n")
  markers_before_df <- read.table(markers_before, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  markers_before_genes <- unique(unlist(strsplit(markers_before_df$V2, ',')))
  existing_genes <- intersect(markers_before_genes, genes_in_seu)
  missing_genes <- setdiff(markers_before_genes, genes_in_seu)
  cat("\nTotal input genes:", length(markers_before_genes), "\n")
  cat("Genes found in Seurat object:", length(existing_genes), "\n")
  cat("Genes missing from Seurat object:", length(missing_genes), "\n")
  if (length(missing_genes) > 0) {
    cat("Missing genes:\n")
    print(missing_genes)
  }
  dedup_markers_df <- unique(markers_before_df[c("V1", "V2")])
  dedup_df_sorted <- dedup_markers_df[order(dedup_markers_df$V1), ]
  if (!is.null(celltype_order)) {
    dedup_df_sorted <- dedup_markers_df %>%
      mutate(V1 = factor(V1, levels = celltype_order)) %>%
      filter(!is.na(V1)) %>%
      arrange(V1)
  }
  markers_list <- setNames(strsplit(dedup_df_sorted$V2, ","), dedup_df_sorted$V1)
  seen <- character(0)
  markers_list <- lapply(markers_list, function(x) {
    new_genes <- x[!x %in% seen]
    seen <<- c(seen, new_genes)
    return(new_genes)
  })
  empty_groups <- names(markers_list)[sapply(markers_list, length) == 0]
  if (length(empty_groups) > 0) {
    message("The following cell types have no marker genes left after removing duplicates:")
    print(empty_groups)
  } else {
    message("All cell types retained at least one marker gene after deduplication.")
  }
  Idents(seu) <- "seurat_clusters"
  p1 <- DotPlot(
    seu,
    features = markers_list,
    cols = c("grey60", "red"),
    cluster.idents = TRUE,
    assay = "RNA",
    dot.scale = 8
  ) +
    RotatedAxis() +
    theme(
      strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
      panel.border = element_rect(color = "black", size = 1),
      panel.spacing = unit(2, "mm"),
      strip.text = element_text(margin = margin(b = 2, t = 2, unit = "mm"), face = "bold", size = 12),
      strip.placement = 'outside',
      axis.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 14)
    ) + scale_x_discrete("Markers") + scale_y_discrete("Seurat_Clusters")
  ggsave(file.path(outpath, paste0(prefix, "cluster_markers.png")), p1,
         width = plot_width, height = plot_height, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "cluster_markers.pdf")), p1,
         width = plot_width, height = plot_height)
  cat("\nPlot : ", prefix, "cluster_markers.png(pdf) has been generated successfully !\n", sep = "")
}

if (!is.null(markers_before) && split) {
  if (use_cached) cat("\n[Info] Cached seu_annotated.rds loaded, but -MB + -S was given. Re-painting per celltype DotPlots only.\n")
  dotplot_dir <- file.path(outpath, paste0(prefix, "markers_dotplot"))
  if (!dir.exists(dotplot_dir)) {
    dir.create(dotplot_dir, recursive = TRUE)
  }
  markers_before_df <- readLines(markers_before)
  markers_list <- list()
  for (line in markers_before_df) {
    split_line <- strsplit(line, "\t")[[1]]
    subset_name <- split_line[1]
    subset_content <- strsplit(split_line[2], ",")[[1]]
    markers_list[[subset_name]] <- unique(subset_content)
  }
  Idents(seu) <- "seurat_clusters"
  for (ct in names(markers_list)) {
    cat("\nCelltype processing: ", ct, " ......\n")
    markers <- markers_list[[ct]]
    existing_genes <- intersect(markers, genes_in_seu)
    missing_genes <- setdiff(markers, genes_in_seu)
    cat("Total input genes:", length(markers), "\n")
    cat("Genes found in Seurat object:", length(existing_genes), "\n")
    cat("Genes missing from Seurat object:", length(missing_genes), "\n")
    if (length(missing_genes) > 0) {
      cat("Missing genes:\n")
      print(missing_genes)
    }
    p <- DotPlot(seu, features = markers, cols = c("grey60", "red"), dot.scale = 8, cluster.idents = TRUE) +
      RotatedAxis() + scale_x_discrete("Markers") + scale_y_discrete("Seurat_Clusters") +
      labs(title = paste0("Markers for ", ct))
    ct_safe <- gsub(" ", "_", ct)
    ggsave(file.path(dotplot_dir, paste0(prefix, ct_safe, "_cluster_markers.png")), p,
           width = plot_width, height = plot_height, bg = "white", dpi = 300)
    ggsave(file.path(dotplot_dir, paste0(prefix, ct_safe, "_cluster_markers.pdf")), p,
           width = plot_width, height = plot_height)
    rm(list = c("markers", "p"))
  }
  cat("Plots saved in directory", paste0(prefix, "markers_dotplot"), "successfully !\n")
}

# 5. Annotate (final mode, -MA)
if (!is.null(markers_after)) {
  cat("\n>>> Final annotation mode (-MA) ...\n")
  if (use_cached) cat("    [Info] Cached seu_annotated.rds exists but -MA was given. Re-mapping cluster -> celltype from the TSV and overwriting rds.\n")

  celltype_df <- read.delim(markers_after, sep = "\t", stringsAsFactors = FALSE)
  cluster_to_celltype <- celltype_df %>%
    select(Cluster, CellType) %>%
    deframe()

  seu$celltype <- unname(cluster_to_celltype[as.character(seu$seurat_clusters)])
  seu$celltype <- as.character(seu$celltype)

  # 报告未映射 cluster
  unmapped <- setdiff(unique(as.character(seu$seurat_clusters)), names(cluster_to_celltype))
  if (length(unmapped) > 0) {
    cat("[Warning] Some clusters are NOT in -MA file and were set to NA:", paste(unmapped, collapse = ", "), "\n")
  }

  p2 <- DimPlot(seu, reduction = "umap", group.by = "celltype", label = TRUE, label.size = 6) +
    labs(title = "UMAP plot group by CellType") +
    theme(legend.position = "right")

  ggsave(file.path(outpath, paste0(prefix, "umap_celltype.png")), p2,
         width = 12, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "umap_celltype.pdf")), p2,
         width = 12, height = 12)
  cat("\nPlot : ", prefix, "umap_celltype.png(pdf) has been generated successfully !\n", sep = "")

  dedup_markers_df <- unique(celltype_df[c("CellType", "Markers")])
  if (!is.null(celltype_order)) {
    dedup_markers_df <- dedup_markers_df %>%
      mutate(CellType = factor(CellType, levels = celltype_order)) %>%
      filter(!is.na(CellType)) %>%
      arrange(CellType)
  } else {
    dedup_markers_df <- dedup_markers_df %>%
      arrange(CellType)
  }

  marker_long <- dedup_markers_df %>%
    mutate(Markers = strsplit(Markers, ",")) %>%
    unnest(Markers) %>%
    mutate(Markers = trimws(Markers))

  marker_long_dedup <- marker_long %>%
    distinct(Markers, .keep_all = TRUE)

  genes_in_seurat <- rownames(seu[["RNA"]])
  missing_genes <- setdiff(unique(marker_long_dedup$Markers), genes_in_seurat)
  if (length(missing_genes) > 0) {
    message("[DotPlot] ", length(missing_genes),
            " markers NOT found in the Seurat object and will be removed:\n",
            paste(missing_genes, collapse = ", "))
  }
  marker_long_dedup <- marker_long_dedup %>% filter(Markers %in% genes_in_seurat)

  markers_list <- marker_long_dedup %>%
    group_by(CellType) %>%
    summarise(Markers = list(Markers), .groups = "drop")
  markers_list <- markers_list[match(dedup_markers_df$CellType, markers_list$CellType), ]
  markers_list <- setNames(markers_list$Markers, markers_list$CellType)

  seu@meta.data$celltype <- factor(seu@meta.data$celltype,
                                  levels = rev(dedup_markers_df$CellType))

  p3 <- DotPlot(
    seu,
    features = markers_list,
    cols = c("grey60", "red"),
    group.by = "celltype",
    assay = "RNA",
    dot.scale = 8
  ) +
    RotatedAxis() +
    theme(
      strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
      panel.border = element_rect(color = "black", size = 1),
      panel.spacing = unit(2, "mm"),
      strip.text = element_text(margin = margin(b = 2, t = 2, unit = "mm"),
                                face = "bold", size = 12),
      strip.placement = 'outside',
      axis.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 14)
    ) +
    scale_x_discrete("Markers") +
    scale_y_discrete("Cell type")
  ggsave(file.path(outpath, paste0(prefix, "dotplot_celltype.png")), p3,
         width = plot_width, height = plot_height, bg = "white", dpi = 300)
  ggsave(file.path(outpath, paste0(prefix, "dotplot_celltype.pdf")), p3,
         width = plot_width, height = plot_height)
  cat("\nPlot : ", prefix, "dotplot_celltype.png(pdf) has been generated successfully !\n", sep = "")

  saveRDS(seu, file = annotated_rds)
  cat("\nFile : ", prefix, "seu_annotated.rds saved successfully !\n", sep = "")

  write_summary("\n---------- Cell type composition ----------\n")
  cat("\nCells per cell type:\n")
  print(table(seu$celltype, useNA = "ifany"))
  tb_ct <- as.data.frame(table(seu$celltype, useNA = "ifany"), stringsAsFactors = FALSE)
  colnames(tb_ct) <- c("celltype", "cells")
  cat(capture.output(write.table(tb_ct, sep = "\t", row.names = FALSE, quote = FALSE)),
      file = summary_con, sep = "\n")
  cat("\n", file = summary_con)
}

close(summary_con)
cat("\nFile : ", prefix, "Annotation_summary.txt saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")