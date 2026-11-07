#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
# =============================================================================
# Pseudobulk differential expression across cell types.
#
# Aggregate raw counts per (cell type x sample), test the group difference at
# sample level with DESeq2, then draw a volcano plot per cell type.
#
# Example:
#   scRNA_pseudobulk_DE.r -I seu_annotated.rds -CC celltype -SC orig.ident \
#       -GC group -G1 KO -G2 WT -P 01. -O 01_pseudobulk
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
    Pseudobulk differential expression: aggregate raw counts per (cell type x sample),
    test the group difference at sample level with DESeq2, and draw one volcano plot
    per cell type.
  Example: scRNA_pseudobulk_DE.r -I seu_annotated.rds -CC celltype -SC orig.ident -GC group -G1 treat -G2 control -P 01. -O 01_pseudobulk
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
  help = 'Path to an annotated Seurat rds containing a counts layer. (Required)'
)
parser$add_argument(
  "-CC",
  "--celltype_column",
  dest = "celltype_column",
  type = "character",
  required = TRUE,
  help = 'Metadata column holding the cell type; one DE test is run per level. (Required)'
)
parser$add_argument(
  "-SC",
  "--sample_column",
  dest = "sample_column",
  type = "character",
  required = TRUE,
  help = 'Metadata column holding the sample / donor ID (the biological replicate). (Required)'
)
parser$add_argument(
  "-GC",
  "--group_column",
  dest = "group_column",
  type = "character",
  required = TRUE,
  help = 'Metadata column holding the grouping variable to compare. (Required)'
)
parser$add_argument(
  "-G1",
  "--group1",
  dest = "group1",
  type = "character",
  required = TRUE,
  help = 'Group level used as ident.1; log2FC > 0 means higher in this group. (Required)'
)
parser$add_argument(
  "-G2",
  "--group2",
  dest = "group2",
  type = "character",
  required = TRUE,
  help = 'Group level used as ident.2 (reference). (Required)'
)

# Optional
parser$add_argument(
  "-A",
  "--assay",
  dest = "assay",
  type = "character",
  default = "RNA",
  help = 'Assay holding the counts layer. (Optional, default: RNA)'
)
parser$add_argument(
  "-MC",
  "--min_cells",
  dest = "min_cells",
  type = "integer",
  default = 20,
  help = 'Minimum number of cells in a (cell type x sample) unit; units below are dropped. (Optional, default: 20)'
)
parser$add_argument(
  "-MS",
  "--min_samples",
  dest = "min_samples",
  type = "integer",
  default = 3,
  help = 'Minimum number of samples per group required for a cell type to be tested. (Optional, default: 3)'
)
parser$add_argument(
  "-p_val",
  dest = "p_threshold",
  type = "double",
  default = 0.05,
  help = 'Raw p-value threshold for filtering. (Optional, default: 0.05)'
)
parser$add_argument(
  "-p_adj",
  dest = "p_adj_threshold",
  type = "double",
  default = 0.05,
  help = 'Adjusted p-value (FDR) threshold for filtering. (Optional, default: 0.05)'
)
parser$add_argument(
  "-lgFC",
  dest = "log2FC_threshold",
  type = "double",
  default = 1,
  help = 'Absolute log2 fold change threshold for filtering. (Optional, default: 1)'
)
parser$add_argument(
  "-top_n",
  dest = "top_n",
  type = "integer",
  default = 5,
  help = 'Number of top genes labelled per cell type and per side. (Optional, default: 5)'
)
parser$add_argument(
  "-W",
  "--plot_width",
  dest = "plot_width",
  type = "integer",
  default = 18,
  help = 'Width of the volcano plot in inches. (Optional, default: 18)'
)
parser$add_argument(
  "-H",
  "--plot_height",
  dest = "plot_height",
  type = "integer",
  default = 12,
  help = 'Height of the volcano plot in inches. (Optional, default: 12)'
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
celltype_column <- opts$celltype_column
sample_column <- opts$sample_column
group_column <- opts$group_column
group1 <- opts$group1
group2 <- opts$group2
assay <- opts$assay
min_cells <- opts$min_cells
min_samples <- opts$min_samples
p_threshold <- opts$p_threshold
p_adj_threshold <- opts$p_adj_threshold
log2FC_threshold <- opts$log2FC_threshold
top_n <- opts$top_n
plot_width <- opts$plot_width
plot_height <- opts$plot_height
prefix <- opts$prefix
outpath <- opts$outpath

if (!file.exists(input_rds)) {
  stop("Input file does not exist: ", input_rds)
}
if (group1 == group2) {
  stop("-G1 and -G2 must be different group levels.")
}

message("\n=====================scRNA Pseudobulk DE======================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2026 OYJH. All Rights Reserved.")
message("===============================================================\n")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if (!require(tidyverse, quietly = TRUE)) install.packages("tidyverse")
  if (!require(Seurat, quietly = TRUE)) install.packages("Seurat")
  if (!require(Matrix, quietly = TRUE)) install.packages("Matrix")
  if (!require(DESeq2, quietly = TRUE)) BiocManager::install("DESeq2")
  if (!require(ggrepel, quietly = TRUE)) install.packages("ggrepel")
  if (!require(RColorBrewer, quietly = TRUE)) install.packages("RColorBrewer")
  if (!require(scales, quietly = TRUE)) install.packages("scales")
})

# 3. Load data and validate metadata columns
cat(">>> Loading Seurat object\n")
seo <- readRDS(input_rds)
md <- seo@meta.data
required_cols <- c(celltype_column, sample_column, group_column)
missing_cols <- required_cols[!required_cols %in% colnames(md)]
if (length(missing_cols) > 0) {
  stop("Column(s) not found in metadata: ", paste(missing_cols, collapse = ", "))
}
cat(">>> Object: ", ncol(seo), " cells x ", nrow(seo), " genes\n", sep = "")

celltype_vec <- as.character(md[[celltype_column]])
sample_vec <- as.character(md[[sample_column]])
group_vec <- as.character(md[[group_column]])

group_levels_present <- unique(group_vec)
if (!group1 %in% group_levels_present || !group2 %in% group_levels_present) {
  stop(
    "Group level(s) not found in column '", group_column, "'. Present: ",
    paste(group_levels_present, collapse = ", ")
  )
}

# 4. Aggregate raw counts per (cell type x sample)
assay_obj <- seo[[assay]]
count_layers <- grep("^counts", Layers(assay_obj), value = TRUE)
if (length(count_layers) > 1) {
  cat(">>> Joining split count layers: ", paste(count_layers, collapse = ", "), "\n", sep = "")
  assay_obj <- JoinLayers(assay_obj, layers = "counts")
}
if (!"counts" %in% Layers(assay_obj)) {
  stop("No 'counts' layer found in assay '", assay, "'. Available: ",
       paste(Layers(assay_obj), collapse = ", "))
}
counts_all <- LayerData(assay_obj, layer = "counts")
cat(">>> Counts matrix: ", nrow(counts_all), " x ", ncol(counts_all), "\n", sep = "")

unit_key <- paste(celltype_vec, sample_vec, sep = "||")
unit_levels <- unique(unit_key)
unit_index <- match(unit_key, unit_levels)
if (length(unit_index) != ncol(counts_all)) {
  stop("Internal error: cell count (", ncol(counts_all), ") does not match metadata rows (",
       length(unit_index), ").")
}
celltypes_of_unit <- vapply(strsplit(unit_levels, "||", fixed = TRUE), `[`, character(1), 1)
samples_of_unit <- vapply(strsplit(unit_levels, "||", fixed = TRUE), `[`, character(1), 2)

indicator <- sparseMatrix(
  i = seq_len(ncol(counts_all)),
  j = unit_index,
  x = 1,
  dims = c(ncol(counts_all), length(unit_levels)),
  repr = "C"
)
pb_counts <- as.matrix(counts_all %*% indicator)
colnames(pb_counts) <- unit_levels
unit_cellcount <- as.integer(tabulate(unit_index, nbins = length(unit_levels)))
pb_group <- group_vec[match(unit_levels, unit_key)]
cat(">>> Pseudobulk units: ", ncol(pb_counts), " (cell type x sample)\n", sep = "")

# 5. Per cell type: DESeq2 on the pseudobulk matrix
cat(">>> Running DESeq2 per cell type (group1 = ", group1, ", group2 = ", group2, ")\n", sep = "")

celltypes <- sort(unique(celltype_vec))
results_list <- list()
unit_report <- data.frame()

for (ct in celltypes) {
  sel <- which(celltypes_of_unit == ct)
  keep <- sel[unit_cellcount[sel] >= min_cells]
  if (length(keep) == 0) {
    cat(">>> Skip ", ct, ": no unit reaches ", min_cells, " cells\n", sep = "")
    next
  }

  mat <- pb_counts[, keep, drop = FALSE]
  grp <- factor(pb_group[keep], levels = c(group2, group1))
  n1 <- sum(grp == group1)
  n2 <- sum(grp == group2)

  unit_report <- rbind(
    unit_report,
    data.frame(
      celltype = ct, n_units = length(keep), n_group1 = n1, n_group2 = n2,
      min_cells = min(unit_cellcount[keep]), max_cells = max(unit_cellcount[keep]),
      stringsAsFactors = FALSE
    )
  )
  if (n1 < min_samples || n2 < min_samples) {
    cat(">>> Skip ", ct, ": group1 units = ", n1, ", group2 units = ", n2,
        " (need ", min_samples, " each)\n", sep = "")
    next
  }

  res <- try({
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = mat,
      colData = data.frame(grp = grp, row.names = colnames(mat)),
      design = ~ grp
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res_df <- as.data.frame(DESeq2::results(dds, contrast = c("grp", group1, group2)))

    expressed <- mat > 0
    is_g1 <- grp == group1
    is_g2 <- grp == group2
    pct1 <- rowMeans(expressed[, is_g1, drop = FALSE])
    pct2 <- rowMeans(expressed[, is_g2, drop = FALSE])

    data.frame(
      gene = rownames(res_df),
      p_val = res_df$pvalue,
      avg_log2FC = res_df$log2FoldChange,
      pct.1 = pct1[match(rownames(res_df), rownames(mat))],
      pct.2 = pct2[match(rownames(res_df), rownames(mat))],
      p_val_adj = res_df$padj,
      cluster = ct,
      stringsAsFactors = FALSE
    )
  }, silent = TRUE)

  if (inherits(res, "try-error")) {
    cat(">>> Skip ", ct, ": DESeq2 failed (", conditionMessage(attr(res, "condition")), ")\n", sep = "")
    next
  }

  res <- res[order(res$avg_log2FC, decreasing = TRUE), ]
  results_list[[ct]] <- res
  cat(">>> ", ct, ": ", nrow(res), " genes tested (", n1, " vs ", n2, " units)\n", sep = "")
}

if (length(results_list) == 0) {
  stop("No cell type passed the filters - nothing to test.")
}

# 6. Combine and write tables
combined <- do.call(rbind, results_list)
rownames(combined) <- NULL
combined <- combined[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster")]

file_n <- 0
plot_n <- 0

raw_file <- file.path(outpath, paste0(prefix, group1, "_vs_", group2, "_DE_raw.csv"))
write.csv(combined, file = raw_file, row.names = FALSE, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", basename(raw_file), " saved successfully !\n", sep = "")

filtered <- combined %>%
  filter(p_val < p_threshold & p_val_adj < p_adj_threshold & abs(avg_log2FC) > log2FC_threshold) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  ungroup()
filtered_file <- file.path(outpath, paste0(prefix, group1, "_vs_", group2, "_DE_filtered.csv"))
write.csv(filtered, file = filtered_file, row.names = FALSE, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", basename(filtered_file), " saved successfully !\n", sep = "")

unit_file <- file.path(outpath, paste0(prefix, group1, "_vs_", group2, "_pseudobulk_units.csv"))
write.csv(unit_report, file = unit_file, row.names = FALSE, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", basename(unit_file), " saved successfully !\n", sep = "")

# 7. Volcano plot
cat(">>> Drawing volcano plot\n")
df <- filtered
if (nrow(df) == 0) {
  cat(">>> No gene passed the thresholds - skipping the volcano plot.\n")
} else {
  df$label <- ifelse(df$avg_log2FC > 0, "sigUp", "sigDown")
  top_df <- df %>%
    group_by(cluster) %>%
    group_by(sign_group = ifelse(avg_log2FC > 0, "positive", "negative"), .add = TRUE) %>%
    arrange(ifelse(sign_group == "positive", -avg_log2FC, avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    dplyr::select(-sign_group)
  non_top_df <- df %>% anti_join(top_df, by = names(df))

  df$cluster <- as.factor(df$cluster)
  cluster_order <- levels(droplevels(df$cluster))
  top_log2FC <- df %>%
    filter(cluster %in% cluster_order) %>%
    group_by(cluster) %>%
    summarise(
      max = max(avg_log2FC, na.rm = TRUE),
      min = min(avg_log2FC, na.rm = TRUE)
    ) %>%
    arrange(match(cluster, cluster_order)) %>%
    mutate(
      same_sign = sign(max) == sign(min),
      max_adj = ifelse(same_sign, ifelse(abs(max) >= abs(min), max, 0), max),
      min_adj = ifelse(same_sign, ifelse(abs(min) > abs(max), min, 0), min)
    )
  bar_top <- data.frame(x = top_log2FC$cluster, y = top_log2FC$max_adj)
  bar_bot <- data.frame(x = top_log2FC$cluster, y = top_log2FC$min_adj)

  p1 <- ggplot() +
    geom_col(data = bar_top, mapping = aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.5) +
    geom_col(data = bar_bot, mapping = aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.5) +
    scale_y_continuous(
      trans = "pseudo_log",
      breaks = function(lims) {
        auto_breaks <- pretty_breaks()(lims)
        sort(unique(c(auto_breaks, -log2FC_threshold, log2FC_threshold)))
      }
    )

  top_df$label <- factor(top_df$label, levels = c("sigUp", "sigDown"))
  top_df$cluster <- factor(top_df$cluster, levels = levels(bar_top$x))
  non_top_df$cluster <- factor(non_top_df$cluster, levels = levels(bar_top$x))
  p2 <- p1 +
    geom_jitter(data = top_df, aes(x = cluster, y = avg_log2FC, color = label), size = 2, width = 0.4) +
    geom_jitter(data = non_top_df, aes(x = cluster, y = avg_log2FC, color = label), size = 1, width = 0.4) +
    scale_color_manual(values = c("sigUp" = "red", "sigDown" = "blue"))

  dfcol <- data.frame(x = c(1:length(unique(df$cluster))), y = 0, label = cluster_order)
  n_clusters <- length(unique(df$cluster))
  if (n_clusters <= 66) {
    col_1 <- brewer.pal(brewer.pal.info["Set1", "maxcolors"], "Set1")
    col_2 <- brewer.pal(brewer.pal.info["Dark2", "maxcolors"], "Dark2")
    col_3 <- brewer.pal(brewer.pal.info["Pastel1", "maxcolors"], "Pastel1")
    col_4 <- brewer.pal(brewer.pal.info["Paired", "maxcolors"], "Paired")
    col_5 <- brewer.pal(brewer.pal.info["Set3", "maxcolors"], "Set3")
    col_6 <- brewer.pal(brewer.pal.info["Accent", "maxcolors"], "Accent")
    col_7 <- brewer.pal(brewer.pal.info["Pastel2", "maxcolors"], "Pastel2")
    mycols <- c(col_1, col_2, col_3, col_4, col_5, col_6, col_7)[1:n_clusters]
  } else {
    warning("Large number of clusters detected (", n_clusters, "). Using viridis color palette.")
    suppressPackageStartupMessages({
      if (!require(viridis, quietly = TRUE)) install.packages("viridis")
    })
    mycols <- viridis(n_clusters, option = "D")
  }
  p3 <- p2 + geom_tile(
    data = dfcol, aes(x = x, y = y), height = log2FC_threshold, color = "black",
    fill = mycols, alpha = 1, show.legend = FALSE
  )
  p4 <- p3 +
    geom_text_repel(
      data = top_df, aes(x = cluster, y = avg_log2FC, label = gene),
      force = 1, point.size = NA, point.padding = 0.3, box.padding = 0.5,
      arrow = NULL, segment.color = NA, fontface = "bold"
    )
  p5 <- p4 +
    labs(x = celltype_column, y = "average log2FC (pseudo_log)") +
    geom_text(data = dfcol, aes(x = x, y = y, label = label), size = 5, color = "white", fontface = "bold")
  p6 <- p5 +
    labs(title = paste0("Volcano Plot of DEGs in ", group1, " vs ", group2)) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, color = "black", face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      axis.line.y = element_line(color = "black", linewidth = 1.2),
      axis.line.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 14, color = "black"),
      panel.grid = element_blank(),
      legend.position = "top",
      legend.direction = "vertical",
      legend.justification = c(1, 0),
      legend.text = element_text(size = 15),
      legend.title = element_blank()
    )

  ggsave(file.path(outpath, paste0(prefix, group1, "_vs_", group2, "_DE_volcano.png")),
         plot = p6, width = plot_width, height = plot_height, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, group1, "_vs_", group2, "_DE_volcano.pdf")),
         plot = p6, width = plot_width, height = plot_height)
  plot_n <- plot_n + 1
  cat("\nPlot ", plot_n, " : ", paste0(prefix, group1, "_vs_", group2, "_DE_volcano.png (pdf)"),
      " has been generated successfully !\n", sep = "")
}

message("\n=========================All Done !============================\n")
