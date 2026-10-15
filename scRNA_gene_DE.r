#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
rm(list = ls())
options(timeout = 1000)

# 1. Parameter settings
suppressPackageStartupMessages({
  if (!require(argparse, quietly = TRUE)) install.packages("argparse")
})

parser <- ArgumentParser(
  description = paste0(
'
    Per-sample x per-celltype average expression of ONE target gene from an
    annotated Seurat object, with three visualizations:

      CSV      sample x celltype average expression (missing combos = 0)
      Heatmap  per-celltype z-score of the sample x celltype table
      Barplot  absolute average expression per sample, one per -CT celltype
      Boxplot  sample-level boxplot of two groups (-GC), one per -CT celltype,
               with points, a Wilcoxon p value and a significance bar

    Cells are NOT used as replicates in the boxplot: each point is one sample.

  Example: scRNA_gene_DE -I 05_annotation/05.seu_annotated.rds -G FTH1 \\
                          -CT Neut,Act-Fib -GO Normal,BAS -O 07_geneDE
'
  )
)
parser$add_argument(
  "-I", "--input",
  dest = "input_file", type = "character", required = TRUE,
  help = "Annotated Seurat rds (e.g. seu_annotated.rds). (Required)"
)
parser$add_argument(
  "-G", "--gene",
  dest = "gene", type = "character", required = TRUE,
  help = "Single target gene symbol (exactly one). (Required)"
)
parser$add_argument(
  "-S", "--sample_col",
  dest = "sample_col", type = "character", default = "sample",
  help = "Metadata column holding the sample (default: sample)."
)
parser$add_argument(
  "-CTC", "--celltype_col",
  dest = "celltype_col", type = "character", default = "celltype",
  help = "Metadata column holding the celltype (default: celltype)."
)
parser$add_argument(
  "-GC", "--group_col",
  dest = "group_col", type = "character", default = "group",
  help = "Metadata column holding the group for boxplots (default: group)."
)
parser$add_argument(
  "-L", "--layer",
  dest = "layer", type = "character", default = "data",
  help = "Seurat layer used for averaging (default: data, log-normalized)."
)
parser$add_argument(
  "-CT", "--celltypes",
  dest = "celltypes", type = "character", default = "",
  help = "Comma-separated celltypes to draw barplots + boxplots for (e.g. Neut,Act-Fib). Empty = CSV + heatmap only."
)
parser$add_argument(
  "-SO", "--sample_order",
  dest = "sample_order", type = "character", default = "",
  help = "Comma-separated sample order for BARPLOTS, left to right (e.g. Ctrl1,Ctrl2,Tr1). Empty = table order. Barplots only; heatmap and boxplots are unaffected."
)
parser$add_argument(
  "-GO", "--group_order",
  dest = "group_order", type = "character", default = "",
  help = "Comma-separated group order for BOXPLOTS, left to right (empty = sorted). The first group is the reference and the one-sided test checks the second > the first."
)
parser$add_argument(
  "--alt",
  dest = "alternative", type = "character", default = "auto",
  choices = c("auto", "two.sided", "greater", "less"),
  help = "Boxplot test direction (default: auto). 'auto' = one-sided 'less': H1 = group2 > group1 (right-most group higher). 'two.sided' = non-directional."
)
parser$add_argument(
  "--signif_bar",
  dest = "signif_bar", type = "character", default = "auto",
  choices = c("auto", "none"),
  help = "Draw the boxplot significance bar (default: auto). The bar always sits above the tallest data point so no sample is clipped. 'none' omits it."
)
parser$add_argument(
  "--bar_gap",
  dest = "bar_gap", type = "numeric", default = 0.5,
  help = "Boxplot bar clearance above the tallest data point, as a multiple of the y-axis padding (default: 0.5)."
)
parser$add_argument(
  "--col_order", dest = "col_order", type = "character", default = "hclust",
  choices = c("hclust", "csv"),
  help = paste0("Heatmap celltype column order (default: hclust). ",
                "hclust = cluster by per-sample expression profile so that ",
                "celltypes with similar patterns sit next to each other ",
                "(columns are z-scored first, then Euclidean distance + ",
                "complete linkage); csv = keep the table order. Only the ",
                "heatmap is affected - the CSV, barplots and boxplots keep ",
                "the table order.")
)
parser$add_argument(
  "--box_coef",
  # NOTE: the default must be a PYTHON literal (argparse is Python-backed), so
  # Inf cannot be used here - R's Inf would raise NameError. 1e9 is
  # effectively Inf for an IQR multiple.
  dest = "box_coef", type = "numeric", default = 1.5,
  help = paste0("Boxplot whisker length as a multiple of the IQR (default: 1.5, ",
                "the classic Tukey rule). The whisker reaches the most extreme ",
                "DATA POINT inside the fence, never the fence value itself. ",
                "Samples outside the fence are drawn as black diamond markers ",
                "and excluded from the jitter layer so each sample appears once. ",
                "Use 1e9 (effectively Inf) to make whiskers span the true ",
                "min/max and flag nothing as an outlier. NOTE: ggplot's coef = 0 ",
                "is NOT base R's range = 0 - it makes the fence equal to Q1/Q3, ",
                "so every point outside the box becomes an outlier and the ",
                "whiskers collapse to zero length.")
)
parser$add_argument(
  "-W", "--width",
  dest = "plot_width", type = "numeric", default = 14,
  help = "Heatmap width in inches (default: 14)."
)
parser$add_argument(
  "-H", "--height",
  dest = "plot_height", type = "numeric", default = 10,
  help = "Heatmap height in inches (default: 10)."
)
parser$add_argument(
  "-BW", "--bar_width",
  dest = "bar_width", type = "numeric", default = 8,
  help = "Barplot width in inches (default: 8)."
)
parser$add_argument(
  "-BH", "--bar_height",
  dest = "bar_height", type = "numeric", default = 6,
  help = "Barplot height in inches (default: 6)."
)
parser$add_argument(
  "-QW", "--box_width",
  dest = "box_width", type = "numeric", default = 5,
  help = "Boxplot width in inches (default: 5)."
)
parser$add_argument(
  "-QH", "--box_height",
  dest = "box_height", type = "numeric", default = 6,
  help = "Boxplot height in inches (default: 6)."
)
parser$add_argument(
  "-P", "--prefix",
  dest = "prefix", type = "character", default = "",
  help = "Prefix for output files (Optional, default: none)."
)
parser$add_argument(
  "--no-plot",
  dest = "no_plot", action = "store_true", default = FALSE,
  help = "Skip all visualizations (CSV output only)."
)
parser$add_argument(
  "-O", "--outpath",
  dest = "outpath", type = "character", default = getwd(),
  help = "Output directory (Optional, default: current directory)."
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help(); quit(status = 1)
}
opts <- parser$parse_args()

input_file   <- opts$input_file
gene         <- opts$gene
sample_col   <- opts$sample_col
celltype_col <- opts$celltype_col
group_col    <- opts$group_col
layer        <- opts$layer
prefix       <- opts$prefix
outpath      <- opts$outpath
no_plot      <- opts$no_plot
plot_width   <- opts$plot_width
plot_height  <- opts$plot_height
bar_width    <- opts$bar_width
bar_height   <- opts$bar_height
box_width    <- opts$box_width
box_height   <- opts$box_height
alternative  <- opts$alternative
signif_bar   <- opts$signif_bar
bar_gap      <- opts$bar_gap
box_coef     <- opts$box_coef

if (grepl("[, \t\n]", gene)) {
  stop("Error: -G/--gene accepts exactly ONE gene symbol; got a list. ",
       "Remove commas/spaces from the value.")
}

# Celltypes to plot (empty = CSV + heatmap only)
celltypes_to_plot <- trimws(unlist(strsplit(opts$celltypes, "[,;]+")))
celltypes_to_plot <- celltypes_to_plot[nchar(celltypes_to_plot) > 0]

# Barplot sample order (empty = table order)
bar_sample_order <- trimws(unlist(strsplit(opts$sample_order, "[,;]+")))
bar_sample_order <- bar_sample_order[nchar(bar_sample_order) > 0]

# Boxplot group order (empty = sorted)
group_order <- trimws(unlist(strsplit(opts$group_order, "[,;]+")))
group_order <- group_order[nchar(group_order) > 0]

message("\n==================scRNA GeneDE Sample x Celltype==================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Input RDS   : ", input_file)
message("Gene        : ", gene)
message("Sample col  : ", sample_col, " | Celltype col: ", celltype_col,
        " | Group col: ", group_col)
message("Layer       : ", layer)
if (length(celltypes_to_plot) > 0) {
  message("Barplot/Boxplot celltypes: ", paste(celltypes_to_plot, collapse = ", "))
} else {
  message("Barplot/Boxplot celltypes: (none, CSV + heatmap only)")
}
message("Copyright \u00A9 2026 OYJH. All Rights Reserved.")
message("\n==================================================================\n")

if (!dir.exists(outpath)) dir.create(outpath, recursive = TRUE)

# Summary file
summary_file <- file.path(outpath, paste0(prefix, "GeneDE_summary.txt"))
summary_con  <- file(summary_file, open = "wt")
write_summary <- function(...) {
  args_list <- list(...)
  multi <- vapply(args_list, length, integer(1)) > 1
  if (any(multi)) {
    which_arg <- which(multi)[1]
    stop("write_summary arg #", which_arg, " is a length-",
         length(args_list[[which_arg]]),
         " vector which would be collapsed into a single line by paste0. ",
         "Use cat(..., sep=\"\\n\", file = summary_con) instead for multi-line content.")
  }
  cat(paste0(..., collapse = ""), file = summary_con)
}
write_summary("===================scRNA GeneDE Summary===================\n")
write_summary("Current time : ", as.character(Sys.time()), "\n")
write_summary("Input RDS    : ", input_file, "\n")
write_summary("Gene         : ", gene, "\n")
write_summary("Sample col   : ", sample_col, " | Celltype col: ", celltype_col, "\n")
write_summary("Group col    : ", group_col, "\n")
write_summary("Layer        : ", layer, "\n")

# 2. Load packages and data
suppressPackageStartupMessages({
  if (!require(tidyverse, quietly = TRUE)) install.packages("tidyverse")
  if (!require(Seurat,    quietly = TRUE)) install.packages("Seurat")
  if (!require(ggplot2,   quietly = TRUE)) install.packages("ggplot2")
  if (!require(reshape2,  quietly = TRUE)) install.packages("reshape2")
})

seu <- readRDS(input_file)
cat("\nLoaded Seurat object from:", input_file, "\n")
cat("Cells:", ncol(seu), " | Genes:", nrow(seu), " | Layers:",
    paste(Layers(seu), collapse = ", "), "\n")

for (cl in c(sample_col, celltype_col)) {
  if (!(cl %in% colnames(seu@meta.data))) {
    stop("Error: metadata column '", cl, "' not found. Available: ",
         paste(colnames(seu@meta.data), collapse = ", "))
  }
}
if (!(layer %in% Layers(seu))) {
  stop("Error: layer '", layer, "' not found in the Seurat object. Available: ",
       paste(Layers(seu), collapse = ", "))
}
if (!(gene %in% rownames(seu))) {
  stop("Error: gene '", gene, "' not found in the Seurat object.")
}

# Keep the ORIGINAL ordering semantics: for a factor column, sort() returns
# level order (not alphabetical). Converting to character first would silently
# reorder the CSV columns and therefore the heatmap columns.
sample_levels   <- as.character(sort(unique(seu[[sample_col,   drop = TRUE]])))
celltype_levels <- as.character(sort(unique(seu[[celltype_col, drop = TRUE]])))
cat("Samples (", length(sample_levels), "): ",
    paste(sample_levels, collapse = ", "), "\n", sep = "")
cat("Celltypes (", length(celltype_levels), "): ",
    paste(celltype_levels, collapse = ", "), "\n", sep = "")

write_summary("Gene present : yes\n")
write_summary("Samples      : ", length(sample_levels), "\n")
write_summary("Celltypes    : ", length(celltype_levels), "\n")

# 3. Compute sample x celltype mean (layer)
cat("\n>>> Computing sample x celltype mean expression (layer:", layer, ") ...\n")
avg_mat <- AverageExpression(seu, features = gene,
                             group.by = c(sample_col, celltype_col),
                             layer = layer, verbose = FALSE)$RNA

df_long <- reshape2::melt(as.matrix(avg_mat),
                          varnames = c("gene", "group"),
                          value.name = "avg_expression")

# Seurat joins group.by levels with "_", and sample/celltype names may contain
# "_" themselves, so split greedily against the known level lists.
parts <- matrix(NA_character_, nrow = nrow(df_long), ncol = 2,
                dimnames = list(NULL, c(sample_col, celltype_col)))
for (i in seq_len(nrow(df_long))) {
  g <- as.character(df_long$group[i])
  matched <- FALSE
  for (s in sample_levels) {
    if (startsWith(g, paste0(s, "_"))) {
      rest <- substring(g, nchar(s) + 2L)
      if (rest %in% celltype_levels) {
        parts[i, 1] <- s
        parts[i, 2] <- rest
        matched <- TRUE
        break
      }
    }
  }
  if (!matched) {
    # fallback: split on the LAST underscore
    idx <- max(gregexpr("_", g, fixed = TRUE)[[1]])
    if (idx > 0) {
      parts[i, 1] <- substring(g, 1, idx - 1L)
      parts[i, 2] <- substring(g, idx + 1L)
    } else {
      parts[i, 1] <- g
      parts[i, 2] <- NA_character_
    }
  }
}
df_long[[sample_col]]   <- parts[, 1]
df_long[[celltype_col]] <- parts[, 2]
df_long$group <- NULL

# Fill missing sample x celltype combos with 0
df_full <- df_long
df_full[["avg_expression"]] <- as.numeric(df_full[["avg_expression"]])
df_full <- tidyr::complete(
  df_full,
  gene = unique(gene),
  .data[[sample_col]],
  .data[[celltype_col]],
  fill = list(avg_expression = 0)
)

# Wide-form: rows = sample, cols = celltype
sub <- df_full[, c(sample_col, celltype_col, "avg_expression")]
sub[["key"]] <- paste(sub[[sample_col]], sub[[celltype_col]], sep = "\r")
full_keys <- as.vector(outer(sample_levels, celltype_levels,
                             function(s, c) paste(s, c, sep = "\r")))
val_vec <- sub[["avg_expression"]][match(full_keys, sub[["key"]])]
val_vec[is.na(val_vec)] <- 0
mat <- matrix(val_vec, nrow = length(sample_levels),
              ncol = length(celltype_levels), byrow = FALSE,
              dimnames = list(sample_levels, celltype_levels))
wide <- data.frame(sample = sample_levels, mat, check.names = FALSE,
                   stringsAsFactors = FALSE)
names(wide)[1] <- sample_col
# Sanitize only characters that would break CSV parsing (comma, quote, newline).
# Do NOT use make.names(): it rewrites "-" to "." and would turn "Act-Fib" into
# "Act.Fib", corrupting celltype names in the CSV and every plot label.
colnames(wide) <- gsub('[",\r\n]', "_", colnames(wide))
# Optional column clustering for the heatmap. Cluster on the z-scored matrix
# so high-expression columns (e.g. Neut) do not dominate the distance metric.
# Runs AFTER wide exists so we have the per-celltype mean values.
if (opts$col_order == "hclust") {
  z_mat <- scale(as.matrix(wide[, celltype_levels, drop = FALSE]))
  d     <- dist(t(z_mat), method = "euclidean")
  hc    <- stats::hclust(d, method = "complete")
  celltype_levels <- hc$labels[hc$order]
  cat("Heatmap col order (hclust, top of tree -> bottom):\n  ",
      paste(celltype_levels, collapse = ", "), "\n", sep = "")
}
csv_out <- file.path(outpath, paste0(prefix, gene, "_sample_celltype_avg.csv"))
write.csv(wide, file = csv_out, row.names = FALSE, quote = FALSE)
cat("\nFile 1 : ", prefix, gene, "_sample_celltype_avg.csv saved successfully !\n", sep = "")
write_summary("\n---------- Output ----------\n")
write_summary("CSV          : ", basename(csv_out), " (", length(sample_levels),
              " samples x ", length(celltype_levels), " celltypes)\n")

# 4. Sample -> group map (for boxplots), taken from the rds metadata
group_map <- unique(data.frame(
  sample = as.character(seu[[sample_col, drop = TRUE]]),
  group  = as.character(seu[[group_col,  drop = TRUE]]),
  stringsAsFactors = FALSE
))
group_map <- group_map[!duplicated(group_map$sample), , drop = FALSE]

if (!no_plot) {
  cat("\n>>> Plotting ...\n")
  # geom_jitter() uses random offsets; fix the seed so a re-run reproduces
  # byte-identical figures (otherwise every run gives slightly different
  # point positions).
  set.seed(42)

  # ---- 4a. Heatmap: per-celltype z-score of the sample x celltype table ----
  z_mat <- as.matrix(wide[, -1, drop = FALSE])
  col_means <- colMeans(z_mat, na.rm = TRUE)
  col_sds   <- apply(z_mat, 2, sd, na.rm = TRUE)
  col_sds[col_sds == 0 | is.na(col_sds)] <- 1
  z_mat <- (z_mat - matrix(col_means, nrow = nrow(z_mat), ncol = ncol(z_mat), byrow = TRUE)) /
           matrix(col_sds,   nrow = nrow(z_mat), ncol = ncol(z_mat), byrow = TRUE)
  z_wide <- data.frame(wide[, 1, drop = FALSE], z_mat,
                       check.names = FALSE, stringsAsFactors = FALSE)

  hm_long <- reshape2::melt(z_wide,
                            id.vars = sample_col,
                            variable.name = celltype_col,
                            value.name = "z_score")
  hm_long[[sample_col]]   <- factor(hm_long[[sample_col]],   levels = sample_levels)
  hm_long[[celltype_col]] <- factor(hm_long[[celltype_col]], levels = celltype_levels)
  hm_long$x <- as.numeric(hm_long[[celltype_col]])
  # Reverse the sample POSITION (not the label text) so the first sample in
  # sample_levels sits at the TOP. Reversing only the labels would mislabel
  # every row (y=1 would be drawn with the last sample's name).
  hm_long$y <- as.numeric(factor(as.character(hm_long[[sample_col]]),
                                 levels = rev(sample_levels)))

  p_hm <- ggplot2::ggplot(hm_long, ggplot2::aes(x = x, y = y, fill = z_score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "navy", mid = "white", high = "firebrick3",
      midpoint = 0, limits = c(-3, 3), name = "Expression"
    ) +
    ggplot2::scale_x_continuous(breaks = seq_along(celltype_levels),
                                labels = celltype_levels, position = "top",
                                expand = ggplot2::expansion(0)) +
    ggplot2::scale_y_continuous(breaks = seq_along(sample_levels),
                                labels = rev(sample_levels),
                                expand = ggplot2::expansion(0)) +
    ggplot2::labs(title = paste0(gene, " Average Expression"), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 0, vjust = 0.05,
                                          size = 9, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.ticks = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(0, "pt"),
      legend.text  = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_text(face = "bold"),
      plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      panel.grid   = ggplot2::element_blank(),
      plot.margin  = ggplot2::margin(t = 5, r = 4, b = 4, l = 4, unit = "mm")
    )
  ggplot2::ggsave(file.path(outpath, paste0(prefix, gene, "_heatmap.png")),
                  p_hm, width = plot_width, height = plot_height,
                  dpi = 300, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(outpath, paste0(prefix, gene, "_heatmap.pdf")),
                  p_hm, width = plot_width, height = plot_height,
                  limitsize = FALSE)
  cat("Plot 1 : ", prefix, gene, "_heatmap.png(pdf) saved successfully !\n", sep = "")
  write_summary("Heatmap      : ", basename(paste0(prefix, gene, "_heatmap.png")),
                " (per-celltype z-score, [-3, 3], legend label = Expression, ",
                plot_width, "x", plot_height, " in)\n")

  # ---- 4b. Barplots + 4c. Boxplots per requested celltype ----
  if (length(celltypes_to_plot) > 0) {
    missing_ct <- setdiff(celltypes_to_plot, colnames(wide))
    if (length(missing_ct) > 0) {
      warning("Celltype(s) not found, skipped: ", paste(missing_ct, collapse = ", "),
              "\n  Available: ", paste(colnames(wide)[-1], collapse = ", "))
    }
    cts_use <- intersect(celltypes_to_plot, colnames(wide))
    if (length(cts_use) == 0) stop("Error: none of the requested celltypes exist.")

    # Barplot sample order: explicit order first, remaining samples appended
    bar_levels <- sample_levels
    if (length(bar_sample_order) > 0) {
      unknown_s <- setdiff(bar_sample_order, sample_levels)
      if (length(unknown_s) > 0) {
        warning("Sample(s) in -SO not found, ignored: ",
                paste(unknown_s, collapse = ", "))
      }
      rest <- setdiff(sample_levels, bar_sample_order)
      bar_levels <- c(intersect(bar_sample_order, sample_levels), rest)
      if (length(rest) > 0) {
        warning("Sample(s) not listed in -SO were appended at the right: ",
                paste(rest, collapse = ", "))
      }
    }
    cat("Barplot sample order (left -> right): ",
        paste(bar_levels, collapse = ", "), "\n", sep = "")

    # Boxplot long table: sample x celltype -> group
    box_long_list <- list()
    for (ct in cts_use) {
      box_long_list[[length(box_long_list) + 1]] <- data.frame(
        sample = sample_levels, celltype = ct,
        value = as.numeric(wide[[ct]]), stringsAsFactors = FALSE
      )
    }
    box_long <- do.call(rbind, box_long_list)
    box_long <- merge(box_long, group_map, by = "sample", all.x = TRUE)
    if (any(is.na(box_long$group))) {
      bad <- unique(box_long$sample[is.na(box_long$group)])
      warning("Sample(s) missing from the group map (dropped): ",
              paste(bad, collapse = ", "))
      box_long <- box_long[!is.na(box_long$group), , drop = FALSE]
    }
    groups_present <- unique(box_long$group)
    if (length(group_order) > 0) {
      unknown_g <- setdiff(group_order, groups_present)
      if (length(unknown_g) > 0) {
        warning("Group(s) in -GO not found, ignored: ",
                paste(unknown_g, collapse = ", "))
      }
      groups_use <- c(intersect(group_order, groups_present),
                      setdiff(groups_present, group_order))
    } else {
      groups_use <- sort(groups_present)
    }
    box_long$group <- factor(box_long$group, levels = groups_use)
    cat("Groups (left -> right): ", paste(groups_use, collapse = ", "), "\n", sep = "")
    cat("Cells are NOT used as replicates: each boxplot point is one sample.\n")

    stats_rows <- list()
    for (ct in cts_use) {
      safe_ct <- gsub("[/\\\\:*?\"<>|,]", "_", ct)

      # ---------- Barplot: absolute average expression per sample ----------
      df_ct <- data.frame(sample = sample_levels, value = wide[[ct]],
                          stringsAsFactors = FALSE)
      df_ct$sample <- factor(df_ct$sample, levels = bar_levels)
      p_bar <- ggplot2::ggplot(df_ct, ggplot2::aes(x = sample, y = value)) +
        ggplot2::geom_col(fill = "steelblue", width = 0.7) +
        # No bottom padding: y = 0 sits exactly on the x axis so bars touch it.
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::labs(title = paste0(gene, " Average Expression in ", ct),
                      x = NULL, y = "Average expression") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(size = 10, face = "bold",
                                               colour = "black", angle = 45, hjust = 1),
          axis.text.y  = ggplot2::element_text(size = 10, face = "bold", colour = "black"),
          axis.title.y = ggplot2::element_text(size = 11, face = "bold", colour = "black",
                                               margin = ggplot2::margin(r = 14)),
          axis.line    = ggplot2::element_line(colour = "black", linewidth = 0.6),
          axis.ticks   = ggplot2::element_line(colour = "black", linewidth = 0.6),
          plot.title   = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
          panel.grid.minor = ggplot2::element_blank(),
          plot.margin  = ggplot2::margin(t = 5, r = 4, b = 4, l = 4, unit = "mm")
        )
      f_bar_png <- file.path(outpath, paste0(prefix, gene, "_", safe_ct, "_barplot.png"))
      f_bar_pdf <- file.path(outpath, paste0(prefix, gene, "_", safe_ct, "_barplot.pdf"))
      ggplot2::ggsave(f_bar_png, p_bar, width = bar_width, height = bar_height,
                      dpi = 300, bg = "white", limitsize = FALSE)
      ggplot2::ggsave(f_bar_pdf, p_bar, width = bar_width, height = bar_height,
                      limitsize = FALSE)
      cat("Barplot : ", basename(f_bar_png), "(png/pdf) saved successfully !\n", sep = "")

      # ---------- Boxplot: one point per sample, two groups ----------
      d <- box_long[box_long$celltype == ct, , drop = FALSE]
      n_per_group <- sapply(groups_use, function(g) sum(d$group == g))
      alt <- if (alternative == "auto") "less" else alternative
      if (all(n_per_group >= 2) && length(groups_use) == 2) {
        # exact = TRUE handles small samples (n1 + n2 <= 50). "auto" = one-sided
        # "less" because a typical hypothesis is "the right-side group is higher".
        wt <- suppressWarnings(wilcox.test(value ~ group, data = d,
                                           alternative = alt, exact = TRUE))
        pv <- wt$p.value
      } else {
        pv <- NA_real_
      }
      p_lab <- if (is.na(pv)) "p = NA"
               else if (pv < 1e-4) "p < 0.0001"
               else sprintf("p = %.4f", pv)

      pal <- c("#4C72B0", "#B22222")[seq_along(groups_use)]
      names(pal) <- groups_use
      y_rng <- range(d$value)
      y_pad <- diff(y_rng) * 0.18

      # Significance bar: a horizontal segment above the boxes with two short
      # ticks pointing DOWN toward the boxes. In ggplot the y axis grows
      # UPWARD, so downward ticks use yend = bar_y - tick_h (verified by pixel
      # measurement of an A/B render: -tick_h => downward, +tick_h => upward).
      if (signif_bar != "none") {
        means <- sapply(groups_use, function(g) mean(d$value[d$group == g]))
        # Anchor above BOTH the tallest data point and the taller group mean,
        # otherwise the y-limit below clips the tallest sample out of the
        # boxplot (ggplot silently drops it: "Removed 1 row ... stat_boxplot").
        bar_y     <- max(y_rng[2], unname(means)) + y_pad * bar_gap
        tick_h    <- y_pad * 0.18
        p_label_y <- bar_y + y_pad * 0.18
      } else {
        p_label_y <- y_rng[2] + y_pad
      }

      # Classic Tukey rule. Flag outliers with the SAME rule ggplot uses
      # (type-7 quantiles, fences at +- coef * IQR) so we can draw them as
      # dedicated markers at the group centre AND drop them from the jitter
      # layer - otherwise the same sample appears twice (jitter point +
      # outlier marker). With large box_coef nothing is flagged, so this is a
      # no-op for the "whiskers to min/max" mode. The rule is computed per
      # group, matching ggplot's per-group stat.
      d$is_outlier <- FALSE
      for (g in groups_use) {
        idx <- which(d$group == g)
        v <- d$value[idx]
        if (length(v) >= 2) {
          qs  <- quantile(v, c(0.25, 0.75))
          iqr <- qs[2] - qs[1]
          d$is_outlier[idx] <- v < qs[1] - box_coef * iqr | v > qs[2] + box_coef * iqr
        }
      }
      n_out <- sum(d$is_outlier)
      if (n_out > 0) {
        cat("  Outlier(s) in ", ct, ": ",
            paste(sprintf("%s=%.2f", as.character(d$group[d$is_outlier]),
                          d$value[d$is_outlier]), collapse = ", "),
            "\n", sep = "")
      }

      p_box <- ggplot2::ggplot(d, ggplot2::aes(x = group, y = value, fill = group)) +
        # Classic Tukey whiskers (coef = 1.5 by default): the whisker reaches the
        # most extreme DATA POINT inside the fence, not the fence itself. With a
        # tiny n, ggplot's type-7 quantile can put Q3 above every non-outlier
        # value, so the upper whisker can collapse to zero length and the top
        # sample is drawn as an outlier marker instead.
        # outlier.shape is set to NA and the outliers are drawn by our own
        # geom_point layer below: the outlier.* parameters were silently ignored
        # by this ggplot2 version (ggplot_build showed shape = 19 and nothing was
        # rendered), so we draw them explicitly to guarantee they appear.
        ggplot2::geom_boxplot(width = 0.5, alpha = 0.55,
                              colour = "black", linewidth = 0.4, coef = box_coef,
                              outlier.shape = NA) +
        ggplot2::geom_jitter(data = d[!d$is_outlier, , drop = FALSE],
                             width = 0.12, size = 1.8, alpha = 0.85,
                             ggplot2::aes(colour = group), shape = 16) +
        # Outliers: hollow black circles on the group centre line, sized to
        # match the jitter points (size = 1.8)
        ggplot2::geom_point(data = d[d$is_outlier, , drop = FALSE],
                            shape = 1, size = 1.8, colour = "black")
      if (signif_bar != "none") {
        # Three annotate() calls (one per segment) so every yend is a scalar;
        # annotate("segment") silently drops all but the first vector yend.
        p_box <- p_box +
          ggplot2::annotate("segment", x = 1, xend = 2,
                            y = bar_y, yend = bar_y, linewidth = 0.6) +
          ggplot2::annotate("segment", x = 1, xend = 1,
                            y = bar_y, yend = bar_y - tick_h, linewidth = 0.6) +
          ggplot2::annotate("segment", x = 2, xend = 2,
                            y = bar_y, yend = bar_y - tick_h, linewidth = 0.6)
      }
      p_box <- p_box +
        ggplot2::annotate("text", x = 1.5, y = p_label_y, label = p_lab,
                          size = 4, fontface = "bold") +
        ggplot2::scale_y_continuous(
          limits = c(min(0, y_rng[1] - y_pad * 0.1), p_label_y + y_pad * 0.18),
          expand = ggplot2::expansion(mult = c(0.01, 0.02))
        ) +
        ggplot2::scale_fill_manual(values = pal) +
        ggplot2::scale_colour_manual(values = pal) +
        ggplot2::labs(x = NULL, y = paste0(gene, " average expression")) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(
          axis.text.x  = ggplot2::element_text(size = 11, face = "bold", colour = "black"),
          axis.text.y  = ggplot2::element_text(size = 10, face = "bold", colour = "black"),
          axis.title.y = ggplot2::element_text(size = 11, face = "bold", colour = "black",
                                               margin = ggplot2::margin(r = 14)),
          axis.line    = ggplot2::element_line(colour = "black", linewidth = 0.6),
          axis.ticks   = ggplot2::element_line(colour = "black", linewidth = 0.6),
          legend.position = "none",
          panel.grid.minor = ggplot2::element_blank(),
          plot.margin  = ggplot2::margin(t = 5, r = 4, b = 4, l = 4, unit = "mm")
        )
      f_box_png <- file.path(outpath, paste0(prefix, gene, "_", safe_ct, "_boxplot.png"))
      f_box_pdf <- file.path(outpath, paste0(prefix, gene, "_", safe_ct, "_boxplot.pdf"))
      ggplot2::ggsave(f_box_png, p_box, width = box_width, height = box_height,
                      dpi = 300, bg = "white", limitsize = FALSE)
      ggplot2::ggsave(f_box_pdf, p_box, width = box_width, height = box_height,
                      limitsize = FALSE)
      cat("Boxplot : ", basename(f_box_png), "(png/pdf) saved successfully !\n", sep = "")

      med <- sapply(groups_use, function(g) median(d$value[d$group == g]))
      mn  <- sapply(groups_use, function(g) mean(d$value[d$group == g]))
      stats_rows[[length(stats_rows) + 1]] <- data.frame(
        gene = gene, celltype = ct,
        group1 = groups_use[1], group2 = groups_use[2],
        n_group1 = n_per_group[1], n_group2 = n_per_group[2],
        median_group1 = sprintf("%.3f", med[1]),
        median_group2 = sprintf("%.3f", med[2]),
        mean_group1   = sprintf("%.3f", mn[1]),
        mean_group2   = sprintf("%.3f", mn[2]),
        test = paste0("wilcoxon_one_sided_", alt),
        n_outliers = n_out,
        p_value = format(pv, scientific = TRUE, digits = 4),
        stringsAsFactors = FALSE
      )
    }

    stats_df <- do.call(rbind, stats_rows)
    stats_f  <- file.path(outpath, paste0(prefix, gene, "_boxplot_stats.csv"))
    write.csv(stats_df, stats_f, row.names = FALSE, quote = FALSE)
    cat("\nFile : ", basename(stats_f), " saved successfully !\n", sep = "")

    write_summary("Barplots     : ", length(cts_use), " celltype(s), ",
                  bar_width, "x", bar_height, " in, order = ",
                  paste(bar_levels, collapse = ","), "\n")
    write_summary("Boxplots     : ", length(cts_use), " celltype(s), ",
                  box_width, "x", box_height, " in, groups = ",
                  paste(groups_use, collapse = ","), "\n")
    write_summary("Stats CSV    : ", basename(stats_f), "\n")
    write_summary("\n---------- Wilcoxon rank-sum test (per sample) ----------\n")
    write_summary("Alternative  : ",
                  if (alternative == "auto") "less (one-sided: group2 > group1)"
                  else alternative, "\n")
    write_summary("Points are samples, not cells. Group sizes: ",
                  paste(groups_use, "=",
                        sapply(groups_use, function(g)
                          sum(box_long$group == g & box_long$celltype == cts_use[1])),
                        collapse = ", "), "\n")
    cat("\n")
    cat(capture.output(print(stats_df, row.names = FALSE)), sep = "\n")
  }
}

close(summary_con)
cat("\nFile 2 : ", prefix, "GeneDE_summary.txt saved successfully !\n", sep = "")
cat("\n============================All Done============================\n")
