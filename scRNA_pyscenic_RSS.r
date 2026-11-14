#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
# =============================================================================
# Regulon specificity (RSS) analysis on a pySCENIC AUCell loom.
#
# Read the loom written by `pyscenic aucell`, match its cells against a Seurat
# object, compute the Regulon Specificity Score per (regulon x cell type), draw
# the summary figures and write the regulon AUC matrix back into the object.
#
# Example:
#   scRNA_pyscenic_RSS.r -I seu_annotated.rds -L 01_pyscenic/01.SCENIC.loom \
#       -AC celltype -TF 'Jun(+),Foxp2(+)' -P 02. -O 02_pyscenic_RSS
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
    Compute the Regulon Specificity Score (RSS) per regulon and cell type from a
    pySCENIC AUCell loom, draw specificity figures, and write the regulon AUC
    matrix back into the Seurat object.
  Example: scRNA_pyscenic_RSS.r -I seu_annotated.rds -L 01.SCENIC.loom -AC celltype -TF "Jun(+),Foxp2(+)" -P 02. -O 02_RSS
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
  help = 'Path to the Seurat rds used to build the loom (annotated). (Required)'
)
parser$add_argument(
  "-L",
  "--input_loom",
  dest = "input_loom",
  type = "character",
  required = TRUE,
  help = 'Loom file written by `pyscenic aucell`. (Required)'
)
parser$add_argument(
  "-AC",
  "--annotation_column",
  dest = "annotation_column",
  type = "character",
  required = TRUE,
  help = 'Metadata column holding the cell type used for RSS. (Required)'
)

# Optional
parser$add_argument(
  "-TF",
  "--target_tfs",
  dest = "target_tfs",
  type = "character",
  default = "",
  help = 'Comma-separated regulon names to plot individually, e.g. "Jun(+),Foxp2(+)". Empty means pick the top -NT regulons per cell type. (Optional, default: empty)'
)
parser$add_argument(
  "-NT",
  "--n_top",
  dest = "n_top",
  type = "integer",
  default = 2,
  help = 'Number of top-RSS regulons taken per cell type for the dotplot and the summary heatmap. (Optional, default: 2)'
)
parser$add_argument(
  "-OS",
  "--one_set",
  dest = "one_set",
  type = "character",
  default = "",
  help = 'Comma-separated cell types for the per-cell-type RSS panel (2x2 grid). Empty means the top 4 cell types by max RSS. (Optional, default: empty)'
)
parser$add_argument(
  "-ZW",
  "--z_threshold",
  dest = "z_threshold",
  type = "double",
  default = 3,
  help = 'zThreshold passed to plotRSS. (Optional, default: 3)'
)
parser$add_argument(
  "-TH",
  "--rss_threshold",
  dest = "rss_threshold",
  type = "double",
  default = 0.1,
  help = 'thr passed to plotRSS; RSS below this is left blank. (Optional, default: 0.1)'
)
parser$add_argument(
  "-W",
  "--width",
  dest = "width",
  type = "double",
  default = 10,
  help = 'Width (inches) for the RSS summary figure only; the multi-panel regulon figures size themselves. (Optional, default: 10)'
)
parser$add_argument(
  "-H",
  "--height",
  dest = "height",
  type = "double",
  default = 12,
  help = 'Height (inches) for the RSS summary figure only; the multi-panel regulon figures size themselves. (Optional, default: 12)'
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
input_loom <- opts$input_loom
annotation_column <- opts$annotation_column
target_tfs <- opts$target_tfs
n_top <- opts$n_top
one_set <- opts$one_set
z_threshold <- opts$z_threshold
rss_threshold <- opts$rss_threshold
width <- opts$width
height <- opts$height
prefix <- opts$prefix
outpath <- opts$outpath

if (!file.exists(input_rds)) {
  stop("Input rds does not exist: ", input_rds)
}
if (!file.exists(input_loom)) {
  stop("Input loom does not exist: ", input_loom)
}
if (n_top < 1) {
  stop("-NT/--n_top must be >= 1.")
}

message("\n====================scRNA pySCENIC RSS=====================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2026 OYJH. All Rights Reserved.")
message("\n===========================================================\n")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if (!require(Seurat, quietly = TRUE)) install.packages("Seurat")
  if (!require(SCENIC, quietly = TRUE)) install.packages("SCENIC")
  if (!require(SCopeLoomR, quietly = TRUE)) install.packages("SCopeLoomR")
  if (!require(gridExtra, quietly = TRUE)) install.packages("gridExtra")
  if (!require(tidyverse, quietly = TRUE)) install.packages("tidyverse")
})

file_n <- 0
plot_n <- 0

save_gg <- function(p, stem, w, h) {
  ggsave(file.path(outpath, paste0(prefix, stem, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, stem, ".pdf")), p, width = w, height = h)
  plot_n <<- plot_n + 1
  cat("\nPlot ", plot_n, " : ", prefix, stem, ".png(pdf) has been generated successfully !\n", sep = "")
}

# Lay n panels out on a near-square grid so the figure stays roughly square
# instead of becoming a long strip when many regulons are plotted.
grid_dims <- function(n, panel) {
  ncol <- ceiling(sqrt(n))
  nrow <- ceiling(n / ncol)
  list(ncol = ncol, nrow = nrow, width = ncol * panel, height = nrow * panel)
}

# 3. Load loom and extract regulons / AUC
cat(">>> Reading loom\n")
loom <- open_loom(input_loom)
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
close_loom(loom)
cat(">>> Regulons: ", length(regulons), " | AUC matrix: ", nrow(regulonAUC), " regulons x ",
    ncol(regulonAUC), " cells\n", sep = "")

# 4. Match cells against the Seurat object
cat(">>> Loading Seurat object\n")
seu <- readRDS(input_rds)
if (!annotation_column %in% colnames(seu@meta.data)) {
  stop("Column not found in metadata: ", annotation_column,
       "; available: ", paste(colnames(seu@meta.data), collapse = ", "))
}
common_cells <- intersect(colnames(seu), colnames(regulonAUC))
cat(">>> Cells shared by object and loom: ", length(common_cells), "\n", sep = "")
if (length(common_cells) < 3) {
  stop("Fewer than 3 shared cells between the Seurat object and the loom; ",
       "check that -I is the object the loom was built from.")
}

seu_subset <- subset(seu, cells = common_cells)
sub_regulonAUC <- regulonAUC[, match(colnames(seu_subset), colnames(regulonAUC))]
if (!identical(colnames(sub_regulonAUC), colnames(seu_subset))) {
  stop("Cell order mismatch between the AUC matrix and the Seurat subset.")
}

cell_annotation <- as.character(seu_subset@meta.data[[annotation_column]])
if (any(is.na(cell_annotation))) {
  stop("NA values found in the annotation column: ", annotation_column)
}
cat(">>> Cell types: ", length(unique(cell_annotation)), "\n", sep = "")
print(table(cell_annotation))

# 5. Compute RSS
cat(">>> Computing RSS\n")
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = cell_annotation)
rss_file <- paste0(prefix, "rss.csv")
write.csv(rss, file.path(outpath, rss_file), row.names = TRUE, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", rss_file, " saved successfully !\n", sep = "")

# 6. Pick regulons and cell types to display
if (nzchar(target_tfs)) {
  tfs_use <- trimws(strsplit(target_tfs, ",")[[1]])
  missing_tfs <- setdiff(tfs_use, rownames(rss))
  if (length(missing_tfs) > 0) {
    cat("!!! Regulon(s) not found and skipped: ", paste(missing_tfs, collapse = ", "), "\n", sep = "")
    tfs_use <- intersect(tfs_use, rownames(rss))
  }
  if (length(tfs_use) == 0) {
    stop("None of the -TF/--target_tfs regulons exist in the loom.")
  }
} else {
  rss_long <- as.data.frame(rss) %>%
    rownames_to_column(var = "TF") %>%
    pivot_longer(cols = -TF, names_to = "CellType", values_to = "Score")
  top_tfs_df <- rss_long %>%
    group_by(CellType) %>%
    arrange(desc(Score)) %>%
    slice_head(n = n_top) %>%
    ungroup()
  tfs_use <- unique(top_tfs_df$TF)
  cat(">>> Auto-selected ", length(tfs_use), " regulon(s) across cell types\n", sep = "")
}

if (nzchar(one_set)) {
  sets_use <- trimws(strsplit(one_set, ",")[[1]])
  sets_use <- intersect(sets_use, colnames(rss))
  if (length(sets_use) == 0) {
    stop("None of the -OS/--one_set cell types exist in the RSS matrix.")
  }
} else {
  cell_max <- apply(rss, 2, max)
  sets_use <- names(sort(cell_max, decreasing = TRUE))[seq_len(min(4, length(cell_max)))]
}

# 7. RSS summary heatmap
cat(">>> Drawing RSS summary heatmap\n")
p1 <- plotRSS(
  rss,
  zThreshold = z_threshold,
  cluster_columns = FALSE,
  order_rows = TRUE,
  thr = rss_threshold,
  varName = "cellType",
  col.low = "#cfd9ed",
  col.mid = "#66CC66",
  col.high = "#FFCC33"
)
save_gg(p1$plot, "RSS_plot", width, height)

# 8. Per-cell-type RSS panels (2x2)
cat(">>> Drawing per-cell-type RSS panels\n")
one_set_plots <- lapply(sets_use, function(ct) plotRSS_oneSet(rss, ct))
n_panel <- length(one_set_plots)
n_col <- if (n_panel <= 1) 1 else 2
grid_obj <- do.call(grid.arrange, c(one_set_plots, ncol = n_col, nrow = ceiling(n_panel / n_col)))
ggsave(file.path(outpath, paste0(prefix, "RSS_oneSet.png")), grid_obj, width = 12, height = 12, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "RSS_oneSet.pdf")), grid_obj, width = 12, height = 12)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "RSS_oneSet.png(pdf) has been generated successfully !\n", sep = "")

# 9. Regulon activity on the object
cat(">>> Writing regulon AUC into the Seurat object\n")
auc_mat <- t(getAUC(sub_regulonAUC))
dup_cols <- intersect(colnames(auc_mat), colnames(seu_subset@meta.data))
if (length(dup_cols) > 0) {
  cat("!!! Metadata already holds column(s) with regulon names, overwriting: ",
      paste(dup_cols, collapse = ", "), "\n", sep = "")
  seu_subset@meta.data <- seu_subset@meta.data[, setdiff(colnames(seu_subset@meta.data), dup_cols), drop = FALSE]
}
seu_subset@meta.data <- cbind(seu_subset@meta.data, auc_mat[rownames(seu_subset@meta.data), , drop = FALSE])

Idents(seu_subset) <- annotation_column
gd2 <- grid_dims(length(tfs_use), 4)
p2 <- VlnPlot(seu_subset, features = tfs_use, pt.size = 0, ncol = gd2$ncol)
save_gg(p2, "regulon_violin", gd2$width, gd2$height)

if (length(Reductions(seu_subset)) > 0) {
  gd3 <- grid_dims(length(tfs_use), 4)
  p3 <- FeaturePlot(seu_subset, features = tfs_use, pt.size = 0, ncol = gd3$ncol)
  save_gg(p3, "regulon_feature", gd3$width, gd3$height)
} else {
  cat("!!! No dimensionality reduction found; skipping the feature plot.\n")
}

# Seurat's DotPlot marks a feature as detected when AUC > 0, but AUCell scores are
# positive in almost every cell, so that fraction saturates near 100% and the dot
# size stops carrying information (measured: median 0.99, half the regulons at 1.0).
# Use the fraction of cells above this regulon's global 75th percentile instead, and
# colour by the mean AUC so the two channels answer different questions:
# colour = how active, size = how consistently active.
auc_all <- getAUC(sub_regulonAUC)
q75 <- apply(auc_all, 1, quantile, probs = 0.75)
ct_vec <- as.character(seu_subset@meta.data[[annotation_column]])
ct_levels <- sort(unique(ct_vec))
ct_factor <- factor(ct_vec, levels = ct_levels)

mean_mat <- sapply(tfs_use, function(r) as.numeric(tapply(auc_all[r, ], ct_factor, mean)))
frac_mat <- sapply(tfs_use, function(r) as.numeric(tapply(auc_all[r, ] > q75[r], ct_factor, mean)))
mean_mat <- matrix(mean_mat, nrow = length(ct_levels), dimnames = list(ct_levels, tfs_use))
frac_mat <- matrix(frac_mat, nrow = length(ct_levels), dimnames = list(ct_levels, tfs_use))

dot_df <- data.frame(
  celltype = rep(rownames(mean_mat), times = ncol(mean_mat)),
  regulon  = rep(colnames(mean_mat), each = nrow(mean_mat)),
  mean_auc = as.vector(mean_mat),
  frac     = as.vector(frac_mat)
)
p4 <- ggplot(dot_df, aes(x = regulon, y = celltype)) +
  geom_point(aes(size = frac, colour = mean_auc)) +
  scale_size_continuous(range = c(0, 6), limits = c(0, 1),
                        breaks = c(0, 0.25, 0.5, 0.75, 1),
                        name = "Fraction of cells\nAUC > global Q75") +
  scale_colour_gradientn(colours = c("#cfd9ed", "#66CC66", "#FFCC33"),
                         name = "Mean AUC") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title = "Regulon activity by cell type")
save_gg(p4, "regulon_dotplot", max(8, 4 + 0.35 * length(tfs_use)), max(5, 0.35 * length(ct_levels) + 3))

# 10. Save the object carrying the regulon activity
rds_file <- paste0(prefix, "seu_SCENIC.rds")
saveRDS(seu_subset, file.path(outpath, rds_file))
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", rds_file, " saved successfully !\n", sep = "")
cat(">>> Regulon activity columns added to metadata: ", ncol(auc_mat), "\n", sep = "")

cat("\n>>> Done. Results in ", outpath, "\n", sep = "")
message("\n=========================All Done !============================\n")
