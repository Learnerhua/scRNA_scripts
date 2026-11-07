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
    A script for gene set scoring of single-cell data.
  Example: scRNA_genesetScore -I test_seo_annotated.rds -A RNA -AC celltype -G geneSet.txt -P test_ -O test
  '
  )
)
#Required
parser$add_argument(
  "-I",
  "--input",
  dest = "input_rds",
  type = "character",
  required = TRUE,
  help = 'Path to an rds file for a seurat object used for scRNA data annotation. (Required)'
)
parser$add_argument(
  "-G",
  "--geneset",
  dest = "geneset",
  type = "character",
  required = TRUE,
  help = 'Path to a text file containing gene names for gene set scoring. One gene per line (Required)'
)
parser$add_argument(
  "-AC",
  "--annotation_column",
  dest = "annotation_column",
  type = "character",
  help = 'Column name in the Seurat object used for GSVA analysis. (Required)')

# Optional
parser$add_argument(
  "-A",
  "--assay",
  dest = "assay",
  type = "character",
  default = "RNA",
  help = 'Assay name to use for gene set scoring. (Default: RNA)'
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
input_rds <- opts$input_rds
geneset <- opts$geneset
annotation_column <- opts$annotation_column
assay <- opts$assay
prefix <- opts$prefix
outpath <- opts$outpath
file_n <- 0
plot_n <- 0

message("\n=======================scRNA GeneSet Score=====================\n")
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

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(GSVA,quietly = TRUE)){BiocManager::install("GSVA")}
  if(!require(AUCell,quietly = TRUE)){BiocManager::install("AUCell")}
})

# 3. Load data
cat("\nLoading Seurat object from: ", input_rds, "...\n")
seu <- readRDS(input_rds)
genelist_raw <- readLines(geneset)
genelist_raw <- trimws(genelist_raw)
genelist_unique <- unique(genelist_raw)
seu_genes <- rownames(seu)
genes_not_in <- setdiff(genelist_unique, seu_genes)
genes_in_seu <- intersect(genelist_unique, seu_genes)
cat("Total genes read in:        ", length(genelist_raw), "\n")
cat("Unique genes after dedup:   ", length(genelist_unique), "\n")
cat("Genes found in Seurat:      ", length(genes_in_seu), "\n")
cat("Genes NOT found in Seurat:  ", length(genes_not_in), "\n\n")
if (length(genes_not_in) > 0) {
  cat("Genes not found in Seurat object:\n")
  cat(paste(genes_not_in, collapse = ", "), "\n")
} else {
  cat("All genes were found in the Seurat object.\n")
}

# 4. GSVA score
cat("\nCalculating GSVA scores...\n")
Idents(seu) <- annotation_column
pb <- AggregateExpression(
  seu,
#  group.by = "celltype",
  assays   = assay,
  slot     = "counts"
)
pb_counts <- pb[[assay]]
pb_norm <- log2(t(t(pb_counts) / colSums(pb_counts)) * 1e6 + 1)
gene_sets <- list(
  MyGeneSet = genes_in_seu
)
gsvapar <- gsvaParam(
  expr = as.matrix(pb_norm),
  geneSets = gene_sets,
  kcdf = "Gaussian",   # pseudo-bulk / log-normalized 表达用 Gaussian
)
gsva_res <- gsva(
  gsvapar,
  verbose = TRUE
)
gsva_res_df <- data.frame(
  celltype = colnames(gsva_res),
  GSVA_score = as.numeric(gsva_res[1, ])
)

gsva_df_plot <- gsva_res_df %>%
  mutate(
    sign = ifelse(GSVA_score >= 0, "Positive", "Negative"),
    abs_score = abs(GSVA_score)
  ) %>%
  arrange(
    desc(sign == "Positive"),          # 正值整体在上
    ifelse(
      sign == "Positive",
      desc(abs_score),                 # 正值：绝对值从大到小
      abs_score                         # 负值：绝对值从小到大
    )
  ) %>%
  mutate(
    celltype = factor(celltype, levels = rev(unique(celltype)))
  )

offset <- 0.01
p1 <- ggplot(gsva_df_plot, aes(x = celltype, y = GSVA_score, fill = sign)) +
  geom_col(width = 0.7) +
  geom_text(
    aes(
      y = ifelse(GSVA_score >= 0, -offset, offset),
      label = celltype,
      color = sign
    ),
    hjust = ifelse(gsva_df_plot$GSVA_score >= 0, 1, 0),
    size = 6,
    fontface = "bold"
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Positive" = "#D73027",
      "Negative" = "#4575B4"
    )
  ) +
  scale_color_manual(
    values = c(
      "Positive" = "#D73027",
      "Negative" = "#4575B4"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.text.x  = element_text(color = "black", size = 12, face = "bold"),
    axis.title.x = element_text(color = "black", size = 14, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "GSVA score"
  )

ggsave(
  filename = file.path(outpath, paste0(prefix, "GSVA_score_barplot.png")), plot = p1, width = 16, height = 12, dpi = 300, bg = "white"
)
ggsave(
  filename = file.path(outpath, paste0(prefix, "GSVA_score_barplot.pdf")), plot = p1, width = 16, height = 12
)
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", paste0(prefix, "GSVA_score_barplot.png (pdf)"), " has been generated  successfully !\n", sep = "")

# 5. ssGSEA score
cat("\nCalculating ssGSEA scores...\n")
ssgseapar <- ssgseaParam(
  expr = as.matrix(pb_norm),
  geneSets = gene_sets
)
gsva_res_ssgsea <- gsva(ssgseapar, verbose = TRUE)
ssGSEA_df <- data.frame(
  celltype = colnames(gsva_res_ssgsea),
  ssGSEA_score = as.numeric(gsva_res_ssgsea[1, ])
)
ssGSEA_df <- ssGSEA_df %>%
  arrange(desc(ssGSEA_score)) %>%
  mutate(celltype = factor(celltype, levels = celltype))  # 保持排序
p2 <- ggplot(ssGSEA_df, aes(x = celltype, y = ssGSEA_score)) +
  geom_col(
    fill = "#9ecae1",   # 统一浅蓝色
    width = 0.8
  ) +
  scale_y_continuous(
    expand = c(0, 0)    # 底部紧贴 x 轴
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Cell type",
    y = "ssGSEA score"
  ) +
  theme(
    # 坐标轴标题
    axis.title.x = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20, face = "bold"),

    # 坐标轴刻度文字
    axis.text.x  = element_text(size = 14, face = "bold", angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14, face = "bold"),

    # 坐标轴线
    axis.line = element_line(size = 1, colour = "#B0B0B0"),

    # 刻度线
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1, colour = "#B0B0B0"),

    # 画布边距
    plot.margin = margin(
      t = 16,  # 上
      r = 16,  # 右
      b = 16,  # 下
      l = 16   # 左
    )
  )
ggsave(
  filename = file.path(outpath, paste0(prefix, "ssGSEA_score_barplot.png")), plot = p2, width = 12, height = 8, dpi = 300, bg = "white"
)
ggsave(
  filename = file.path(outpath, paste0(prefix, "ssGSEA_score_barplot.pdf")), plot = p2, width = 12, height = 8
)
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", paste0(prefix, "ssGSEA_score_barplot.png (pdf)"), " has been generated  successfully !\n", sep = "")

# 6. Save score results
cat("\nSaving score results...\n")
score_merged_df <- ssGSEA_df %>%
  left_join(gsva_res_df, by = "celltype")
write.csv(
  score_merged_df,
  file = file.path(outpath, paste0(prefix, "GSVA_ssGSEA_score.csv")),
  quote = FALSE, 
  row.names = FALSE
)
file_n <- file_n + 1
cat("File ", file_n, " : ", paste0(prefix, "GSVA_ssGSEA_score.csv"), " saved successfully !\n", sep = "")


# 7. AddModuleScore
cat("\nCalculating AddModuleScore...\n")
seu <- AddModuleScore(
  object = seu,
  features = list(genes_in_seu),
  assay = assay,
  name = "AddModuleScore"
)
p3 <- FeaturePlot(
  object = seu,
  features = "AddModuleScore1",
  reduction = "umap",
  cols = c("lightgrey", "red")
)+
  theme(plot.title = element_blank())
ggsave(
  filename = file.path(outpath, paste0(prefix, "AddModuleScore_UMAP.png")), plot = p3, width = 8, height = 6, dpi = 300, bg = "white"
)
ggsave(
  filename = file.path(outpath, paste0(prefix, "AddModuleScore_UMAP.pdf")), plot = p3, width = 8, height = 6
  )
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", paste0(prefix, "AddModuleScore_UMAP.png (pdf)"), " has been generated  successfully !\n", sep = "")

# 8. AUCell
cat("\nCalculating AUCell scores...\n")
seu <- NormalizeData(
  object = seu,
  assay = assay,
  normalization.method = "LogNormalize",
  scale.factor = 1e6
)
expr_mat <- GetAssayData(
  object = seu,
  assay = assay,
  layer = "data"
)
cells_rankings <- AUCell_buildRankings(
  expr_mat,
  plotStats = FALSE
)
gene_sets <- list(
  GeneSet = genes_in_seu
)
cells_auc <- AUCell_calcAUC(
  gene_sets,
  cells_rankings
)
auc_matrix <- as.matrix(getAUC(cells_auc))
seu$AUCell_score <- auc_matrix["GeneSet", colnames(seu)]
p4 <- FeaturePlot(
  object = seu,
  features = "AUCell_score",
  reduction = "umap",
  cols = c("lightgrey", "red")
)+
  theme(plot.title = element_blank())
ggsave(
  filename = file.path(outpath, paste0(prefix, "AUCell_score_UMAP.png")), plot = p4, width = 8, height = 6, dpi = 300, bg = "white"
)
ggsave(
  filename = file.path(outpath, paste0(prefix, "AUCell_score_UMAP.pdf")), plot = p4, width = 8, height = 6
)
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", paste0(prefix, "AUCell_score_UMAP.png (pdf)"), " has been generated  successfully !\n", sep = "")

# 9. Save Seurat object
cat("\nSaving Seurat object with gene set scores...\n")
saveRDS(
  seu,
  file = file.path(outpath, paste0(prefix, "seu_genesetScore.rds"))
)
file_n <- file_n + 1
cat("File ", file_n, " : ", paste0(prefix, "seu_genesetScore.rds"), " saved successfully !\n", sep = "")

cat("\n============================All Done============================\n")
