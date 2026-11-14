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
    "
    A script for gene set scoring of single-cell data. Example: scRNA_hdWGCNA -I seu_annotated.rds -AC celltype -SC orig.ident -TC 'PT(S2)' -A RNA -SP 3 -T 8 -O test -P test_
  "
  )
)

#Required
parser$add_argument(
    "-I",
    "--input",
    dest = "input_rds",
    type = "character",
    required = TRUE,
    help = 'Path to an rds file for a seurat object used for scRNA hdWGCNA analysis. (Required)'
)
parser$add_argument(
    "-AC",
    "--annotation_column",
    dest = "annotation_column",
    type = "character",
    required = TRUE,
    help = 'Column name containing annotation information in the Seurat object used for hdWGCNA analysis. (Required)'
)
parser$add_argument(
    "-SC",
    "--sample_column",
    dest = "sample_column",
    type = "character",
    required = TRUE,
    help = 'Column name containing sample information in the Seurat object used for GSVA analysis. (Required)'
)
parser$add_argument(
    "-TC",
    "--target_cell",
    dest = "target_cell",
    type = "character",
    required = TRUE,
    help = 'The target cell type for hdWGCNA analysis. (Required)'
)
parser$add_argument(
    "-SP",
    "--soft_power",
    dest = "soft_power",
    type = "integer",
    default = 3,
    help = "Soft thresholding power for WGCNA.  You can check the plot 'softPowers.png' for reference (Default: 3)"
)
parser$add_argument(
    "-A",
    "--assay",
    dest = "assay",
    type = "character",
    default = "RNA",
    help = 'Assay name in the Seurat object to be used for hdWGCNA analysis. Default is "RNA".'
)
parser$add_argument(
    "-T",
    "--threads",
    dest = "threads",
    type = "integer",
    default = 8,
    help = "Number of threads to use for hdWGCNA analysis. (Default: 8)"
)

#Optional
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

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_rds <- opts$input_rds
annotation_column <- opts$annotation_column
sample_column <- opts$sample_column
target_cell <- opts$target_cell
soft_power <- opts$soft_power
threads <- opts$threads
assay <- opts$assay
prefix <- opts$prefix
outpath <- opts$outpath
file_n <- 0
plot_n <- 0

message("\n======================scRNA hdWGCNA Analysis====================\n")
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
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(hdWGCNA,quietly = TRUE)){devtools::install_github('smorabit/hdWGCNA', ref='dev')}
  if(!require(UCell,quietly = TRUE)){install.packages("UCell")}
  if(!require(ggradar,quietly = TRUE)){install.packages("ggradar")}
})
allowWGCNAThreads(nThreads = threads)

# 3.Load data
cat("\nLoading data from: ", input_rds, "...\n")
seu <- readRDS(input_rds)
seu <- NormalizeData(seu, normalization.method = "LogNormalize")
# 3.1 定义 WGCNA 分析框架
cat("\nSetting up Seurat object for WGCNA analysis...\n")
seu_1 <- SetupForWGCNA(
  seu,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = paste0(target_cell, "_hdWGCNA")
)
# 3.2 构建 metacells（核心降噪步骤）
cat("\nConstructing metacells...\n")
seu_2 <- MetacellsByGroups(
  seurat_obj = seu_1,
  group.by = c(annotation_column, sample_column),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = annotation_column,
  min_cells = 100
)
cat("\nNormalizing and running dimensionality reduction on metacells...\n")
seu_2 <- NormalizeMetacells(seu_2)
seu_2 <- ScaleMetacells(seu_2, features = VariableFeatures(seu_2))
seu_2 <- RunPCAMetacells(seu_2, features = VariableFeatures(seu_2))
seu_2 <- RunHarmonyMetacells(seu_2, group.by.vars = sample_column)
seu_2 <- RunUMAPMetacells(seu_2, reduction = 'harmony', dims = 1:30)
p1 <- DimPlotMetacells(seu_2, group.by = annotation_column) +
    theme(plot.title = element_blank())
ann_col <- gsub(" ", "_", annotation_column)
ggsave(file.path(outpath, paste0(prefix, "metacell_", ann_col, "umap.png")), p1, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "metacell_", ann_col, "umap.pdf")), p1, width = 10, height = 8)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "metacell_", ann_col, "umap.png (pdf) has been generated successfully !\n", sep = "")
p2 <- DimPlotMetacells(seu_2, group.by = sample_column) + 
    theme(plot.title = element_blank())
sam_col <- gsub(" ", "_", sample_column)
ggsave(file.path(outpath, paste0(prefix, "metacell_", sam_col, "umap.png")), p2, width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "metacell_", sam_col, "umap.pdf")), p2, width = 10, height = 8, dpi = 300, bg = "white")
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "metacell_", sam_col, "umap.png (pdf) has been generated successfully !\n", sep = "")
# 3.3 构建 WGCNA 专用表达矩阵
cat("\nConstructing WGCNA expression matrix...\n")
seu_3 <- SetDatExpr(
  seu_2,
  group_name = target_cell,
  group.by = annotation_column,
  assay = assay,
  layer = "data"
)
# 3.4 软阈值扫描
cat("\nTesting soft-thresholding powers for WGCNA...\n")
seu_4 <- TestSoftPowers(
  seu_3,
  powers = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  networkType = "unsigned"
)
plot_list1 <- PlotSoftPowers(seu_4)
p3 <- wrap_plots(plot_list1, ncol=2)
ggsave(file.path(outpath, paste0(prefix, "softPowers.png")), p3, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "softPowers.pdf")), p3, width = 12, height = 8)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "softPowers.png (pdf) has been generated successfully !\n", sep = "")
# 3.5 构建网络
cat("\nConstructing WGCNA network...\n")
seu_5 <- ConstructNetwork(
  seu_4,
  soft_power   = soft_power,
  networkType  = "unsigned",
  tom_name     = paste0(target_cell, "_TOM"),
  tom_outdir   = file.path(outpath, "TOM"),
  overwrite_tom = TRUE
)
png(file.path(outpath, paste0(prefix, "Dendrogram.png")), width = 12, height = 8, res = 300, units = "in")
PlotDendrogram(seu_5, main='HdWGCNA Dendrogram')
dev.off()
pdf(file.path(outpath, paste0(prefix, "Dendrogram.pdf")), width = 12, height = 8)
PlotDendrogram(seu_5, main='HdWGCNA Dendrogram')
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "Dendrogram.png (pdf) has been generated successfully !\n", sep = "" )
# 3.6 获取模块特征值(模块在细胞中的打分)
cat("\nCalculating module eigengenes...\n")
seu_6 <- ModuleEigengenes(
  seu_5, 
  group.by = sample_column 
)
# 3.7 计算模块连接性（kMEs，基因在模块中的打分）
cat("\nCalculating module connectivity (kMEs)...\n")
seu_7 <- ModuleConnectivity(
  seu_6,
  group.by = annotation_column, 
  group_name = target_cell
)
seu_7 <- ResetModuleNames(
  seu_7,
  new_name = paste0(target_cell, "_NEW_")
)
p4 <- PlotKMEs(seu_7, ncol=4)
ggsave(file.path(outpath, paste0(prefix, "plotKMEs.png")), p4, width = 12, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "plotKMEs.pdf")), p4, width = 12, height = 8)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "plotKMEs.png (pdf) has been generated successfully !\n", sep = "")
# 3.8 获取模块内部的基因信息（基因和模块之间的关系）
cat("\nExtracting module gene information...\n")
modules <- GetModules(seu_7) %>% 
  subset(module != 'grey')
write.csv(modules, file.path(outpath, paste0(prefix, "gene_modules.csv")), quote = FALSE, row.names = TRUE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "gene_modules.csv saved successfully !\n", sep = "" )
hub_df <- GetHubGenes(seu_7)
write.csv(hub_df, file.path(outpath, paste0(prefix, "top10_hubgenes.csv")), quote = FALSE, row.names = TRUE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "top10_hubgenes.csv saved successfully !\n", sep = "" )
# 3.9 模块前25个基因（按kME）在seu对象中的得分情况
cat("\nCalculating module expression scores...\n")
seu_8 <- ModuleExprScore(
  seu_7,
  n_genes = 25,
  method='UCell'
)
MEs <- GetMEs(seu_8, harmonized = FALSE)
write.csv(MEs, file.path(outpath, paste0(prefix, "MEs.csv")), quote = FALSE, row.names = TRUE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "MEs.csv saved successfully !\n", sep = "" )
hMEs <- GetMEs(seu_8, harmonized = TRUE)
write.csv(hMEs, file.path(outpath, paste0(prefix, "hMEs.csv")), quote = FALSE, row.names = TRUE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "hMEs.csv saved successfully !\n", sep = "" )

# 4 可视化
cat("\nGenerating visualizations...\n")
# 4.1 特征得分
plot_list2 <- ModuleFeaturePlot(
  seu_8,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
p5 <- wrap_plots(plot_list2)
ggsave(file.path(outpath, paste0(prefix, "cell_module.png")), p5, width = 18, height = 16, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "cell_module.pdf")), p5, width = 18, height = 16)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "cell_module.png (pdf) has been generated successfully !\n", sep = "")
plot_list3 <- ModuleFeaturePlot(
  seu_8,
  features='scores', # plot the hub gene scores
  order='shuffle', # order so cells are shuffled
  ucell = TRUE # depending on Seurat vs UCell for gene scoring
)
p6 <- wrap_plots(plot_list3)
ggsave(file.path(outpath, paste0(prefix, "cell_hubGenes.png")), p6, width = 18, height = 16, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "cell_hubGenes.pdf")), p6, width = 18, height = 16)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "cell_module.png (pdf) has been generated successfully !\n", sep = "")
target_cells <- rownames(
  seu_8@meta.data[seu_8@meta.data[[annotation_column]] == target_cell, , drop = FALSE]
)
p7 <- ModuleRadarPlot(
  seu_8,
  group.by = sample_column,
  barcodes = target_cells, 
  axis.label.size = 4,
  grid.label.size = 4
)
ggsave(file.path(outpath, paste0(prefix, target_cell, "_radar.png")), p7,  width = 18, height = 16, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, target_cell, "_radar.pdf")), p7,  width = 18, height = 16)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", target_cell, "_radar.png (pdf) has been generated successfully !\n", sep = "")
# 4.2模块之间的关系
png(file.path(outpath, paste0(prefix, "module_correlation.png")), width = 8, height = 8, res = 300, units = "in")
ModuleCorrelogram(seu_8)
dev.off()
pdf(file.path(outpath, paste0(prefix, "module_correlation.pdf")), width = 8, height = 8)
ModuleCorrelogram(seu_8)
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "module_correlation.png (pdf) has been generated successfully !\n", sep = "")
# 4.3 细胞类型与模块的相关性气泡图
mods <- levels(modules$module); mods <- mods[mods != 'grey']
seu_8@meta.data <- cbind(seu_8@meta.data, hMEs)
p8 <- DotPlot(seu_8, features=mods, group.by = annotation_column)
p8 <- p8 +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')
ggsave(file.path(outpath, paste0(prefix, "module_celltype_dotplot.png")), p8, width = 18, height = 16, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "module_celltype_dotplot.pdf")), p8, width = 18, height = 16)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", prefix, "module_celltype_dotplot.png (pdf) has been generated successfully !\n", sep = "")

# 5. Save data
cat("\nSaving Seurat object with hdWGCNA results...\n")
saveRDS(seu_8, file.path(outpath, paste0(prefix, "seu_hdWGCNA.rds")))
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "seu_hdWGCNA.rds saved successfully !\n", sep = "")
message("\n=========================All Done !============================\n")







