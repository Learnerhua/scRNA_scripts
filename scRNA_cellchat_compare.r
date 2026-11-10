#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
rm(list = ls())
options(timeout = 1000)
options(future.globals.maxSize = 100 * 1024^3)

# 1. Parameter settings
# ---- Define argument parser ----
suppressPackageStartupMessages({
  if (!require(argparse, quietly = TRUE)) {
    install.packages("argparse")
    library(argparse)
  }
})

parser <- ArgumentParser(description = paste0(
  "
A script for cell communication comparison analysis of two single cell transcriptomic data.
Example: scRNA_cellchat_compare -IC cellchat_control.rds -IT cellchat_treatment.rds -GN control,treatment -P CT_vs_TR -O /home/user/data/output
"
))

# ---- Required arguments ----
parser$add_argument("-IC", "--input_C",
                    help = "Path to an RDS file containing a cellchat object for single cell transcriptomics data from control group. (required)",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-IT", "--input_T",
                    help = "Path to an RDS file containing a cellchat object for single cell transcriptomics data from treatment group. (required)",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-GN", "--group_names",
                    help = "Names of the two groups (control and treatment), separated by comma. (e.g., control,treatment) (required)",
                    type = "character",
                    required = TRUE
)
# ---- Optional arguments ----
parser$add_argument("-TI", "--target_idents",
                    help = "Target cell identities to include in the comparison, separated by comma. If not specified, all cell identities will be used. (optional)",
                    type = "character"
)
parser$add_argument("-P", "--prefix",
                    help = "Prefix for output files (optional, default: none)",
                    type = "character",
                    default = "")
parser$add_argument("-O", "--outpath",
                    help = "Output directory (absolute path required). (optional, default: current working directory)",
                    type = "character",
                    default = getwd())

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_C <- opts$input_C
input_T <- opts$input_T
group_names <- unlist(strsplit(opts$group_names, split = ","))
if (length(group_names) != 2) {
  stop("Please provide exactly two group names separated by a comma.")
}
target_idents <- opts$target_idents
if (!is.null(target_idents)) {
  target_idents <- unlist(strsplit(target_idents, split = ","))
}
prefix <- opts$prefix
outpath <- opts$outpath
file_n <- 0
plot_n <- 0

message("\n=================scRNA Cell Communication Anaysis===============\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}
setwd(outpath)

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(CellChat,quietly = TRUE)){install.packages("CellChat")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(ComplexHeatmap,quietly = TRUE)){install.packages("ComplexHeatmap")}
})

# 3. merge CellChat objects
if (!file.exists(file.path(outpath, paste0(prefix, "merged_cellchat.rds")))) {
  cellchat_C <- readRDS(input_C)
  cat("\n>>> Control CellChat object loaded. Idents information:\n")
  print(table(cellchat_C@idents))
  idents_C <- names(table(cellchat_C@idents))
  cellchat_T <- readRDS(input_T)
  cat("\n>>> Treatment CellChat object loaded. Idents information:\n")
  print(table(cellchat_T@idents))
  idents_T <- names(table(cellchat_T@idents))
  if(!identical(sort(idents_C), sort(idents_T))){
    stop("The cell identities in the two CellChat objects are not identical. Please check!")
  }
  cat("\n>>> Merging CellChat objects...\n")
  if(!is.null(target_idents)){
    object.list <- list(
    cellchat_1 = subsetCellChat(cellchat_C, idents.use = target_idents),
    cellchat_2 = subsetCellChat(cellchat_T, idents.use = target_idents))
  } else {
    object.list <- list(
    cellchat_1 = cellchat_C,
    cellchat_2 = cellchat_T
    )
  }
  cellchat_merged <- mergeCellChat(object.list, add.names = group_names)
  cat(">>> Merged CellChat object created. Idents information:\n")
  print(table(cellchat_merged@idents$joint))
  saveRDS(cellchat_merged, file = file.path(outpath, paste0(prefix, "merged_cellchat.rds")))
  file_n <- file_n + 1
  cat("\nFile ", file_n, " : ", paste0(prefix,"merged_cellchat.rds"), " saved successfully !\n", sep = "")
} else {
  cellchat_merged <- readRDS(file.path(outpath, paste0(prefix, "merged_cellchat.rds")))
  cat("\n>>> Merged CellChat object loaded from existing file. Idents information:\n")
  table(cellchat_merged@idents$joint)
  # object.list is needed downstream (getMaxWeight, group sizes, per-group
  # heatmaps), so rebuild it from the two input objects.
  cellchat_C <- readRDS(input_C)
  cellchat_T <- readRDS(input_T)
  if(!is.null(target_idents)){
    object.list <- list(
      cellchat_1 = subsetCellChat(cellchat_C, idents.use = target_idents),
      cellchat_2 = subsetCellChat(cellchat_T, idents.use = target_idents))
  } else {
    object.list <- list(
      cellchat_1 = cellchat_C,
      cellchat_2 = cellchat_T)
  }
}

# 4. Visualization
## 4.1 Barplot
p1 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2))
p2 <- compareInteractions(cellchat_merged, show.legend = F, group = c(1,2), measure = "weight")
plot1 <- p1 + p2
ggsave(filename = file.path(outpath, paste0(prefix, "comparison_barplot.png")), plot = plot1, width = 8, height = 8)
ggsave(filename = file.path(outpath, paste0(prefix, "comparison_barplot.pdf")), plot = plot1, width = 8, height = 8)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"comparison_barplot.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.2 Circle plot
circle_plot <- function() {
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, label.edge = F)
  netVisual_diffInteraction(cellchat_merged, weight.scale = T, measure = "weight", label.edge = F)
}
png(filename = file.path(outpath, paste0(prefix, "comparison_circleplot.png")), res = 300, width = 12, height = 6, units = "in")
circle_plot()
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "comparison_circleplot.pdf")), width = 12, height = 6)
circle_plot()
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"comparison_circleplot.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.3 Heatmap plot
heatmap_weight_plot <- function() {
  par(mfrow = c(1,1), xpd=TRUE)
  netVisual_heatmap(cellchat_merged, measure = "weight")
}
png(filename = file.path(outpath, paste0(prefix, "weight_comparison_heatmap.png")), res = 300, width = 8, height = 8, units = "in")
heatmap_weight_plot()
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "weight_comparison_heatmap.pdf")), width = 8, height = 8)
heatmap_weight_plot()
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"weight_comparison_heatmap.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.4 Overall circleplot
weight.max <- getMaxWeight(object.list, attribute = c("idents", "weight"))
count.max <- getMaxWeight(object.list, attribute = c("idents", "count"))
groupsize_C <- as.numeric(table(object.list[[1]]@idents))
groupsize_T <- as.numeric(table(object.list[[2]]@idents))

overall_ciecleplot <- function(){
par(mfrow = c(2,2), xpd=TRUE)
netVisual_circle(object.list[[1]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = count.max[2], edge.width.max = 12, vertex.weight = groupsize_C, title.name = paste0("Number of interactions - ", names(object.list)[1]))
netVisual_circle(object.list[[2]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = count.max[2], edge.width.max = 12, vertex.weight = groupsize_T,title.name = paste0("Number of interactions - ", names(object.list)[2]))
netVisual_circle(object.list[[1]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, vertex.weight = groupsize_C, title.name = paste0("Strength of interactions - ", names(object.list)[1]))
netVisual_circle(object.list[[2]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, vertex.weight = groupsize_T, title.name = paste0("Strength of interactions - ", names(object.list)[2]))
}
png(filename = file.path(outpath, paste0(prefix, "overall_circleplot.png")), res = 300, width = 16, height = 16, units = "in")
overall_ciecleplot()
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "overall_circleplot.pdf")), width = 12, height = 12)
overall_ciecleplot()
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"overall_circleplot.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.5 2D space (cellgroup comparison)
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
 gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = group_names[i], weight.MinMax = weight.MinMax)
}
plot2 <- wrap_plots(gg, ncol = 2) 
ggsave(filename = file.path(outpath, paste0(prefix, "role_2d.png")), plot = plot2, width = 12, height = 6)
ggsave(filename = file.path(outpath, paste0(prefix, "role_2d.pdf")), plot = plot2, width = 12, height = 6)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"role_2d.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.6 pathway comparison
cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")
cellchat_merged <- netEmbedding(cellchat_merged, type = "functional")
cellchat_merged <- netClustering(cellchat_merged, type = "functional")
png(filename = file.path(outpath, paste0(prefix, "pathway_functionalClustering.png")), res = 300, width = 8, height = 8, units = "in")
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "pathway_functionalClustering.pdf")), width = 8, height = 8)
netVisual_embeddingPairwise(cellchat_merged, type = "functional", label.size = 3.5)
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"pathway_functionalClustering.png (pdf)"), " has been generated successfully!\n", sep = "")

png(filename = file.path(outpath, paste0(prefix, "pathway_functionalSimilarity.png")), res = 300, width = 8, height = 8, units = "in")
rankSimilarity(cellchat_merged, type = "functional")
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "pathway_functionalSimilarity.pdf")), width = 8, height = 8)
rankSimilarity(cellchat_merged, type = "functional")
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"pathway_functionalSimilarity.png (pdf)"), " has been generated successfully!\n", sep = "")

p3 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",stacked = T, do.stat = TRUE, cutoff.pvalue = 0.05)
p4 <- rankNet(cellchat_merged, mode = "comparison", measure = "weight",stacked = F, do.stat = TRUE, cutoff.pvalue = 0.05)
plot3 <- p3 + p4
ggsave(filename = file.path(outpath, paste0(prefix, "weight_DEG_pathways.png")), plot = plot3, width = 12, height = 18)
ggsave(filename = file.path(outpath, paste0(prefix, "weight_DEG_pathways.pdf")), plot = plot3, width = 12, height = 18)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"weight_DEG_pathways.png (pdf)"), " has been generated successfully!\n", sep = "")

i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
p5 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = group_names[i], width = 10, height = 30)
p6 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = group_names[i+1], width = 10, height = 30)
png(filename = file.path(outpath, paste0(prefix, "outgoing_pathway_heatmap.png")), res = 300, width = 12, height = 20, units = "in")
draw(p5 + p6, ht_gap = unit(1, "cm"))
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "outgoing_pathway_heatmap.pdf")), width = 12, height = 20)
draw(p5 + p6, ht_gap = unit(1, "cm"))
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"outgoing_pathway_heatmap.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.7 LRs comparison
png(filename = file.path(outpath, paste0(prefix, "LRs_comparison_bubble.png")), res = 300, width = 8, height = 12, units = "in")
netVisual_bubble(cellchat_merged, sources.use = 1, targets.use = c(1:2), comparison = c(1, 2), angle.x = 45)
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "LRs_comparison_bubble.pdf")), width = 8, height = 12)
netVisual_bubble(cellchat_merged, sources.use = 1, targets.use = c(1:2), comparison = c(1, 2), angle.x = 45)
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"LRs_comparison_bubble.png (pdf)"), " has been generated successfully!\n", sep = "")

p7 <- netVisual_bubble(cellchat_merged, sources.use = 1, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", group_names[2]), angle.x = 45, remove.isolate = T)
p8 <- netVisual_bubble(cellchat_merged, sources.use = 1, targets.use = c(1:2),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Increased signaling in ", group_names[1]), angle.x = 45, remove.isolate = T)

png(filename = file.path(outpath, paste0(prefix, "LRs_comparison_bubble_changed.png")), res = 300, width = 16, height = 8, units = "in")
p7 + p8
dev.off()
pdf(file = file.path(outpath, paste0(prefix, "LRs_comparison_bubble_changed.pdf")), width = 16, height = 8)
p7 + p8
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"LRs_comparison_bubble_changed.png (pdf)"), " has been generated successfully!\n", sep = "")

# 4.8 Lrs expression comparison
p_dist <- rankSimilarity(cellchat_merged, type = "functional")
p_dist_data <- as.data.frame(p_dist$data)
pathway <- p_dist_data %>% 
  arrange(desc(dist)) %>%
  head(1) %>%
  pull(name) %>% 
  as.character()
cellchat_merged@meta$datasets = factor(cellchat_merged@meta$datasets, levels = c(group_names[1], group_names[2])) # set factor level
p9 <- plotGeneExpression(cellchat_merged, signaling = pathway, split.by = "datasets", colors.ggplot = T, type = "violin") + 
  scale_fill_manual(name = "Dataset", 
                    values = setNames(c("red", "blue"), group_names),
                    labels = c(group_names[1], group_names[2]))
ggsave(filename = file.path(outpath, paste0(prefix, pathway, "_LRs_expression_comparison.png")), plot = p9, width = 16, height = 10, dpi = 300, bg = "white")
ggsave(filename = file.path(outpath, paste0(prefix, pathway, "_LRs_expression_comparison.pdf")), plot = p9, width = 8, height = 6)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix,"LRs_expression_comparison_",pathway,".png (pdf)"), " has been generated successfully!\n", sep = "")
message("\n=========================All Done !============================\n")

