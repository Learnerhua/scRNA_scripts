#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
rm(list = ls())
options(timeout = 1000)

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
A script for cell communication analysis of single cell transcriptomic data.
Example: scRNA_cellchat_single -I test.rds -AC celltype -S mouse -TP truncatedMean -TM 0.1 -T 40 -P test_ -O test
"
))

# ---- Required arguments ----
parser$add_argument("-I", "--input_rds",
                    help = "Path to an RDS file containing a Seurat object for single cell transcriptomics. (required)",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-AC", "--annotation_column",
                    help = "Column name in metadata used as cell identity label (required, e.g., 'celltype').",
                    type = "character",
                    required = TRUE
)
parser$add_argument("-S", "--species",
                    help = "Species used in the analysis (required, 'human' or 'mouse').",
                    type = "character",
                    required = TRUE
)

# ---- Optional arguments ----
parser$add_argument("-TP", "--type",
                    help = "Methods for computing the average gene expression per cell group. By default = 'triMean', producing fewer but stronger interactions; When setting ‘type = 'truncatedMean', a value should be assigned to 'trim', producing more interactions.",
                    type = "character",
                    default = "triMean"
)
parser$add_argument("-TM", "--trim",
                    help = "the fraction (0 to 0.25) of observations to be trimmed from each end of x before the mean is computed.",
                    type = "double",
                    default = 0.1
)
parser$add_argument("-NP", "--nPatterns",
                    help = "Number of patterns to identify in the cell-cell communication network which can be infereed from the selectK plot.",
                    type = "integer"
)
parser$add_argument("-MC", "--min_cell",
                    help = "The minmum number of cells required in each cell group for cell-cell communication (default: 10).",
                    type = "integer",
                    default = 10
)
parser$add_argument("-T", "--threads",
                    help = "Number of threads to use in this analysis (default: 24).",
                    type = "double",
                    default = 24
)
parser$add_argument("-P", "--prefix",
                    help = "Prefix for output files (optional, default: none)",
                    type = "character",
                    default = "")
parser$add_argument("-O", "--outpath",
                    help = "Output directory (optional, default: current directory)",
                    type = "character",
                    default = getwd())

# ---- Parse arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}


opts <- parser$parse_args()
input_rds <- opts$input_rds
annotation_column <- opts$annotation_column
species <- opts$species
type <- opts$type
trim <- opts$trim
min_cell <- opts$min_cell
nPatterns <- opts$nPatterns
threads <- opts$threads
prefix <- opts$prefix
outpath <- opts$outpath
file_n <- 0
plot_n <- 0
options(future.globals.maxSize = threads * 5 * 1024^3)  # Set max size for future package

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

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(CellChat,quietly = TRUE)){install.packages("CellChat")}
  if(!require(NMF,quietly = TRUE)){install.packages("NMF")}
  if(!require(ggalluvial,quietly = TRUE)){install.packages("ggalluvial")}
})

if (!file.exists(file.path(outpath, paste0(prefix,"cellchat.rds")))){
    # 3. Load data
    seu <- readRDS(input_rds)
    cat("\nSeurat objecct loaded.\n")
    seu <- NormalizeData(seu)
    seu@meta.data$samples <- seu@meta.data$orig.ident
    seu@meta.data$samples <- as.factor(seu@meta.data$samples)
    cellchat <- createCellChat(object = seu, meta = seu@meta.data, group.by = annotation_column, datatype = "RNA")

    # 4. Computing the communication probability
    if (species == "human") {
        CellChatDB <- CellChatDB.human
    } else if (species == "mouse") {
        CellChatDB <- CellChatDB.mouse
    }else{
        stop("Species must be either 'human' or 'mouse'.")
    }
    CellChatDB.use <- subsetDB(CellChatDB)
    cellchat@DB <- CellChatDB.use
    cellchat <- subsetData(cellchat)
    future::plan("multicore", workers = threads) 
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)
    if (species == "mouse") {
        cellchat<- smoothData(cellchat, adj = PPI.mouse)
    } else if (species == "human") {
        cellchat<- smoothData(cellchat, adj = PPI.human)
    }
    cellchat <- computeCommunProb(
        cellchat,
        type = type, trim = trim,
        population.size = TRUE,
        raw.use = FALSE
    )
    cellchat <- filterCommunication(cellchat, min.cells = min_cell) 
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    net_df <- subsetCommunication(cellchat, slot.name = "net")
    write.csv(net_df, file.path(outpath, paste0(prefix, "communication_genes.csv")), row.names = FALSE, quote = FALSE)
    file_n <- file_n + 1
    cat("\nFile ", file_n, " : ", paste0(prefix, "communication_genes.csv"), " saved successfully !\n", sep = "")
    netP_df <- subsetCommunication(cellchat, slot.name = "netP")
    write.csv(netP_df, file.path(outpath, paste0(prefix, "communication_pathways.csv")), row.names = FALSE, quote = FALSE)
    file_n <- file_n + 1
    cat("\nFile ", file_n, " : ", paste0(prefix, "communication_pathways.csv"), " saved successfully !\n", sep = "")
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    saveRDS(cellchat, file.path(outpath, paste0(prefix,"cellchat.rds")))
    file_n <- file_n + 1
    cat("\nFile ", file_n, " : ", paste0(prefix,"cellchat.rds"), " saved successfully !\n", sep = "")
} else {
    cellchat <- readRDS(file.path(outpath, paste0(prefix,"cellchat.rds")))
    cat("\nCellChat object loaded from existing file.\n")
    netP_df <- subsetCommunication(cellchat, slot.name = "netP")
}

# 5. Visualization
# 5.1 Circle plot
groupSize <- as.numeric(table(cellchat@idents))
circle_plot <- function(){
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction strength")
}
png(file.path(outpath, paste0(prefix, "circle_overall.png")), res=300, units="in", width = 16, height = 8)
circle_plot()
dev.off()
pdf(file.path(outpath, paste0(prefix,"circle_overall.pdf")), width = 16, height = 8)
circle_plot()
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "circle_overall.png (pdf) has been generated successfully !\n", sep = "")

mat_count <- cellchat@net$count
circle_split_count_plot <- function(){
  par(mfrow = c(ceiling(nrow(mat_count) / 3), 3),xpd = TRUE)
for (i in 1:nrow(mat_count)) {
  mat2 <- matrix(0, nrow = nrow(mat_count), ncol = ncol(mat_count),
                 dimnames = dimnames(mat_count))
  mat2[i, ] <- mat_count[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                  #  edge.weight.max = max(mat_count),
                   title.name = rownames(mat_count)[i])
}
}
png(file.path(outpath, paste0(prefix, "circle_split_count.png")), width = 12, height = ceiling(nrow(mat_count)/3) * 4, res = 300, units = "in")
circle_split_count_plot()
dev.off()
pdf(file.path(outpath, paste0(prefix, "circle_split_count.pdf")), width = 12, height = ceiling(nrow(mat_count)/3) * 4)
circle_split_count_plot()
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "circle_split_count.png (pdf) has been generated successfully !\n", sep = "")

mat_weight <- cellchat@net$weight
circle_split_weight_plot <- function(){
  par(mfrow = c(ceiling(nrow(mat_weight) / 3), 3),xpd = TRUE)
for (i in 1:nrow(mat_weight)) {
  mat2 <- matrix(0, nrow = nrow(mat_weight), ncol = ncol(mat_weight),
                 dimnames = dimnames(mat_weight))
  mat2[i, ] <- mat_weight[i, ]
  netVisual_circle(mat2,
                   vertex.weight = groupSize,
                   weight.scale = TRUE,
                  #  edge.weight.max = max(mat_weight),
                   title.name = rownames(mat_weight)[i])
}
}
png(file.path(outpath, paste0(prefix, "circle_split_weight.png")), width = 12, height = ceiling(nrow(mat_weight)/3) * 4, res =  300, units = "in")
circle_split_weight_plot()
dev.off()
pdf(file.path(outpath, paste0(prefix, "circle_split_weight.pdf")), width = 12, height = ceiling(nrow(mat_weight)/3) * 4)
circle_split_weight_plot()
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "circle_split_weight.png (pdf) has been generated successfully !\n", sep = "")

# 5.2 Pathway diaplay
top9_pathways <- netP_df %>%
  filter(pval < 0.05) %>%   
  arrange(desc(prob)) %>% 
  distinct(pathway_name, .keep_all = TRUE) %>% 
  slice_head(n = 9) %>% 
  pull(pathway_name)
if(length(top9_pathways) < 9) {
  warning(
    "Less than 9 significant pathways found. Only ",
    length(top9_pathways),
    " pathways will be plotted.")
}
SignalPath_circle_plot <- function(){
  par(mfrow = c(3, 3), xpd=TRUE)
  for (pathway in top9_pathways){
    par(xpd = TRUE)
    netVisual_aggregate(
      object = cellchat,  
      signaling = pathway,    
      signaling.name = pathway,
      layout = "chord",    
      vertex.weight = groupSize,
      weight.scale = TRUE,   
      edge.weight.max = max(cellchat@net$weight), 
      vertex.label.cex = 1.2,  
      title.space = 0.05, 
      pt.title = 16   
    )
  }
}
png(file.path(outpath, paste0(prefix, "top9_signalPath_chord.png")), res=300, units="in", width = 15, height = 15)
SignalPath_circle_plot()
dev.off()
pdf(file.path(outpath, paste0(prefix, "top9_signalPath_chord.pdf")), width = 15, height = 15)
SignalPath_circle_plot()
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "top9_signalPath_chord.png (pdf) has been generated successfully !\n", sep = "")

# 5.3 Displaying specific pathway
pathways.show <- top9_pathways[1]
pathway_contribution <- function() {
  netAnalysis_contribution(
  object = cellchat,
  signaling = pathways.show,
  title = paste0("Contribution of each L-R pair in ",
  pathways.show, " signaling pathway"),
  font.size = 12,
  font.size.title = 16)
}
png(file.path(outpath, paste0(prefix, pathways.show, "_LR_distribution.png")), res=300, units="in", width = 16, height = 16)
pathway_contribution()
dev.off()
pdf(file.path(outpath, paste0(prefix, pathways.show, "_LR_distribution.pdf")), width = 16, height = 16)
pathway_contribution()
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, pathways.show, "_LR_distribution.png (pdf) has been generated successfully !\n", sep = "")

p1 <- plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = T)
ggsave(file.path(outpath, paste0(prefix, pathways.show, "_LR_expr.png")), p1, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, pathways.show, "_LR_expr.pdf")), p1, width = 8, height = 8, bg = "white")
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, pathways.show, "_LR_expr.png (pdf) has been generated successfully !\n", sep = "")

LR_df <- netAnalysis_contribution(
  object = cellchat,
  signaling = pathways.show,
  return.data = TRUE
)
df_sorted <- LR_df$LR.contribution[order(LR_df$LR.contribution$contribution, decreasing = TRUE), ]


png(file.path(outpath, paste0(prefix, pathways.show,"_role_network.png")), res=300,units="in",width = 6, height = 4)
par(xpd = TRUE)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, font.size = 10 )
dev.off()
pdf(file.path(outpath, paste0(prefix, pathways.show,"_role_network.pdf")), width = 6, height = 4)
par(xpd = TRUE)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, font.size = 10 )
dev.off()

p2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
ggsave(file.path(outpath, paste0(prefix, pathways.show,"_role_scatter.png")), p2, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, pathways.show,"_role_scatter.pdf")), p2, width = 8, height = 8, bg = "white")
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, pathways.show, "_role_network.png (pdf) has been generated successfully !\n", sep = "")

for (LR.show in rownames(df_sorted)[1:3]){
  png(file.path(outpath, paste0(prefix, LR.show, "_circle.png")), res=300, units="in", width = 8, height = 8)
  par(xpd = TRUE)
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", vertex.weight = groupSize, weight.scale = TRUE)
  dev.off()
  pdf(file.path(outpath, paste0(prefix, LR.show, "_circle.pdf")), width = 8, height = 8)
  par(xpd = TRUE)
  netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle", vertex.weight = groupSize, weight.scale = TRUE)
  dev.off()
  plot_n <- plot_n + 1
  cat("Plot ", plot_n, " : ", prefix, LR.show, "_circle.png (pdf) has been generated successfully !\n", sep = "")
}

# 5.4 Bubble plot for overall interactions
p3 <- netVisual_bubble(cellchat, sources.use = 1:length(levels(cellchat@idents)), targets.use = 1:length(levels(cellchat@idents)), remove.isolate = T)
ggsave(file.path(outpath, paste0(prefix, "allLRs_bubble.png")), p3, width = 24 ,height = 24 ,limitsize = F, bg = "white", dpi = 300)
ggsave(file.path(outpath, paste0(prefix, "allLRs_bubble.pdf")), p3, width = 24 ,height = 24 ,limitsize = F)
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "allLRs_bubble.png (pdf) has been generated successfully !\n", sep = "")

# 5.5 Find communication patterns
if (is.null(nPatterns) &&
    (!file.exists(file.path(outpath, paste0(prefix, "selectK_outcoming.png"))) || 
    !file.exists(file.path(outpath, paste0(prefix, "selectK_outcoming.pdf"))))
    ){
  cat("\nGenerating selectK plots for outgoing communication patterns, it will take a while ...\n")
  png(file.path(outpath, paste0(prefix, "selectK_outcoming.png")), res=300, units="in", width = 12, height = 6)
  print(selectK(cellchat, pattern = "outgoing"))
  dev.off()
  pdf(file.path(outpath, paste0(prefix, "selectK_outcoming.pdf")), width = 12, height = 6)
  print(selectK(cellchat, pattern = "outgoing"))
  dev.off()
  plot_n <- plot_n + 1
  cat("Plot ", plot_n, " : ", prefix, "selectK_outcoming.png (pdf) has been generated successfully !\n", sep = "")
} else{
  cat("Plot : ", prefix, "selectK_outcoming.png (pdf) already exists, skipping...\n", sep = "")
}

if (is.null(nPatterns)) {
  stop("Please check the selectK plots and specify the number of patterns to identify using the -NP argument.")
}

png(file.path(outpath, paste0(prefix, "outgoing_patterns_heatmap.png")), res=300, units="in", width = 16, height = 16)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 10, height = 20)
dev.off()
pdf(file.path(outpath, paste0(prefix, "outgoing_patterns_heatmap.pdf")), width = 16, height = 16)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, width = 10, height = 20)
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "outgoing_patterns_heatmap.png (pdf) has been generated successfully !\n", sep = "")

png(file.path(outpath, paste0(prefix, "outgoing_patterns_river.png")), res=300, units="in", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
pdf(file.path(outpath, paste0(prefix, "outgoing_patterns_river.pdf")), width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "outgoing_patterns_river.png (pdf) has been generated successfully !\n", sep = "")

png(file.path(outpath, paste0(prefix, "outgoing_patterns_dot.png")), res=300, units="in", width = 12, height = 12)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()
pdf(file.path(outpath, paste0(prefix, "outgoing_patterns_dot.pdf")), width = 12, height = 12)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()
plot_n <- plot_n + 1
cat("Plot ", plot_n, " : ", prefix, "outgoing_patterns_dot.png (pdf) has been generated successfully !\n", sep = "")
message("\n=========================All Done !============================\n")
