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
    A script for pseudotime and trajectory analysis of single-cell and spatial transcriptomic data based on Monocle3.
    Example: scRNA_spRNA_Monocle3 -I seo_nor_umap.rds -DT SP -C celltype -G Hepatocyte -R1 2e-4 -R2 0.5 -T 40 -P test. -O test
  '
  )
)

#Required
parser$add_argument(
  "-I",
  "--input_rds",
  dest = "input_rds",
  type = "character",
  required = TRUE,
  help = "Path to an RDS file containing a Seurat object for pseudotime and trajectory analysis of single-cell or spatial transcriptomic data."
)
parser$add_argument(
  "-DT",
  "--data_type",
  dest = "data_type",
  type = "character",
  required = TRUE,
  help = "Data type of the input file. Use 'SC' for single-cell data, or 'SP' for spatial transcriptomic data. 
          If 'SP' is specified, the metadata must contain 'x' and 'y' columns to provide spatial annotation information."
)


#Optioanl
parser$add_argument(
  '-C',
  '--cluster_column',
  dest = 'cluster_column',
  type = "character",
  help = "The metadata column name used to specify cell clusters or annotations for trajectory root cell selection. 
          Must be used together with -G/--group."
)
parser$add_argument(
  '-G',
  '--group',
  dest = 'group',
  type = "character",
  help = "The group value within the specified metadata column (-C) that defines the root cell population for pseudotime analysis. Must be used together with -C/--cluster_column."
)
parser$add_argument(
  '-R1',
  '--resolution1',
  dest = 'res1',
  help = 'resolution for clustering in Monocle3.(Optional, default: 1e-5)',
  type = "double",
  default = 1e-5
)
parser$add_argument(
  '-R2',
  '--resolution2',
  dest = 'res2',
  help = 'resolution for module analysis in Monocle3.(Optional, default: 0.5)',
  type = "double",
  default = 0.5
)

parser$add_argument(
  '-T',
  '--threads',
  dest = 'threads',
  type = "double",
  default = 10,
  help = "Number of threads to use for Moran's I test and module analysis (default: 10). Adjust this parameter to control parallel computation and speed."
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
data_type <- opts$data_type
cluster_column <- opts$cluster_column
group <- opts$group
res1 <- opts$res1
res2 <- opts$res2
threads <- opts$threads
prefix <- opts$prefix
outpath <- opts$outpath

if (is.null(cluster_column) != is.null(group)) {
  cat("Error: Both -C/--cluster_column and -G/--group must be provided together for root cell specification.\n")
  quit(status = 1)
}

message("\n========================Monocle3 anaysis=======================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 KMHD. All Rights Reserved.")
message("\n===============================================================\n")

if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(monocle3,quietly = TRUE)){install.packages("monocle3")}
  if(!require(SeuratWrappers,quietly = TRUE)){install.packages("SeuratWrappers")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
  if(!require(pheatmap,quietly = TRUE)){install.packages("pheatmap")}
})

# 3. Load data
n_file <- 0
n_plot <- 0
seo <- readRDS(input_rds)
if (!("umap" %in% names(seo@reductions))){
  warning("UMAP information not found in the Seurat object. Re-normalizing the data, identifying variable features, scaling, extracting the top 30 principal components, and running UMAP again. 
          It is recommended to input data that has already undergone UMAP dimensionality reduction to save time.")
  seo <- seo %>% 
    NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
    ScaleData(features = VariableFeatures(seo))
  seo <- RunPCA(seo, ndims = 30, features = VariableFeatures(object = seo))
  seo <- RunUMAP(seo, dims = 1:30, verbose = FALSE)
  saveRDS(seo, file = file.path(outpath, paste0(prefix, "seo_nor_umap.rds")))
  n_file <- n_file + 1
  cat("\nFile ", n_file, " : ", prefix, "seo_nor_umap.rds saved successfully!\n", sep = "")
}

seo_cds <- as.cell_data_set(seo)
seo_cds <- cluster_cells(seo_cds, resolution = res1)
cat("\nThe clustering results at resolution", res1, ":\n")
print(table(clusters(seo_cds)))
cat("\nThe partition results at resolution", res1, ":\n")
print(table(partitions(seo_cds)))
colData(seo_cds)$monocle_cluster <- clusters(seo_cds)
colData(seo_cds)$monocle_partition <- partitions(seo_cds)
#table(colData(seo_cds)$monocle_cluster, colData(seo_cds)$celltype)
cat("\nLearning trajectory graph...\n")
seo_cds <- learn_graph(seo_cds, use_partition = TRUE, verbose = FALSE)
saveRDS(seo_cds, file.path(outpath, paste0(prefix, "seo_cds.rds")))
n_file <- n_file + 1
cat("\nFile ", n_file, " : ", prefix, "seo_cds.rds saved successfully!\n", sep = "")

p1 <- plot_cells(seo_cds,
                 reduction_method = "UMAP",
                 color_cells_by = "cluster",
                 group_cells_by = "cluster",
                 label_cell_groups = TRUE,
                 group_label_size = 6,
                 label_groups_by_cluster= FALSE,
                 label_leaves= FALSE,
                 label_branch_points= FALSE,
                 label_roots = FALSE,
                 trajectory_graph_color = "black")
ggsave(file.path(outpath, paste0(prefix, "umap_cluster.png")), p1, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "umap_cluster.pdf")), p1, width = 12, height = 10)
n_plot <- n_plot + 1
cat("\nPlot ", n_plot, " : ", prefix, "umap_cluster.png(pdf) has been generated successfully !\n", sep = "")

# 4. Specify root cells.
if (is.null(cluster_column) || is.null(group)) {
  warning("No root cell specification provided. Selecting cells from cluster '1' for pseudotime analysis. 
          It is recommended to specify root cells using -C/--cluster_column and -G/--group options.\n")
  seo_cds <- order_cells(seo_cds, root_cells = colnames(seo_cds[,colData(seo_cds)$monocle_cluster == 1]))
} else {
  cat("\nSpecifying root cells from group '", group, "' in metadata column '", cluster_column, "' for pseudotime analysis...\n", sep = "")
  seo_cds <- order_cells(seo_cds, root_cells = colnames(seo_cds[,colData(seo_cds)[[cluster_column]] == group]))
}
colData(seo_cds)$pseudotime <- principal_graph_aux(seo_cds)$UMAP$pseudotime

p2 <- plot_cells(seo_cds,
                 reduction_method = "UMAP",
                 color_cells_by = "pseudotime",
                 group_cells_by = "cluster",
                 label_cell_groups = FALSE,
                 label_groups_by_cluster= FALSE,
                 label_leaves= FALSE,
                 label_branch_points= FALSE,
                 label_roots = FALSE,
                 trajectory_graph_color = "black")+
  scale_color_viridis_c(option = "viridis", name = "Pseudotime")
ggsave(file.path(outpath, paste0(prefix, "umap_pseudotime.png")), p2, width = 12, height = 10, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "umap_pseudotime.pdf")), p2, width = 12, height = 10)
n_plot <- n_plot + 1
cat("\nPlot ", n_plot, " : ", prefix, "umap_pseudotime.png(pdf) has been generated successfully !\n", sep = "")

# 5. Moran's I test for spatial autocorrelation of gene expression
cat("\nPerforming Moran's I test for spatial autocorrelation of gene expression...\n")
moran <- graph_test(seo_cds, neighbor_graph = "principal_graph", cores = threads)
moran <- rownames_to_column(moran, var = "gene")
write.csv(moran, file.path(outpath, paste0(prefix, "moran_result.csv")), row.names = FALSE, quote = FALSE)
n_file <- n_file + 1
cat("\nFile ", n_file, " : ", prefix, "moran_result.csv saved successfully!\n", sep = "")
moran_DEGs <- moran %>%
  filter(p_value < 0.05 & q_value < 0.1 ) %>%
  arrange(desc(morans_I)) %>%
  slice_head(n = 12)  %>% 
  pull(gene)
if (length(moran_DEGs) > 0){
  p3 <- plot_genes_in_pseudotime(seo_cds[moran_DEGs,],
                                 color_cells_by = "pseudotime",
                                 label_by_short_name = FALSE,
                                 min_expr = 0,
                                 ncol = 3,
                                 cell_size = 0.1,
                                 vertical_jitter = 0.2,
                                 horizontal_jitter = 0.1
  ) +
    scale_color_viridis_c(option = "viridis", name = "Pseudotime")
  ggsave(file.path(outpath, paste0(prefix, "moran_DEGs_pseudotime.png")), p3, width = 9, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "moran_DEGs_pseudotime.pdf")), p3, width = 9, height = 12)
  n_plot <- n_plot + 1
  cat("\nPlot ", n_plot, " : ", prefix, "moran_DEGs_pseudotime.png(pdf) has been generated successfully !\n", sep = "")
} else {
  cat("\nNo significant genes found in Moran's I test with the given thresholds.\n")
}

# 6. Module analysis
moran_DEGs_all <- moran %>%
  filter(p_value < 0.05 & q_value < 0.1 & morans_I > 0.1) %>%
  pull(gene)
cat("\nNumber of significant genes for module analysis:", length(moran_DEGs_all), "\n")
if (length(moran_DEGs_all) > 30){
  seo_cds <- preprocess_cds(seo_cds, method = "PCA")
  cat("\nPerforming module analysis on significant genes...\n")
  gene_module <- find_gene_modules(seo_cds[moran_DEGs_all,], resolution = res2, cores = threads)
  cat("\nNumber of gene modules identified:", length(unique(gene_module$module)), "\n")
  write.csv(gene_module, file.path(outpath, paste0(prefix, "gene_modules.csv")), row.names = FALSE, quote = FALSE)
  n_file <- n_file + 1
  cat("\nFile ", n_file, " : ", prefix, "gene_modules.csv saved successfully!\n", sep = "")
  cellGroup <- tibble(cell = row.names(colData(seo_cds)), cell_group = colData(seo_cds)[[cluster_column]])
  agg_mat <- aggregate_gene_expression(seo_cds, gene_module, cellGroup)
  row.names(agg_mat) <- str_c("Module ", row.names(agg_mat))
  png(file.path(outpath, paste0(prefix, "module_heatmap.png")), width = 1800, height = 1200, res = 300)
  pheatmap(agg_mat, scale="column"  ,clustering_method = "ward.D2")
  dev.off()
  pdf(file.path(outpath, paste0(prefix, "module_heatmap.pdf")), width = 12, height = 8)
  pheatmap(agg_mat, scale="column"  ,clustering_method = "ward.D2")
  dev.off()
  n_plot <- n_plot + 1
  cat("\nPlot ", n_plot, " : ", prefix, "module_heatmap.png(pdf) has been generated successfully !\n", sep = "")
} else {
  cat("\nNot enough significant genes for module analysis. At least 30 genes with morans_I > 0.1, p_value < 0.05, and q_value < 0.1 are required.\n")
}

# 7. Spatial pseudotime plot for spatial transcriptomic data
if (data_type == "SP") {
  if (!all(c("x", "y") %in% colnames(colData(seo_cds)))) {
    cat("Error: For spatial transcriptomic data, the metadata must contain 'x' and 'y' columns for spatial coordinates.\n")
    quit(status = 1)
  }
  plot_data <- data.frame(
    x = colData(seo_cds)$x,
    y = colData(seo_cds)$y,
    pseudotime = colData(seo_cds)$pseudotime
  )
  p4 <- ggplot(plot_data, aes(x = x, y = y, color = pseudotime)) +
    geom_point(size = 0.5) + scale_color_viridis_c(option = "viridis", name = "Pseudotime")  +
    theme_classic() +
    labs(
      title = "Spatial Pseudotime Distribution",
      x = "X Coordinate",
      y = "Y Coordinate"
    ) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.4, "cm"),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  ggsave(file.path(outpath, paste0(prefix, "spatial_pseudotime.png")), p4, width = 12, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, "spatial_pseudotime.pdf")), p4, width = 12, height = 10)
  n_plot <- n_plot + 1
  cat("\nPlot ", n_plot, " : ", prefix, "spatial_pseudotime.png(pdf) has been generated successfully !\n", sep = "")
}
message("\n=========================All Done !============================\n")




