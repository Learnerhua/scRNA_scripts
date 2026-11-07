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
A script for celltype proportion analysis of single-cell data.
Example: scRNA_spRNA_cellPro -I seo_annotated.rds -SC orig.ident -CC celltype -GC Group -GO WT,KO -EC Unknown -P pre_ -O output_dir
"
))
# ---- Required arguments ----
parser$add_argument("-I", "--input_rds",
                  help = "Path to an RDS file containing a Seurat object for celltype proportion analysis.",
                  type = "character",
                  required = TRUE
                  )
parser$add_argument("-SC", "--sample_column",
                  help = "Column name for sample identity in metadata (e.g., orig.ident)",
                  type = "character",
                  required = TRUE
                  )
parser$add_argument("-CC", "--celltype_column",
                  help = "Column name for celltype identity in metadata (e.g., celltype)",
                  type = "character",
                  required = TRUE
                  )
parser$add_argument("-SCC", "--sub_celltype_column",
                    help = "Column name for sub-celltype identity in metadata (e.g., sub_celltype)",
                    type = "character"
)
parser$add_argument("-GC", "--group_column",
                  help = "Column name for group identity in metadata (e.g., Group)",
                  type = "character",
                  required = TRUE
                  )
# ---- Optional arguments ----
parser$add_argument("-GO", "--group_order",
                  help = "Order of groups (comma or space separated, e.g., WT,KO)"
                  )
parser$add_argument("-EC", "--end_celltype",
                  help = "Celltype to be placed at the end in plots (optional, default: NULL)",
                  type = "character",
                  default = NULL)
parser$add_argument("-BW", "--bar_plot_width",
                  help = "Width of the bar plot (optional, default: 10)",
                  type = "numeric",
                  default = 10)
parser$add_argument( "-BH", "--bar_plot_height",
                  help = "Height of the bar plot (optional, default: 12)",
                  type = "numeric",
                  default = 12)
parser$add_argument("-SW", "--sankey_plot_width",
                  help = "Width of the sankey plot (optional, default: 12)",
                  type = "numeric",
                  default = 12)
parser$add_argument("-SH", "--sankey_plot_height",
                  help = "Height of the sankey plot (optional, default: 6)",
                  type = "numeric",
                  default = 6)
parser$add_argument("-DW", "--DE_plot_width",
                    help = "Width of the celltype DE plot (optional, default: 14)",
                    type = "numeric",
                    default = 14)
parser$add_argument("-DH", "--DE_plot_height",
                    help = "Height of the celltype DE plot (optional, default: 8)",
                    type = "numeric",
                    default = 8)
parser$add_argument("-SDW", "--sub_DE_plot_width",
                    help = "Width of the sub-celltype DE plot (optional, default: 36)",
                    type = "numeric",
                    default = 36)
parser$add_argument("-SDH", "--sub_DE_plot_height",
                    help = "Height of the sub-celltype DE plot (optional, default: 28)",
                    type = "numeric",
                    default = 28)
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
sample_column <- opts$sample_column
celltype_column <- opts$celltype_column
sub_celltype_column <- opts$sub_celltype_column
group_column <- opts$group_column
end_celltype <- opts$end_celltype
bar_plot_width <- opts$bar_plot_width
bar_plot_height <- opts$bar_plot_height
sankey_plot_width <- opts$sankey_plot_width
sankey_plot_height <- opts$sankey_plot_height
DE_plot_width <- opts$DE_plot_width
DE_plot_height <- opts$DE_plot_height
sub_DE_plot_width <- opts$sub_DE_plot_width
sub_DE_plot_height <- opts$sub_DE_plot_height
prefix <- opts$prefix
outpath <- opts$outpath

if (!is.null(opts$group_order)) {
  group_order <- unlist(strsplit(opts$group_order, "[, ]+"))
} else {
  group_order <- NULL
}


message("\n===================Celltype Proportion anaysis=================\n")
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
  if(!require(RColorBrewer,quietly = TRUE)){install.packages("RColorBrewer")}
  if(!require(ggalluvial,quietly = TRUE)){install.packages("ggalluvial")}
  if(!require(ggpubr,quietly = TRUE)){install.packages("ggpubr")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
})


# 3.Load data
file_n <- 0
plot_n <- 0
seo <- readRDS(input_rds)
meta_data <- seo@meta.data
if (!is.null(group_order)){
meta_data <- meta_data %>% 
  filter(!!sym(group_column) %in% group_order)
}
get_celltype_count <- function(data, sample_col, celltype_col) {
  data %>%
    count(!!sym(sample_col), !!sym(celltype_col)) %>%
    ungroup() %>%
    pivot_wider(
      names_from = !!sym(celltype_col),
      values_from = n,
      values_fill = 0
    ) %>%
    rename(Sample = !!sym(sample_col)) %>%
    select(Sample, sort(colnames(.)[-1]))
}
celltype_count_sample <- get_celltype_count(meta_data, sample_column, celltype_column)
write.csv(celltype_count_sample, file.path(outpath, paste0(prefix, celltype_column , "_count_per_sample.csv")), row.names = F, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, celltype_column, "_count_per_sample.csv"), " saved successfully !\n", sep = "")

get_celltype_prop <- function(data, sample_col, celltype_col) {
  data %>%
    count(!!sym(sample_col), !!sym(celltype_col)) %>%
    group_by(!!sym(sample_col)) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(
      names_from = !!sym(celltype_col),
      values_from = prop,
      values_fill = 0
    ) %>%
    rename(Sample = !!sym(sample_col)) %>%
    select(Sample, sort(colnames(.)[-1]))
}
celltype_prop_sample <- get_celltype_prop(meta_data, sample_column, celltype_column)
write.csv(celltype_prop_sample, file.path(outpath, paste0(prefix, celltype_column
                                                          , "_proportion_per_sample.csv")), row.names = F, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, celltype_column, "_proportion_per_sample.csv"), " saved successfully !\n", sep = "")

celltype_prop <- meta_data %>%
  group_by(!!sym(sample_column), !!sym(celltype_column)) %>%
  summarise(
    cell_count = n(),
    .groups = "drop"
  ) %>%
  group_by(!!sym(sample_column)) %>%
  mutate(
    proportion = cell_count / sum(cell_count)
  ) %>%
  ungroup() %>%
  rename(sample = !!sym(sample_column)) %>%
  left_join(
    meta_data %>% distinct(!!sym(sample_column), !!sym(group_column)) %>%
      rename(sample = !!sym(sample_column), Group = !!sym(group_column)),
    by = "sample"
  )%>%
  relocate(Group, .before = sample)
write.csv(celltype_prop, file.path(outpath, paste0(prefix, celltype_column, "_proportion_stat.csv")), row.names = F, quote = F)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, celltype_column, "_proportion_stat.csv"), " saved successfully !\n", sep = "")

# 3. Bar plot
global_prop <- meta_data %>%
  count(!!sym(celltype_column)) %>%
  mutate(prop = n / sum(n)) %>%
  arrange(desc(prop))

if (!is.null(end_celltype)) {
  celltype_levels <- global_prop %>%
    arrange(if_else(!!sym(celltype_column) == end_celltype, Inf, -prop)) %>%
    pull(!!sym(celltype_column))
} else {
  celltype_levels <- global_prop %>%
    arrange(-prop) %>%
    pull(!!sym(celltype_column))
}

n_colors <- length(celltype_levels)
palette_colors <- brewer.pal(min(n_colors, 12), "Set3")
if (n_colors > 12) {
  palette_colors <- colorRampPalette(palette_colors)(n_colors)
}
names(palette_colors) <- celltype_levels

if(!is.null(group_order)) {
  meta_data <- meta_data %>%
    mutate(!!sym(celltype_column) := factor(!!sym(celltype_column), levels = celltype_levels)) %>% 
    mutate(!!sym(group_column) := factor(!!sym(group_column), levels = group_order))
} else {
  meta_data <- meta_data %>%
    mutate(!!sym(celltype_column) := factor(!!sym(celltype_column), levels = celltype_levels))
}

facet_theme <- theme(
  strip.background = element_rect(fill = "grey90", colour = "black", linewidth = 0.8),
  strip.text = element_text(size = 13, face = "bold", colour = "black"),
  strip.text.x = element_text(margin = margin(5, 5, 5, 5))
)

p1 <- ggplot(meta_data, aes(x = !!sym(sample_column), fill = !!sym(celltype_column))) +
  geom_bar(position = "fill", width = 0.9, colour = "black", linewidth = 0.75) +
  facet_wrap(vars(!!sym(group_column)), scales = "free_x") +
  labs(y = "Percentage", x = "Sample", title = "Celltype Composition per Sample") +
  theme_classic(base_size = 14) +
  facet_theme +
  theme(
    panel.grid = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(lineheight = 0.8, face = "bold", hjust = 0.5, size = 20),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(size = 1.2),
    axis.title = element_text(face = "bold",size = 16),
    axis.title.x = element_text(vjust = -6),
    axis.title.y = element_text(vjust = 6),
    axis.text.x = element_text(face = "bold",angle = 45,hjust = 1,size = 12, color = "black"),
    axis.text.y = element_text(size = 12,face = "bold", color = "black"),
    axis.line = element_line(colour = "black", linewidth = 1),
    legend.title = element_text(size=14, color = "black", face="bold"),
    plot.margin = unit(rep(1,4), "cm")
  )+
  scale_y_continuous(expand = c(0,0.01),labels = scales::percent)+
  scale_fill_manual(values = palette_colors)

ggsave(file.path(outpath, paste0(prefix, "bar.png")), p1, dpi = 300, width = bar_plot_width, height = bar_plot_height, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "bar.pdf")), p1, width = bar_plot_width, height = bar_plot_height)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix, "bar.png (pdf)"), " has been generated successfully !\n", sep = "")


# 4. Sankey Plot ----------------------------------------------------------------
celltype_count_group <- meta_data %>%
  count(!!sym(group_column), !!sym(celltype_column)) %>%
  pivot_wider(
    names_from = !!sym(celltype_column),
    values_from = n,
    values_fill = 0 
  ) %>%
  rename(Group = !!sym(group_column)) %>%
  select(Group, sort(colnames(.)[-1]))

write.csv(celltype_count_group, file.path(outpath, paste0(prefix, celltype_column, "_count_per_group.csv")), row.names = F, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, celltype_column, "_count_per_group.csv"), " saved successfully !\n", sep = "")

celltype_prop_group <- meta_data %>%
  count(!!sym(group_column), !!sym(celltype_column)) %>%
  group_by(!!sym(group_column)) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(
    names_from = !!sym(celltype_column),
    values_from = prop,
    values_fill = 0 
  ) %>%
  rename(Group = !!sym(group_column)) %>%
  select(Group, sort(colnames(.)[-1]))
write.csv(celltype_prop_group, file.path(outpath, paste0(prefix, celltype_column, "_proportion_per_group.csv")), row.names = F, quote = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", paste0(prefix, celltype_column, "_proportion_per_group.csv"), " saved successfully !\n", sep = "")

color_max <- brewer.pal.info["Paired", "maxcolors"]
nCluster <- ncol(celltype_prop_group) - 1
base_colors <- brewer.pal(min(nCluster, color_max), "Paired")
if (nCluster > color_max) {
  palette_colors <- colorRampPalette(base_colors)(nCluster)
} else {
  palette_colors <- base_colors
}
names(palette_colors) <- sort(colnames(celltype_prop_group)[-1])

sankey_plot <- function(prop_data, group_order = NULL) {
  # 转换数据为长格式
  plot_data <- prop_data %>%
    pivot_longer(
      cols = -Group,
      names_to = "CellType",
      values_to = "Freq"
    )
  
  # 设置因子水平
  if(!is.null(group_order)) {
    plot_data <- plot_data %>%
      filter(Group %in% group_order) %>%
      mutate(Group = factor(Group, levels = group_order))
  }
  
  plot_data$CellType <- factor(plot_data$CellType, 
                               levels = rev(celltype_levels))
  orig_levels <- celltype_levels
  
  # 绘制桑基图
  ggplot(plot_data,
         aes(x = Group, stratum = CellType, alluvium = CellType,
             y = Freq,
             fill = CellType, label = CellType)) +
    scale_fill_manual(values = palette_colors, breaks = orig_levels) +  # 使用之前定义的颜色
    geom_flow(stat = "alluvium", lode.guidance = "leftright",
              color = "white", size = 1.5) +
    geom_stratum(color = "transparent") +
    ylab("Frequency (%)") +
    theme_classic() +
    coord_flip() +
    scale_y_continuous(labels = scales::percent)+
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "black", size = 12, face = "bold"),
      axis.text.y = element_text(color = "black",size = 14, face = "bold"),
      axis.title = element_text(size = 14, face = "bold"),
      axis.line = element_line(color = "black", size = 1),
      axis.ticks.x = element_line(size=1),
      axis.ticks.length = unit(-4, "pt"),
      axis.ticks.y = element_blank(),
      legend.position = "right"
    )
}

p2 <- sankey_plot(celltype_prop_group, rev(group_order))
ggsave(file.path(outpath, paste0(prefix, "sankey.png")), p2 , width = sankey_plot_width, height = sankey_plot_height, dpi = 300, bg = "white")
ggsave(file.path(outpath, paste0(prefix, "sankey.pdf")), p2 , width = sankey_plot_width, height = sankey_plot_height)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, " : ", paste0(prefix, "sankey.png (pdf)"), " has been generated successfully !\n", sep = "")


# 5.cellPro DE
if (!is.null(group_order) && length(group_order) == 2){
  celltype_prop <- celltype_prop %>%
    filter(Group %in% group_order) %>%
    mutate(Group = factor(Group, levels = group_order)) %>% 
    mutate(celltype = factor(!!sym(celltype_column), levels = celltype_levels))
  p3 <- ggplot(celltype_prop, aes(x = celltype, y = proportion, fill = Group))+
    stat_boxplot(
      geom = "errorbar", 
      position = position_dodge(width = 0.8), 
      width = 0.3, 
      size = 0.8,    
      color = "black"
    ) +
    geom_boxplot(
      position = position_dodge(width = 0.8),  
      width = 0.6,                            
      outlier.shape = NA,                     
      size = 0.3,                             
      color = "black"                        
    ) +
    geom_point(
      color = "black",
      position = position_jitterdodge(jitter.width = 0.2),
      size = 1,
    )+
    scale_fill_manual(
      values = setNames(
        c("#90BFF9", "#FFA0A0"),
        group_order
      ))+
    stat_compare_means(
      aes(group = Group),                     
      label = "p.signif",                     
      method = "wilcox.test",                     
      hide.ns = F
    )+
    labs(
      x = "celltype",
      y = "Proportion",
      fill = "Group"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_text(face = "bold")
    )
  ggsave(file.path(outpath, paste0(prefix, celltype_column, "_DE.png")), p3, width = DE_plot_width, heigh = DE_plot_height, dpi = 300 ,bg = "white")
  ggsave(file.path(outpath, paste0(prefix, celltype_column, "_DE.pdf")), p3, width = DE_plot_width, heigh = DE_plot_height)
  plot_n <- plot_n + 1
  cat("\nPlot ", plot_n, " : ", paste0(prefix, celltype_column, "_DE.png (pdf)"), " has been generated successfully !\n", sep = "")
}
# 6.subcelltype DE
if (!is.null(sub_celltype_column) && !is.null(group_order) && length(group_order) == 2){
  sub_celltype_prop <- meta_data %>%
    group_by(!!sym(sample_column), !!sym(sub_celltype_column)) %>%
    summarise(
      cell_count = n(),
      .groups = "drop" 
    ) %>%
    group_by(!!sym(sample_column)) %>%
    mutate(
      proportion = cell_count / sum(cell_count)
    ) %>%
    ungroup() %>%
    rename(sample = !!sym(sample_column)) %>%
    left_join(
      meta_data %>%
        distinct(!!sym(sample_column), !!sym(group_column), !!sym(sub_celltype_column), !!sym(celltype_column)) %>%
        rename(sample = !!sym(sample_column), Group = !!sym(group_column)),
      by = c("sample", sub_celltype_column)
    ) %>%
    relocate(Group, .before = sample) %>%
    relocate(!!sym(celltype_column), .before = !!sym(sub_celltype_column))
  write.csv(sub_celltype_prop, file.path(outpath, paste0(prefix, sub_celltype_column, "_proportion_stat.csv")), row.names = F, quote = F)
  file_n <- file_n + 1
  cat("\nFile ", file_n, " : ", paste0(prefix, sub_celltype_column, "_proportion_stat.csv"), " saved successfully !\n", sep = "")

  sub_plots_list <- list()
  for (ct in celltype_levels){
    cat("\nProcessing：", ct, "......\n")
    plot_data <- sub_celltype_prop %>% 
      filter(Group %in% group_order) %>%
      mutate(Group = factor(Group, levels = group_order)) %>% 
      filter(!!sym(celltype_column) == ct)
    p <- ggplot(plot_data, aes(x = !!sym(sub_celltype_column), y = proportion, fill = Group))+
      stat_boxplot(
        geom = "errorbar", 
        position = position_dodge(width = 0.8),
        width = 0.3,  
        size = 0.8,   
        color = "black"
      ) +
      geom_boxplot(
        position = position_dodge(width = 0.8), 
        width = 0.6,                            
        outlier.shape = NA,                   
        size = 0.3,                         
        color = "black"                      
      ) +
      geom_point(
        color = "black",
        position = position_jitterdodge(jitter.width = 0.2),
        size = 1
      )+
      scale_fill_manual(
        values = setNames(
          c("#90BFF9", "#FFA0A0"),
          group_order
        ) 
      )+
      stat_compare_means(
        aes(group = Group),                     
        label = "p.signif",                     
        method = "wilcox.test",                     
        hide.ns = F
      )+
      labs(
        title = ct,
        y = "Proportion",
        fill = "Group"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(face = "bold"),
        plot.margin = margin(b = 10)
      )
    sub_plots_list[[ct]] <- p
  }
  p4 <- wrap_plots(sub_plots_list, ncol = 3)
  ggsave(file.path(outpath, paste0(prefix, sub_celltype_column, "_DE.png")), p4, width = sub_DE_plot_width, heigh = sub_DE_plot_height, dpi = 300 ,bg = "white")
  ggsave(file.path(outpath, paste0(prefix, sub_celltype_column, "_DE.pdf")), p4, width = sub_DE_plot_width, heigh = sub_DE_plot_height)
  plot_n <- plot_n + 1
  cat("\nPlot ", plot_n, " : ", paste0(prefix, sub_celltype_column, "_DE.png (pdf)"), " has been generated successfully !\n", sep = "")
}
message("\n=========================All Done !============================\n")





