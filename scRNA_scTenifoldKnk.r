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
    A script for scTenifoldKnk-based single-cell in silico gene knockout analysis. Example: 
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
    help = 'Path to an rds file for a seurat object used for scRNA scTenifoldKnk analysis. (Required)'
)
parser$add_argument(
    "-TG",
    "--target_gene",
    dest = "target_gene",
    type = "character",
    required = TRUE,
    help = 'The target gene for scTenifoldKnk analysis. (Required)'
)
parser$add_argument(
    "-T",
    "--threads",
    dest = "threads",
    type = "integer",
    default = 8,
    help = "Number of threads to use for scTenifoldKnk analysis. (Default: 8)"
)
parser$add_argument(
    "-p_val",
    dest = "p_threshold",
    type = "double",
    default = 0.05,
    help = "Adjusted p-value threshold for significant differentially regulated genes. (Default: 0.05)"
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
target_gene <- opts$target_gene
threads <- opts$threads
prefix <- opts$prefix
outpath <- opts$outpath
p_threshold <- opts$p_threshold
file_n <- 0
plot_n <- 0

message("\n================ scRNA scTenifoldKnk Analysis ====================\n")
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
    if(!require(scTenifoldKnk,quietly = TRUE)){remotes::install_github('cailab-tamu/scTenifoldKnk')}
    if(!require(enrichR,quietly = TRUE)){install.packages("enrichR")}
    if(!require(igraph,quietly = TRUE)){install.packages("igraph")}
    if(!require(ggrepel,quietly = TRUE)){install.packages("ggrepel")}
})
setEnrichrSite("https://maayanlab.cloud/Enrichr")
set.seed(123)

# 3. Load data (prepare count matrix and gene list for scTenifoldKnk)
cat("Loading data from ", input_rds, " ...\n")
seu <- readRDS(input_rds)
count_matrix <- GetAssayData(seu, layer = "counts")
if (!target_gene %in% rownames(count_matrix)) {
    stop(paste("Error: Target gene", target_gene, "not found in the count matrix."))
}
HVGs <- VariableFeatures(seu)
input_data <- as.data.frame(count_matrix[unique(c(target_gene, HVGs)),])
cat("Input data prepared with ", nrow(input_data), " genes and ", ncol(input_data), " cells.\n")

# 4. Run scTenifoldKnk
cat("\nRunning scTenifoldKnk analysis for target gene: ", target_gene, " ...\n")
ko_analysis_result <- scTenifoldKnk(
    countMatrix = input_data,
    gKO = target_gene,      #需要敲除的基因
    qc = TRUE,              #是否进行QC
    qc_mtThreshold = 0.1,   #mt的阈值
    qc_minLSize = 300,      #文库阈值(细胞测到的基因总数)
    nc_nNet = 10,           #子网络数量
    nc_nCells = 500,        #每个网络中随机抽取的细胞数
    nCores = threads        #并行计算的线程数
)
diff_regulation_df <- ko_analysis_result$diffRegulation
write.csv(diff_regulation_df, file=file.path(outpath, paste0(prefix, "diffRegulation_raw.csv")), quote=F, row.names=F)
file_n <- file_n + 1
cat("\nFile ", file_n, ": ", paste0(prefix, "diffRegulation_raw.csv"), " saved successfully !\n", sep="")

# 保存完整分析结果，便于之后换 q 值重画 KO_effect_plot 而不用重跑
saveRDS(ko_analysis_result, file = file.path(outpath, paste0(prefix, "scTenifoldKnk_result.rds")))
file_n <- file_n + 1
cat("\nFile ", file_n, ": ", paste0(prefix, "scTenifoldKnk_result.rds"), " saved successfully !\n", sep="")
diff_regulation_df <- diff_regulation_df[diff_regulation_df$gene != target_gene, ]
significant_diff_table <- diff_regulation_df[diff_regulation_df$p.adj < p_threshold,]
if (nrow(significant_diff_table) == 0) {
    stop(paste("No significantly differentially regulated genes found with adjusted p-value <", p_threshold, ".\n"))
} else {
    cat(nrow(significant_diff_table), " significantly differentially regulated genes found with adjusted p-value <", p_threshold, ".\n")
    write.csv(significant_diff_table, file=file.path(outpath, paste0(prefix, "diffRegulation_sig.csv")), quote=F, row.names=F)
    file_n <- file_n + 1
    cat("\nFile ", file_n, ": ", paste0(prefix, "diffRegulation_sig.csv"), " saved successfully !\n", sep="")
}

# 5. Visualization
# 5.1 富集分析网络图
png(file.path(outpath, paste0(prefix, "KO_effect_plot.png")), width=16, height=12, res=300, units="in")
plotKO(
ko_analysis_result,
gKO = target_gene,
q = 0.99,
annotate = TRUE
)
dev.off()
pdf(file.path(outpath, paste0(prefix, "KO_effect_plot.pdf")), width=16, height=12)
plotKO(
ko_analysis_result,
gKO = target_gene,
q = 0.99,
annotate = TRUE
)
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, ": ", paste0(prefix, "KO_effect_plot.png (pdf)"), " has been generated successfully !\n", sep="")

# 5.2 FC值条形图
top20_diff_genes <- significant_diff_table %>% 
    arrange(desc(FC)) %>% 
    head(20)
p1 <- ggplot(top20_diff_genes, aes(x=reorder(gene, FC), y=FC)) +
  geom_bar(stat='identity', fill='#5A9BD4') +
  coord_flip() + 
  scale_y_continuous(expand = c(0, 0.5)) +
  labs(title="Top 20 Differentially Regulated Genes", x="Gene", y="FC") +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size=12, face="bold"),
        axis.text.x = element_text(size=12, face="bold"),
        axis.title = element_text(size=14, face="bold")
  )
ggsave(filename=file.path(outpath, paste0(prefix, "diff_genes_barplot.png")), plot=p1, width=12, height=10, dpi=300, bg="white")
ggsave(filename=file.path(outpath, paste0(prefix, "diff_genes_barplot.pdf")), plot=p1, width=12, height=10)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, ": ", paste0(prefix, "diff_genes_barplot.png (pdf)"), " has been generated successfully !\n", sep="")

# 5.3 Z-score与-log10(p-value)散点图
plot_df <- diff_regulation_df
plot_df$log_pval <- -log10(diff_regulation_df$p.adj)
label_genes <- subset(plot_df, abs(Z) > 2 & p.adj < p_threshold)
p2 <- ggplot(plot_df, aes(x=Z, y=log_pval)) +
        geom_point(alpha=0.5) +
        geom_hline(yintercept=-log10(p_threshold), linetype="dashed", color="red") +
        geom_vline(xintercept=c(2), linetype="dashed", color="blue") +
        geom_text_repel(data=label_genes, aes(label=gene),
                  size=4, max.overlaps=50) +
        labs(title="Z vs -log10(p-value)",
            x="Z-score", y="-log10(p-value)") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size=12, face="bold"),
            axis.title = element_text(size=14, face="bold")
)
ggsave(filename=file.path(outpath, paste0(prefix, "Z_logp_scatter.png")), plot=p2, width=12, height=10, dpi=300, bg="white")
ggsave(filename=file.path(outpath, paste0(prefix, "Z_logp_scatter.pdf")), plot=p2, width=12, height=10)
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, ": ", paste0(prefix, "Z_logp_scatter.png (pdf)"), " has been generated successfully !\n", sep="")

# 5.4 显著基因占比饼图
diff_regulation_df$significant <- ifelse(diff_regulation_df$p.adj < p_threshold, "Significant", "Not significant")
significant_gene_count <- table(diff_regulation_df$significant)
significant_gene_df <- as.data.frame(significant_gene_count)
colnames(significant_gene_df) <- c("category", "count")
significant_gene_df$percentage <- paste0(round(significant_gene_df$count / sum(significant_gene_df$count) * 100, 1), "%")
# 绘制饼图
p3 <- ggplot(significant_gene_df, aes(x="", y=count, fill=category)) +
        geom_bar(stat="identity", width=1) +
        geom_text(aes(label=percentage), position=position_stack(vjust=0.5), size=4) + 
        coord_polar("y", start=0) + 
        scale_fill_manual(values=c("Significant"="red", "Not significant"="lightgray")) + 
        labs(title="Proportion of Significant Genes", fill="") +
        theme_minimal() + 
        theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text = element_blank(),
                panel.grid = element_blank(), legend.position = "right", plot.title = element_text(hjust = 0.5))
# 输出饼图
png(file=file.path(outpath, paste0(prefix, "pie.png")), width=6, height=5, units = "in", res=300)
print(p3)
dev.off()
pdf(file=file.path(outpath, paste0(prefix, "pie.pdf")), width=6, height=5)
print(p3)
dev.off()
plot_n <- plot_n + 1
cat("\nPlot ", plot_n, ": ", paste0(prefix, "pie.png (pdf)"), " has been generated successfully !\n", sep="")

message("\n=========================All Done !============================\n")