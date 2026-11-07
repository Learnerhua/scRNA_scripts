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
A script for enrichment analysis of single-cell data.
Example: scRNA_spRNA_enrichment -I seo_annotated.rds -DER KO_vs_WT_DE_raw.csv -CC celltype -GC Group -G WT,KO -A RNA -S mouse -IC Unknown -P 15w_ -O test
"
))

parser$add_argument(
  "-I",
  "--input_rds",
  help = "Path to an RDS file containing a Seurat object for enrichment analysis (the same RDS file in DE analysis).",
  type = "character",
  required = TRUE
)
parser$add_argument(
  "-DER",
  "--DE_result",
  help = "Path to the raw DE result (csv file).",
  type = "character",
  required = TRUE
)
parser$add_argument(
  "-CC",
  "--celltype_column",
  help = "Column name for celltype identity in metadata of the seurat object provided. (e.g., celltype)",
  type = "character",
  required = TRUE
)
parser$add_argument(
  "-GC",
  "--group_column",
  help = "Column name for group identity in metadata of the seurat object provided. (e.g., Group)",
  type = "character",
  required = TRUE
)
parser$add_argument("-G", "--groups", help = "The compared groups used in DE analysis (comma or space separated, e.g., WT,KO)")
parser$add_argument("-A", "--assay", help = "Expressed assay used in DE analysis (default : RNA)", default = "RNA")
parser$add_argument(
  "-p",
  "--p_threshold",
  help = "Adjusted p value threshold for filtering DEGs used in GO & KEGG analysis (default : 0.05)",
  type = "double",
  default = 0.05
)
parser$add_argument(
  "-log2FC",
  "--log2FC_threshold",
  help = "Log2FC threshold for filtering DEGs used in GO & KEGG analysis (default : 0)",
  type = "double",
  default = 0
)
parser$add_argument(
  "-GKP",
  "--GO_KEGG_p_threshold",
  help = "P value threshold in GO & KEGG analysis (default : 0.05)",
  type = "double",
  default = 0.05
)
parser$add_argument(
  "-GKQ",
  "--GO_KEGG_q_threshold",
  help = "Q value threshold in GO & KEGG analysis (default : 0.2)",
  type = "double",
  default = 0.2
)
parser$add_argument('-GM',
                    '--GSEA_metric',
                    default = "avg_log2FC",
                    help = 'The ranking metric used to order genes in the GSEA analysis. This should be the name of the column from the file provided by the "-DER" parameter. (default : avg_log2FC)'
)
parser$add_argument(
  "-GSEA_P",
  "--GSEA_p_threshold",
  help = "P value threshold in GSEA analysis (default : 0.05)",
  type = "double",
  default = 0.05
)
parser$add_argument(
  "-GSEA_Q",
  "--GSEA_q_threshold",
  help = "Q value threshold in GSEA analysis (default : 0.25)",
  type = "double",
  default = 0.25
)
parser$add_argument('-NE',
                    '--GSEA_NES',
                    default = 1,
                    type = 'double',
                    help = 'The absolute NES (Normalized Enrichment Score) threshold used to filter GSEA results. Only gene sets with |NES| greater than  this value will be retained. (default : 1).')
parser$add_argument('-S',
                    '--species',
                    default = "human",
                    help = 'one of "human" or "mouse" (default: human).'
)
parser$add_argument('-IC',
                    '--ignore_celltype',
                    help = 'Specific celltype that should be excluded in the enrichment analysis. (e.g., Unknown)'
)
parser$add_argument("-P",
                    "--prefix",
                    help = "Prefix for output files (optional, default: none)",
                    type = "character",
                    default = "")
parser$add_argument(
  "-O",
  "--outpath",
  help = "Output directory (optional, default: current directory)",
  type = "character",
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
DE_result <- opts$DE_result
celltype_column <- opts$celltype_column
group_column <- opts$group_column
groups <- unlist(strsplit(opts$groups, "[, ]+"))
assay <- opts$assay
p_threshold <- opts$p_threshold
log2FC_threshold <- opts$log2FC_threshold
GO_KEGG_p_threshold <- opts$GO_KEGG_p_threshold
GO_KEGG_q_threshold <- opts$GO_KEGG_q_threshold
GSEA_p_threshold <- opts$GSEA_p_threshold
GSEA_q_threshold <- opts$GSEA_q_threshold
GSEA_NES <- opts$GSEA_NES
GSEA_metric <- opts$GSEA_metric
species <- opts$species
prefix <- opts$prefix
outpath <- opts$outpath
if (!is.null(opts$ignore_celltype)){
  ignore_celltype <- unlist(strsplit(opts$ignore_celltype, "[, ]+"))
}else{
  ignore_celltype <- NULL
}

message("\n======================ScRNA Enrichment anaysis==================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2025 OYJH. All Rights Reserved.")
message("\n===============================================================\n")

result_out <- file.path(outpath, "result")
if (!dir.exists(result_out)){
  dir.create(result_out, recursive = TRUE)
}

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
})


# 3. BG genes
seo <- readRDS(input_rds)
meta_data <- seo@meta.data
seo_subset <- seo[, meta_data[[group_column]] %in% groups]
meta_data_subset <- seo_subset@meta.data
celltypes <- unique(meta_data_subset[[celltype_column]])
if (!is.null(ignore_celltype)){
  celltypes <- celltypes[!celltypes %in% ignore_celltype]
}
cat("\nExtracting background genes for each celltype ......\n")
for (ct in celltypes){
  seo_ct <- seo_subset[, meta_data_subset[[celltype_column]] == ct]
  count_mat <- GetAssayData(seo_ct, assay = assay, layer = "counts")
  keep_genes <- rowSums(count_mat > 1) > 0
  genes_to_write <- names(keep_genes)[keep_genes]
  if (length(genes_to_write) < 10){
    cat("Background genes number in celltype", ct, "less than 10, discard it !\n")
  }else{
    cat("Background genes number in celltype", ct,":",length(genes_to_write), "\n")
    write.csv(data.frame(Gene = genes_to_write), file.path(outpath, paste0(prefix, ct, "_bg_genes.csv")), row.names = F, quote = F)
  }
}

# 4. DEGs
df <- read.csv(DE_result)
df_filtered <- df %>% 
  filter(p_val_adj < p_threshold & abs(avg_log2FC) > log2FC_threshold)
cat("\nExtracting  DEGs for each celltype ......\n")
for (ct in celltypes){
  DEGs <- df_filtered %>% 
    filter(cluster == ct) %>% 
    pull(gene)
  if (length(DEGs) < 10){
    cat("DEGs number in celltype", ct, "less than 10, discard it !\n")
  }else{
    cat("DEGs number in celltype", ct,":",length(DEGs), "\n")
    writeLines(DEGs, file.path(outpath, paste0(prefix, ct, "_DEGs.txt")))}
}

# 5. GSEA
cat("\nPreparing GSEA input file ......\n")
for (ct in celltypes){
  df_subset <- df %>% 
    filter(cluster == ct)
  rownames(df_subset) <- df_subset$gene
  df_subset <- df_subset[, setdiff(colnames(df_subset), "gene")]
  cat("GSEA input gene number for celltype", ct, ":", nrow(df_subset), "...\n")
  write.csv(df_subset, file = file.path(outpath, paste0(prefix, ct, "_GSEA_input.csv")), row.names = T, quote = F)
}

# 6.Run GO & KEGG analysis
cat("\nStarting GO & KEGG analysis ......\n")
for (ct in celltypes){
  bg_file <- file.path(outpath, paste0(prefix, ct, "_bg_genes.csv"))
  gene_file <- file.path(outpath, paste0(prefix, ct, "_DEGs.txt"))
  if (file.exists(bg_file) && file.exists(gene_file)) {
    cat("Processing celltype", ct, "...\n")
    cmd <- paste(
      "enrichment_analysis",
      "-G", shQuote(gene_file),
      "-E", shQuote(bg_file),
      "-S", species,
      "-P", shQuote(paste0(prefix, ct, "_")),
      "-p", GO_KEGG_p_threshold,
      "-Q", GO_KEGG_q_threshold,
      "-O", shQuote(file.path(outpath, "result")),
      ">", shQuote(file.path(outpath, paste0(prefix, ct, "_run.log"))),
      "2>&1 &"
    )
    # Save command to shell script for manual rerun
    sh_file <- file.path(outpath, paste0(prefix, ct, "_run.sh"))
    writeLines(cmd, sh_file)
    cat("Script saved to:", sh_file, "\n")
    # Execute command in background
    system(cmd)
  }
}

# 7. GSEA analysis
cat("\nStarting GSEA analysis ......\n")
for (ct in celltypes){
  cat("Executing GSEA analysis for", ct, "...\n")
  cmd <- paste(
    "enrichment_analysis",
    "-GSEA", shQuote(file.path(outpath, paste0(prefix, ct, "_GSEA_input.csv"))),
    "-S", species,
    "-P", shQuote(paste0(prefix, ct, "_")),
    "-NE", GSEA_NES,
    "-GM", GSEA_metric,
    "-p", GSEA_p_threshold,
    "-Q", GSEA_q_threshold,
    "-O", shQuote(result_out),
    ">", shQuote(file.path(outpath, paste0(prefix, ct, "_GSEA_run.log"))),
    "2>&1 &"
  )
  # Save command to shell script for manual rerun
  sh_file <- file.path(outpath, paste0(prefix, ct, "_GSEA_run.sh"))
  writeLines(cmd, sh_file)
  cat("Script saved to:", sh_file, "\n")
  # Execute command in background
  system(cmd)
}

message("\n=========================All Done !============================\n")



