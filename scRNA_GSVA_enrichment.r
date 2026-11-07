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
  Example: scRNA_GSVA_enrichment -I seu_annotated.rds -AC celltype -SC orig.ident -GC Group -G WT,KO -S mouse -P test_ -O test
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
    help = 'Path to an rds file for a seurat object used for scRNA enrichment analysis based on GSVA method. (Required)'
)

parser$add_argument(
    "-AC",
    "--annotation_column",
    dest = "annotation_column",
    type = "character",
    required = TRUE,
    help = 'Column name containing annotation information in the Seurat object used for GSVA analysis. (Required)'
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
    "-GC",
    "--group_column",
    dest = "group_column",
    type = "character",
    required = TRUE,
    help = 'Column name containing group information in the Seurat object used for GSVA analysis. (Required)'
)
parser$add_argument(
    "-G", 
    "--groups",
    dest = "groups",
    type = "character",
    help = "Comma-separated group names for comparison. Note: The first group is considered the CONTROL/BASELINE group, and the second group is the TREATMENT/EXPERIMENTAL group. Example: WT,KO for wild-type vs knockout, or Untreated,Treated for untreated vs treated samples. (Required)"
)
parser$add_argument(
    "-S", 
    "--species",
    dest = "species",
    type = "character",
    help = " 'human' or 'mouse'. (Required)"
)


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
group_column <- opts$group_column
groups <- unlist(strsplit(opts$groups, split = ","))
assay <- opts$assay
species <- opts$species
prefix <- opts$prefix
outpath <- opts$outpath
file_n <- 0
plot_n <- 0

message("\n======================scRNA GSVA Enrichment====================\n")
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
  if(!require(GSVA,quietly = TRUE)){BiocManager::install("GSVA")}
  if(!require(msigdbr,quietly = TRUE)){BiocManager::install("msigdbr")}
  if(!require(limma,quietly = TRUE)){BiocManager::install("limma")}
  if(!require(ComplexHeatmap,quietly = TRUE)){BiocManager::install("ComplexHeatmap")}
})

# 3. Load data
seu <- readRDS(input_rds)
if (species == "human") {
    species <- "Homo sapiens"
    db_species = "HS"
    category = "C"
}else if (species == "mouse") {
    species <- "Mus musculus"
    db_species = "MM"
    category = "M"
}else{
    stop("Species must be either 'Human' or 'Mouse'")
}
msig_go <- msigdbr(
  species = species,
  db_species = db_species,
  collection = paste0(category, "5"),
  subcollection = "GO:BP"
)
go_sets <- split(msig_go$gene_symbol, msig_go$gs_name)
msig_kegg <- msigdbr(
  species = species,
  collection = "C2",
  subcollection = "CP:KEGG_MEDICUS"
)
kegg_sets <- split(
  msig_kegg$gene_symbol,
  msig_kegg$gs_name
)
msig_reactome <- msigdbr(
  species = species,
  db_species = db_species,
  collection = paste0(category, "2"),
  subcollection = "CP:REACTOME"
)
reactome_sets <- msig_reactome %>%
  split(x = .$gene_symbol, f = .$gs_name)
gene_sets <- c(go_sets, kegg_sets, reactome_sets)
gene_sets <- lapply(gene_sets, unique)

# 4. GSVA score calculation per celltype pseudo-bulk
pb <- AggregateExpression(
  seu,
  assays   = assay,
  slot     = "counts",
  group.by = c(sample_column, annotation_column)
)
pb_counts <- pb[[assay]]
pb_norm <- log2(
  t(t(pb_counts) / colSums(pb_counts)) * 1e6 + 1
)
celltypes <- unique(sub(".*_", "", colnames(pb_norm)))
gsva_results <- list()
for (ct in celltypes) {
  cat("\nRunning GSVA for celltype: ", ct, " ...\n", sep = "")
  # 找到该 celltype 在 pseudo-bulk 矩阵中的列
  ct_cols <- colnames(pb_norm)[endsWith(colnames(pb_norm), paste0("_", ct))]
  expr_mat <- pb_norm[, ct_cols, drop = FALSE]
  # GSVA 参数（适合 log-CPM pseudo-bulk）
  gsvapar <- gsvaParam(
    expr = as.matrix(expr_mat),
    geneSets = gene_sets,
    kcdf = "Gaussian",
    minSize = 10,
    maxSize = 300
  )
  gsva_res <- gsva(
    gsvapar,
    verbose = FALSE
  )
  gsva_results[[ct]] <- gsva_res
}
df_all <- lapply(names(gsva_results), function(ct) {
    gsva_mat <- gsva_results[[ct]]
    df <- as.data.frame(gsva_mat)
    df <- tibble::rownames_to_column(df, "pathway")
    df$celltype <- ct
    df <- df %>%
    dplyr::select(
    celltype,
    pathway,
    dplyr::everything()
    )
    sample_cols <- setdiff(
    colnames(df),
    c("celltype", "pathway", "geneset")
    )
    sample_cols_new <- sub("_[^_]+$", "", sample_cols)
    sample_cols_new <- gsub("-", "_", sample_cols_new)
    colnames(df)[colnames(df) %in% sample_cols] <- sample_cols_new
    gene_sets <- attr(gsva_mat, "geneSets")
    geneset_df <- data.frame(
    pathway = names(gene_sets),
    geneset = vapply(
    gene_sets,
    function(x) paste(x, collapse = ","),
    character(1)
    )
    )
    df <- dplyr::left_join(df, geneset_df, by = "pathway")
})  %>%  dplyr::bind_rows()
write.csv(df_all, file.path(outpath, paste0(prefix,"GSVA_score.csv")), quote = FALSE, row.names=FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ",prefix,"GSVA_score.csv saved successfully !\n", sep = "")

# 5. DE Analysis between groups per celltype
sample_meta <- seu@meta.data  %>% 
  dplyr::select(
    sample = all_of(sample_column),
    group  = all_of(group_column)
  )  %>% 
  dplyr::distinct()
stopifnot(
  nrow(sample_meta) == length(unique(seu@meta.data[[sample_column]]))
)
deg_results <- list()
for (ct in names(gsva_results)) {
    cat("\nPerforming DE analysis for celltype: ", ct, " ...\n", sep = "")
    gsva_mat <- gsva_results[[ct]]   # pathways × samples
    samples <- colnames(gsva_mat)
    sample_ids <- sub("_[^_]+$", "", samples)
    sample_ids_std <- gsub("-", "_", sample_ids)
    sample_meta$sample_std <- gsub("-", "_", sample_meta$sample)
    group <- sample_meta$group[
        match(sample_ids_std, sample_meta$sample)
    ]
    group <- factor(group, levels = groups)

    # 样本数检查
    tab <- table(group)
    if (any(tab < 3)) {
    message("  Skip ", ct, ": insufficient samples per group")
    next
    }

    # design matrix
    design <- model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    fit <- lmFit(gsva_mat, design)
    contrast_str <- paste0(groups[2], " - ", groups[1])
    contrast     <- makeContrasts(contrasts = contrast_str, levels = design)

    fit2 <- contrasts.fit(fit, contrast)
    fit2 <- eBayes(fit2)

    deg <- topTable(
    fit2,
    number = Inf,
    adjust.method = "BH"
    )

    deg$pathway  <- rownames(deg)
    deg$celltype <- ct
    deg$contrast <- contrast_str
    deg_results[[ct]] <- deg
}

deg_df_all <- do.call(rbind, deg_results)
write.csv(deg_df_all, file.path(outpath, paste0(prefix, "GSVA_DE_result.csv")), quote=FALSE, row.names = FALSE)
file_n <- file_n + 1
cat("\nFile ", file_n, " : ", prefix, "GSVA_DE_result.csv saved successfully !\n", sep = "")

# 6. Heatmap visualization per celltype
get_db <- function(pathway) {
  sub("_.*", "", pathway)
}
db_levels <- c("GOBP", "KEGG", "REACTOME")
select_pathways <- function(df, top_n = 5) {

  df <- df %>%
    mutate(
      database  = get_db(pathway),
      database  = factor(database, levels = db_levels),
      direction = ifelse(logFC >= 0, "Up", "Down"),
      direction = factor(direction, levels = c("Up", "Down")),
      abs_logFC = abs(logFC)
    )

  ## 1️⃣ 显著性筛选逻辑（保持你之前的规则）
  if (any(df$adj.P.Val < 0.05, na.rm = TRUE)) {
    message("✔ Using adj.P.Val < 0.05")
    df_filt <- df %>% filter(adj.P.Val < 0.05)

  } else if (any(df$P.Value < 0.05, na.rm = TRUE)) {
    message("⚠ No adj.P.Val < 0.05, fallback to P.Value < 0.05")
    df_filt <- df %>% filter(P.Value < 0.05)

  } else {
    message("⚠ No significant pathways, fallback to |logFC| ranking")
    df_filt <- df
  }

  ## 2️⃣ 每个数据库：Up 取前5，Down 取前5
  df_top <- df_filt %>%
    group_by(database, direction) %>%
    arrange(desc(abs_logFC)) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    arrange(database, direction, desc(abs_logFC))

  return(df_top)
}
db_colors <- c(
  GOBP     = "#66c2a5",  # 浅绿色
  KEGG     = "#fc8d62",  # 浅橙色
  REACTOME = "#ffd92f"   # 浅黄色
)

celltypes <- names(gsva_results)

for (ct in celltypes) {
    cat("\nGenerating heatmap for celltype: ", ct, " ...\n", sep = "")

    ## 4.1 差异分析结果
    df_deg <- deg_results[[ct]]
    if (is.null(df_deg) || nrow(df_deg) == 0) {
    message("⚠ No DEG results for ", ct, ", skipped")
    next
    }

    ## 4.2 通路筛选
    df_sel <- select_pathways(df_deg, top_n = 10)
    if (nrow(df_sel) == 0) {
    message("⚠ No pathway selected for ", ct, ", skipped")
    next
    }

    ## 4.3 提取 GSVA 矩阵
    gsva_mat <- gsva_results[[ct]]
    pathways <- df_sel$pathway

    common_pw <- intersect(pathways, rownames(gsva_mat))
    if (length(common_pw) == 0) {
    message("⚠ No matching pathways in GSVA matrix for ", ct)
    next
    }

    mat <- gsva_mat[common_pw, , drop = FALSE]

    ## 4.4 行顺序：数据库 → |logFC|
    row_info <- df_sel %>%
    filter(pathway %in% common_pw) %>%
    arrange(database, direction, desc(abs_logFC))

    row_order <- row_info$pathway
    row_split <- row_info$database

    mat <- mat[row_order, , drop = FALSE]

    ## 4.5 列顺序：Control 在左，Treatment 在右
    # 提取 sample id（去掉 _celltype 后缀）
    sample_ids <- sub(paste0("_", ct, "$"), "", colnames(mat))

    sample_meta <- unique(seu@meta.data[, c(sample_column, group_column)])
    colnames(sample_meta) <- c("sample", "group")
    sample_meta$sample <- gsub("_", "-", sample_meta$sample)

    group_vec <- sample_meta$group[
    match(sample_ids, sample_meta$sample)
    ]

    # 控制顺序：groups[1] 在左，groups[2] 在右
    col_order <- order(factor(group_vec, levels = groups))

    mat <- mat[, col_order, drop = FALSE]
    column_labels <- sub(paste0("_", ct, "$"), "", colnames(mat))

    ## 4.6 行注释（数据库）
    row_db <- get_db(rownames(mat))
    ha_row <- rowAnnotation(
    Database = row_db,
    col = list(Database = db_colors),
    show_annotation_name = FALSE,
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold")
    )
    ## 4.6 绘图
    ht <- Heatmap(
        mat,
        name = "GSVA",
        right_annotation = ha_row,
        row_split = row_split,
        row_gap = unit(3, "mm"),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_order = row_order,
        column_order = colnames(mat),
        row_names_side = "right",
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 16, fontface = "bold"),
        column_names_rot = 45,
        row_title = NULL,
        column_title = NULL,
        column_labels = column_labels,  # <- 这里指定显示的列名
        width = unit(14, "in"),
        height = unit(16, "in"),
        rect_gp = gpar(
            col = "black",   # 边框颜色
            lwd = 0.8        # 边框粗细（0.5–1 推荐）
        )
)

    png(
        filename = file.path(outpath, paste0(prefix, "GSVA_heatmap_", ct, ".png")),
        width = 30,
        height = 20,
        units = "in",
        res = 300
    )
    draw(ht, padding = unit(c(2, 2, 2, 10), "in"), heatmap_legend_side = "left", annotation_legend_side = "top")
    dev.off()
    pdf(
        file = file.path(outpath, paste0(prefix, "GSVA_heatmap_", ct, ".pdf")),
        width = 30,
        height = 20
    )
    draw(ht, padding = unit(c(2, 2, 2, 10), "in"), heatmap_legend_side = "left", annotation_legend_side = "top")
    dev.off()
    plot_n <- plot_n + 1
    cat("Plot ", plot_n, " : ", paste0(prefix, "GSVA_heatmap_", ct, ".png (pdf)"), " has been generated successfully !\n", sep = "")
}
message("\n=========================All Done !============================\n")


