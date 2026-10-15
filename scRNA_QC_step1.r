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
    A script for per-sample QC of single-cell transcriptomic data (based on Seurat v5).
    Each sample is read independently; QC cutoffs are derived per sample from its own
    quantiles (qL/qH on nFeature_RNA and nCount_RNA, plus a hard cap on percent.mt);
    then QC-passing samples are combined into one Seurat object. Gene-level min.cells
    filtering is applied AFTER combining so that low-frequency genes shared across
    samples are not dropped prematurely.
  Example: scRNA_QC_step1 -I input_list.csv -DT 10X -GC 2 -qL 0.01 -qH 0.99 -MaxM 10 -MC 3 -P 01. -O 01_QC
  '
  )
)
#Required
parser$add_argument(
  "-I",
  "--input",
  dest = "input_file",
  type = "character",
  required = TRUE,
  help = 'Path to a CSV file containing three columns: "path" (file path), "sample" (sample name), and "group" (group name). (Required)'
)
#Optional
parser$add_argument(
  '-DT',
  '--data_type',
  dest = 'data_type',
  help = 'The type of the input data, one of "10X" or "h5". (Optional, default: 10X).',
  type = "character",
  default = "10X"
)
parser$add_argument(
  '-GC',
  '--gene_column',
  dest = 'gene_column',
  help = 'Column number in genes.tsv/features.tsv containing gene names for 10X data (Optional, default: 2).',
  type = "numeric",
  default = 2
)
parser$add_argument(
  '-qL',
  '--quantile_low',
  dest = 'q_low',
  help = 'Per-sample lower quantile for nFeature_RNA / nCount_RNA cutoffs. (Optional, default: 0.01)',
  type = "numeric",
  default = 0.01
)
parser$add_argument(
  '-qH',
  '--quantile_high',
  dest = 'q_high',
  help = 'Per-sample upper quantile for nFeature_RNA / nCount_RNA cutoffs. (Optional, default: 0.99)',
  type = "numeric",
  default = 0.99
)
parser$add_argument(
  '-MaxM',
  '--max_mt',
  dest = 'max_mt',
  help = 'Hard upper bound for percent.mt. (Optional, default: 10)',
  type = "numeric",
  default = 10
)
parser$add_argument(
  '-S',
  '--species',
  dest = 'species',
  help = 'Species for RBC (red blood cell) gene filtering. Must be one of: human, mouse, mammal. Must be used together with -MaxHB. (Optional)',
  type = "character",
  choices = c("human", "mouse", "mammal")
)
parser$add_argument(
  '-MaxHB',
  '--max_HBs',
  dest = 'max_HBs',
  help = 'Cells with hemoglobin genes percentage higher than this value will be discarded. Must be used together with -S. (Optional)',
  type = "numeric"
)
parser$add_argument(
  '-MC',
  '--min_cells',
  dest = 'min_cells',
  help = 'Filter genes to those detected in at least this many cells across ALL samples after merging. (Optional, default: 3)',
  type = "integer",
  default = 3
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
parser$add_argument(
  '-W',
  '--width',
  dest = 'plot_width',
  help = 'width of QC plots in inches (Optional, default: 16)',
  type = "numeric",
  default = 16
)
parser$add_argument(
  '-H',
  '--height',
  dest = 'plot_height',
  help = 'height of QC plots in inches (Optional, default: 16)',
  type = "numeric",
  default = 16
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)
}

opts <- parser$parse_args()
input_file <- opts$input_file
data_type <- opts$data_type
gene_column <- opts$gene_column
q_low <- opts$q_low
q_high <- opts$q_high
max_mt <- opts$max_mt
species <- opts$species
max_HBs <- opts$max_HBs
min_cells <- opts$min_cells
prefix <- opts$prefix
outpath <- opts$outpath
plot_width <- opts$plot_width
plot_height <- opts$plot_height

if (q_low < 0 || q_low > 1 || q_high < 0 || q_high > 1 || q_low >= q_high) {
  stop("Error: -qL / -qH must be within [0,1] and -qL must be smaller than -qH.")
}

# 检查 -MaxHB 不能单独出现（物种信息不可缺失）
if (is.null(species) & !is.null(max_HBs)) {
  stop("Error: '-MaxHB' requires '-S' to be specified. Please provide the species.")
}
# 检查 -S 不能单独出现（必须有 -MaxHB 阈值）
if (!is.null(species) & is.null(max_HBs)) {
  stop("Error: '-S' requires '-MaxHB' to be specified. Please provide the hemoglobin percentage cutoff.")
}
do_rbc <- !is.null(species) & !is.null(max_HBs)

message("\n=======================scRNA Per-sample QC========================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("Copyright © 2026 OYJH. All Rights Reserved.")
message("\n==================================================================\n")

# 设置输出目录
if (!dir.exists(outpath)) {
  dir.create(outpath, recursive = TRUE)
}

# Summary 文件（终端输出的结构化记录）
summary_file <- file.path(outpath, paste0(prefix, "QC_summary.txt"))
summary_con <- file(summary_file, open = "wt")
write_summary <- function(...) {
  # 检测：任一入参是长度 >1 的向量时会被 paste0 collapse 压成一行。
  # 多行表格请改用 cat(..., sep="\n") 直接写 summary_con。
  args_list <- list(...)
  multi <- vapply(args_list, length, integer(1)) > 1
  if (any(multi)) {
    which_arg <- which(multi)[1]
    stop("write_summary arg #", which_arg,
         " is a length-", length(args_list[[which_arg]]),
         " vector which would be collapsed into a single line by paste0. ",
         "Use cat(..., sep=\"\\n\", file = summary_con) instead for multi-line content.")
  }
  txt <- paste0(..., collapse = "")
  cat(txt, file = summary_con)
}
write_summary("===================scRNA Per-sample QC Summary====================\n")
write_summary("Current time: ", as.character(Sys.time()), "\n")
write_summary("Input CSV  : ", input_file, "\n")
write_summary("Data type  : ", data_type, " | gene.column: ", gene_column, "\n")
write_summary("QC quantiles: nFeature/nCount [", q_low, ", ", q_high, "] | percent.mt <=", max_mt,
              if (do_rbc) paste0(" | percent.HB <=", max_HBs, " (species: ", species, ")") else "",
              " | min.cells(after QC): ", min_cells, "\n")

# 2. Load required packages
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(Seurat,quietly = TRUE)){install.packages("Seurat")}
  if(!require(patchwork,quietly = TRUE)){install.packages("patchwork")}
})

# 3. Load input list & validate
input_list <- read.csv(input_file, header = TRUE)
required_cols <- c("path", "sample", "group")
missing_cols <- setdiff(required_cols, colnames(input_list))
if (length(missing_cols) > 0) {
  stop("Error: input CSV is missing required column(s): ", paste(missing_cols, collapse = ", "),
       ". The CSV must contain three columns: path, sample, group.")
}
input_list$path   <- trimws(as.character(input_list$path))
input_list$sample <- trimws(as.character(input_list$sample))
input_list$group  <- trimws(as.character(input_list$group))
if (any(input_list$sample == "")) {
  stop("Error: empty sample name(s) found in the input CSV. Please check the 'sample' column.")
}
if (any(duplicated(input_list$sample))) {
  stop("Error: duplicated sample names found in the input CSV: ",
       paste(unique(input_list$sample[duplicated(input_list$sample)]), collapse = ", "))
}
if (data_type == "10X") {
  missing_paths <- input_list$path[!dir.exists(input_list$path)]
  if (length(missing_paths) > 0) {
    stop("Error: the following 10X data directories do not exist:\n",
         paste(missing_paths, collapse = "\n"))
  }
} else if (data_type == "h5") {
  missing_paths <- input_list$path[!file.exists(input_list$path)]
  if (length(missing_paths) > 0) {
    stop("Error: the following h5 files do not exist:\n",
         paste(missing_paths, collapse = "\n"))
  }
}

# 4. Per-sample read + QC
# 定义红细胞基因集（仅当指定 -S 时使用，从原 scRNA_QC_step2.r 移植）
if (do_rbc) {
  human_RBC <- c("HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1", "GYPA")
  mouse_RBC <- c("Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt", "Alas2", "Gypa")
  mammal_RBC <- c("ALAS2", "FECH", "SLC4A1", "ANK1", "SPTA1", "SPTB", "BPGM")
  rbc_genes <- switch(species,
    "human"  = human_RBC,
    "mouse"  = mouse_RBC,
    "mammal" = mammal_RBC
  )
}

seu_list <- list()
threshold_table <- list()
qc_rows_before <- list()   # QC metrics BEFORE filtering (for plots)
qc_rows_after  <- list()   # QC metrics AFTER filtering  (for plots)

for (i in 1:nrow(input_list)) {
  path <- input_list$path[i]
  sample_name <- input_list$sample[i]
  group <- input_list$group[i]
  cat("\nProcessing sample: ", sample_name, " ......\n", sep = "")
  if (data_type == "10X"){
    data <- Read10X(data.dir = path, gene.column = gene_column)
  }else if (data_type == "h5"){
    data <- Read10X_h5(path)
  }else{
    stop("Error: the input file format cannot be recognized. Only 10X and h5 format files are supported. Please check the file path.")
  }
  seu <- CreateSeuratObject(counts = data, project = sample_name,
                            min.features = 0, min.cells = 0)
  seu$sample <- sample_name
  seu$group  <- group
  # percent.mt: (?i)^mt- covers MT- (human) / mt- (mouse) naming
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "(?i)^mt-")

  # percent.HB: 红细胞基因（仅当指定 -S 时）
  if (do_rbc) {
    rbc_genes_matched <- CaseMatch(rbc_genes, rownames(seu))
    if (length(rbc_genes_matched) == 0) {
      cat("    Warning: None of the RBC genes were found in the data. 'percent.HB' will be set to 0 for all cells.\n")
      seu[["percent.HB"]] <- 0
    } else {
      message("    Species: ", species, " | RBC genes matched: ", paste(rbc_genes_matched, collapse = ", "))
      seu[["percent.HB"]] <- PercentageFeatureSet(seu, features = rbc_genes_matched)
    }
  }

  n_in <- ncol(seu)
  qc_rows_before[[i]] <- data.frame(sample = sample_name,
                                    nFeature_RNA = seu$nFeature_RNA,
                                    nCount_RNA   = seu$nCount_RNA,
                                    percent.mt   = seu$percent.mt,
                                    percent.HB   = if (do_rbc) seu$percent.HB else NA)

  if (n_in < 10) {
    # Too few cells to estimate quantiles — keep all, only mt cap applies
    threshold_table[[i]] <- data.frame(sample = sample_name, n_in = n_in, n_keep = n_in,
                                       ng_lo = NA, ng_hi = NA, tc_lo = NA, tc_hi = NA,
                                       mt_max = max_mt, hb_max = if (do_rbc) max_HBs else NA)
    seu_list[[sample_name]] <- seu
    cat("    [", i, "/", nrow(input_list), "] ", sample_name, ": n=", n_in,
        " (too few cells, kept all)\n", sep = "")
    qc_rows_after[[i]] <- data.frame(sample = sample_name,
                                     nFeature_RNA = seu$nFeature_RNA,
                                     nCount_RNA   = seu$nCount_RNA,
                                     percent.mt   = seu$percent.mt,
                                     percent.HB   = if (do_rbc) seu$percent.HB else NA)
    next
  }
  ng_q <- quantile(seu$nFeature_RNA, probs = c(q_low, q_high))
  tc_q <- quantile(seu$nCount_RNA,   probs = c(q_low, q_high))
  ng_lo <- as.numeric(ng_q[1]); ng_hi <- as.numeric(ng_q[2])
  tc_lo <- as.numeric(tc_q[1]); tc_hi <- as.numeric(tc_q[2])

  keep_mask <- seu$nFeature_RNA >= ng_lo & seu$nFeature_RNA <= ng_hi &
               seu$nCount_RNA   >= tc_lo & seu$nCount_RNA   <= tc_hi &
               seu$percent.mt   <= max_mt
  if (do_rbc) {
    keep_mask <- keep_mask & seu$percent.HB <= max_HBs
  }
  n_keep <- sum(keep_mask)
  threshold_table[[i]] <- data.frame(sample = sample_name, n_in = n_in, n_keep = n_keep,
                                     ng_lo = ng_lo, ng_hi = ng_hi,
                                     tc_lo = tc_lo, tc_hi = tc_hi,
                                     mt_max = max_mt, hb_max = if (do_rbc) max_HBs else NA)
  cat("    [", i, "/", nrow(input_list), "] ", sample_name, ": n=", n_in, " → ", n_keep,
      " (ng=[", ng_lo, ",", ng_hi, "], tc=[", tc_lo, ",", tc_hi, "], mt<=", max_mt,
      if (do_rbc) paste0(", hb<=", max_HBs) else "",
      ")\n", sep = "")
  seu_filtered <- subset(seu, cells = colnames(seu)[keep_mask])
  seu_list[[sample_name]] <- seu_filtered
  qc_rows_after[[i]] <- data.frame(sample = sample_name,
                                   nFeature_RNA = seu_filtered$nFeature_RNA,
                                   nCount_RNA   = seu_filtered$nCount_RNA,
                                   percent.mt   = seu_filtered$percent.mt,
                                   percent.HB   = if (do_rbc) seu_filtered$percent.HB else NA)
}

# 5. Print per-sample threshold table
cat("\nPer-sample QC thresholds:\n")
tt <- do.call(rbind, threshold_table)
rownames(tt) <- NULL
print(tt, row.names = FALSE)
write_summary("\n---------- Per-sample QC thresholds ----------\n")
cat(capture.output(write.table(tt, sep = "\t", row.names = FALSE, quote = FALSE)),
    file = summary_con, sep = "\n")
cat("\n", file = summary_con)

cat("\nNumber of cells per sample before QC:\n")
tb_before <- table(unlist(lapply(qc_rows_before, function(x) rep(as.character(x$sample[1]), nrow(x)))))
print(tb_before)
cat("\nNumber of cells per sample after QC:\n")
tb_after <- table(unlist(lapply(qc_rows_after, function(x) rep(as.character(x$sample[1]), nrow(x)))))
print(tb_after)
write_summary("---------- Cells per sample: before / after QC ----------\n")
all_samp <- names(seu_list)
# 直接从 qc_rows 统计，避免 0 细胞样本在 table() 中消失导致对齐错位
cnt_before2 <- unlist(lapply(qc_rows_before, nrow))
names(cnt_before2) <- unlist(lapply(qc_rows_before, function(x) as.character(x$sample[1])))
cnt_after2  <- unlist(lapply(qc_rows_after, nrow))
names(cnt_after2)  <- unlist(lapply(qc_rows_after,  function(x) as.character(x$sample[1])))
for (s in all_samp) {
  write_summary(s, "\tbefore=", ifelse(s %in% names(cnt_before2), cnt_before2[s], 0),
                "\tafter=",  ifelse(s %in% names(cnt_after2),  cnt_after2[s],  0), "\n")
}
write_summary("\n")

# 6. Plot QC metrics before/after (violin + scatter, 2x2 grid when 4 panels)
plot_qc <- function(qc_rows, outpath, prefix, stage) {
  df <- do.call(rbind, qc_rows)
  df$sample <- factor(df$sample, levels = unique(df$sample))

  # violin: nFeature_RNA / nCount_RNA / percent.mt (+ percent.HB when -S given)
  vln_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  if (do_rbc) vln_features <- c(vln_features, "percent.HB")

  # scatter pairs: (x, y) combinations covering all QC metrics
  sc_pairs <- list(c("nCount_RNA", "nFeature_RNA"),
                   c("nCount_RNA", "percent.mt"),
                   c("nFeature_RNA", "percent.mt"))
  if (do_rbc) sc_pairs <- c(sc_pairs, list(c("nCount_RNA", "percent.HB")))

  gg_list <- lapply(vln_features, function(f) {
    ggplot(df, aes(x = sample, y = .data[[f]], fill = sample)) +
      geom_violin(scale = "width") + theme_bw() + ylim(0, NA) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(title = f, x = "", y = "")
  })
  combo_vln <- Reduce(`+`, gg_list) + plot_layout(ncol = 2)
  ggsave(file.path(outpath, paste0(prefix, stage, "_QC_vlnplot.png")), combo_vln,
         width = plot_width, height = plot_height, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, stage, "_QC_vlnplot.pdf")), combo_vln,
         width = plot_width, height = plot_height)
  cat("Plot : ", prefix, stage, "_QC_vlnplot.png(pdf) has been successfully generated !\n", sep = "")

  sc_list <- lapply(sc_pairs, function(pr) {
    ggplot(df, aes(x = .data[[pr[1]]], y = .data[[pr[2]]], color = sample)) +
      geom_point(size = 0.3, alpha = 0.5) + theme_bw() +
      labs(x = pr[1], y = pr[2])
  })
  combo_sc <- Reduce(`+`, sc_list) + plot_layout(ncol = 2, guides = "collect")
  ggsave(file.path(outpath, paste0(prefix, stage, "_QC_scatter.png")), combo_sc,
         width = plot_width, height = plot_height, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix, stage, "_QC_scatter.pdf")), combo_sc,
         width = plot_width, height = plot_height)
  cat("Plot : ", prefix, stage, "_QC_scatter.png(pdf) has been successfully generated !\n", sep = "")
}

plot_qc(qc_rows_before, outpath, prefix, stage = "before")
plot_qc(qc_rows_after,  outpath, prefix, stage = "after")

# 7. Combine QC-passed samples
cat("\nCombining ", length(seu_list), " QC-filtered samples ...\n", sep = "")
if (length(seu_list) == 1) {
  seu_merged <- seu_list[[1]]
  cat("Only one sample detected, skipping combination (no add.cell.ids prefix added).\n")
} else {
  seu_merged <- merge(seu_list[[1]], y = seu_list[-1],
                      add.cell.ids = names(seu_list))
}
seu_merged <- JoinLayers(seu_merged)

# 8. Gene-level filter (min.cells) uniformly across combined samples
if (min_cells > 0) {
  counts_all <- GetAssayData(seu_merged, layer = "counts")
  keep_genes <- Matrix::rowSums(counts_all > 0) >= min_cells
  n_drop <- sum(!keep_genes)
  if (n_drop > 0) {
    cat("\nFiltering genes detected in fewer than ", min_cells,
        " cells across all samples: removing ", n_drop, " genes.\n", sep = "")
  }
  seu_merged <- subset(seu_merged, features = rownames(seu_merged)[keep_genes])
}

cat("\nTotal cells after QC: ", ncol(seu_merged),
    " | Total genes: ", nrow(seu_merged), "\n", sep = "")
cat("\nNumber of cells per sample after QC:\n")
tb_merged <- table(seu_merged$sample)
print(tb_merged)
write_summary("\n---------- Combined object ----------\n")
write_summary("Total cells: ", ncol(seu_merged), " | Total genes: ", nrow(seu_merged), "\n")
write_summary("Cells per sample after QC:\n")
for (s in names(tb_merged)) {
  write_summary(s, "\t", tb_merged[s], "\n")
}
write_summary("\n")

saveRDS(seu_merged, file.path(outpath, paste0(prefix, "seu_afterQC.rds")))
cat("\nFile 1 : ", prefix, "seu_afterQC.rds saved successfully !\n", sep = "")
close(summary_con)
cat("\nFile 2 : ", prefix, "QC_summary.txt saved successfully !\n", sep = "")
cat("\n==========================All Done!===============================\n")
