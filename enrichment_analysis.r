#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/Rscript
rm(list = ls())
options(timeout = 1000)

# 1.参数设置
suppressPackageStartupMessages({
  if (!require(argparse, quietly = TRUE)) {
    install.packages("argparse")
  }
})
parser <- ArgumentParser(
  description = paste0(
    '
  Script for enrichment analysis for RNA data.
  Example: enrichment_analysis -G gene.txt -E expr_data.csv -P 02. -p 0.05 -Q 0.2 -O /home/user/results
  '
  )
)
parser$add_argument('-G',
                    '--genes',
                    dest = 'gene_file',
                    help = 'Target genes file used for enrichment analysis (GO & KEGG). Format: one gene per line (e.g., "TP53\nBRCA1\nEGFR"). Minimum requirement: At least 10 genes. Conflict with parameter "-GSEA".')
parser$add_argument(
  '-E',
  '--expr',
  dest = 'expr_data',
  help = 'expression file used to extract background genes. Format : a file with row names as gene symbols.
                    If not provided, all genes in the database will be used as background by default.'
)
parser$add_argument('-S',
                    '--species',
                    dest = 'species',
                    help = 'one of "human" or "mouse" (default: human).',
                    default = "human")
parser$add_argument(
  '-p',
  '--pvalue',
  dest = 'p_threshold',
  type = "double",
  default = 0.05,
  help = 'P-value cutoff for significant terms or pathways (default: 0.05).'
)
parser$add_argument(
  '-Q',
  '--qvalue',
  dest = 'q_threshold',
  type = "double",
  default = 0.2,
  help = 'Q-value (FDR) cutoff for significant terms or pathways (default: 0.2). When performing GSEA analysis, it is recommended to set this value to 0.25.')
parser$add_argument('-GSEA',
                    '--GSEA_file',
                    dest = 'gsea_file',
                    help = 'The input file for GSEA, where the row names are gene symbols. This parameter must be used together with "-GM", and cannot be used together with parameter "-G".')
parser$add_argument('-NE',
                    '--NES_threshold',
                    dest = 'nes_threshold',
                    default = 1,
                    type = 'double',
                    help = 'The absolute NES (Normalized Enrichment Score) threshold used to filter GSEA results. Only gene sets with |NES| greater than  this value will be retained. Default is 1.')
parser$add_argument('-GM',
                    '--GSEA_metric',
                    dest = 'gsea_me',
                    help = 'The ranking metric used to order genes in the GSEA analysis. This should be the name of the column from the file provided by the "-GSEA" parameter. This parameter must be specified together with "-GSEA".'
)
parser$add_argument(
  '-M',
  '--MSigDB',
  dest = 'MSigDB',
  help = 'Use MSigDB KEGG gene sets (version v2025.1) instead of online KEGG for pathway analysis. This option is recommended when network connection to the KEGG database is unstable or unavailable. (default: FALSE)',
  action = 'store_true',
  default = FALSE
)
parser$add_argument(
  '-P',
  '--prefix',
  dest = 'prefix',
  help = 'prefix for output files (optional)',
  type = "character",
  default = ""
)
parser$add_argument(
  '-O',
  '--outpath',
  dest = 'outpath',
  help = 'output directory (default: current directory)',
  default = getwd()
)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  parser$print_help()
  quit(status = 1)  # 退出并返回错误码
}

opts <- parser$parse_args()
genefile <- opts$gene_file
species <- opts$species
internal_data <- opts$internal_data
gsea_file <- opts$gsea_file
gsea_me <- opts$gsea_me
nes_threshold <- opts$nes_threshold
MSigDB <- opts$MSigDB
expr_data <- opts$expr_data
p_threshold <- opts$p_threshold
q_threshold <- opts$q_threshold
prefix <- opts$prefix
outpath <- opts$outpath

if (!is.null(genefile) && !is.null(gsea_file)){
    stop('The parameters "-G" and "-GSEA" cannot be specified simultaneously!')
}
if (!is.null(expr_data) && is.null(genefile)){
  stop('The parameter "-E" cannot be specified without parameter "-G"')
}
if (is.null(gsea_file) != is.null(gsea_me)) {
  stop('Both "-GSEA" and "-GM" parameters must be specified together.')
}
if (species == "human"){
  orgdb <- "org.Hs.eg.db"
  organism <- "hsa"
}else if(species == "mouse"){
  orgdb <- "org.Mm.eg.db"
  organism <- "mmu"
}else{
  stop('Please specify the species parameter using "human" or "mouse" !')
}

message("\n======================Enrichment Analysis======================\n")
message("Current time: ", Sys.time())
message("Operation system: ", R.version$os)
message("Platform: ", R.version$platform)
message("Working directory: ", getwd())
message("Output directory: ", outpath)
message("\n===============================================================\n")

# 设置输出目录
output_dir <- outpath
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. 加载R包
suppressPackageStartupMessages({
  if(!require(tidyverse,quietly = TRUE)){install.packages("tidyverse")}
  if(!require(clusterProfiler,quietly = TRUE)){BiocManager::install("clusterProfiler")}
  if(!require(enrichplot,quietly = TRUE)){install.packages("enrichplot")}
  if(!require(orgdb, quietly = TRUE, character.only = TRUE)){BiocManager::install(orgdb)}
  if(!require(GseaVis,quietly = TRUE)){install.packages("GseaVis")}
})

if (!is.null(genefile)){
  suppressWarnings({
    genes <- readLines(genefile)
  })
}
if (!is.null(expr_data)){
  expr <- read.csv(expr_data, row.names = 1)
  bg_genes <- rownames(expr)
}

# 3 GO分析
if (!is.null(genefile)){
  # 3.1 函数定义
  GO_analysis <- function(ontology) {
    # 根据 bg_gene 是否存在决定是否传递 universe 参数
    if (!is.null(bg_genes)) {
      ego <- enrichGO(
        gene          = genes,
        universe      = bg_genes,      # 使用自定义背景基因
        OrgDb         = orgdb,
        keyType       = "SYMBOL",
        ont           = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_threshold,
        qvalueCutoff  = q_threshold,
        readable      = TRUE
      )
    } else {
      ego <- enrichGO(
        gene          = genes,        # 不传递 universe，使用默认背景
        OrgDb         = orgdb,
        keyType       = "SYMBOL",
        ont           = ontology,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_threshold,
        qvalueCutoff  = q_threshold,
        readable      = TRUE
      )
    }

    # 统一处理结果
    if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
    ego <- dropGO(ego, level = 1:3)
    if (nrow(ego@result) == 0) return(NULL)
    ego <- pairwise_termsim(ego)
    ego <- simplify(ego, cutoff = 0.9, by = "p.adjust", select_fun = min)
    if (nrow(ego@result) == 0) return(NULL)
    return(ego)
  }
  # 3.2 函数运行
  onts <- c("BP", "MF", "CC")
  ego_list <- lapply(onts, GO_analysis)
  names(ego_list) <- onts
  ego_list <- Filter(Negate(is.null), ego_list)
  if (length(ego_list) == 0) {
    warning("No significant terms are enriched in GO analysis !")
  } else{
    combined_ego <- merge_result(ego_list)
    combined_ego_df <- as.data.frame(combined_ego)
    write.csv(combined_ego_df, file.path(outpath, paste0(prefix, "GO_result.csv")), row.names = FALSE)
    cat("\nFile 1 : ", prefix, "GO_result.csv saved successfully !\n", sep = "")
    #可视化
    plot_data <- combined_ego_df %>%
      # 计算 GeneRatio 的实际数值（如 "5/100" → 0.05）
      mutate(GeneRatio = sapply(strsplit(GeneRatio, "/"),
                                function(x) as.numeric(x[1]) / as.numeric(x[2]))) %>%
      # 按 Cluster 分组，并直接选取 Count 最大的前5行
      group_by(Cluster) %>%
      slice_max(Count, n = 5) %>%
      ungroup()  # 取消分组
    p1 <- ggplot(plot_data,
                 aes(x = GeneRatio,
                     y = reorder(Description, GeneRatio),
                     size = Count,
                     color = -log10(p.adjust))) +
      geom_point() +
      facet_grid(Cluster ~ ., scales = "free_y", space = "free_y") +
      scale_color_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
      scale_size_continuous(range = c(4, 10), name = "Gene Count") +
      labs(x = "Gene Ratio", y = "",
           title = "GO Enrichment Dotplot") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold", size = 14),
            panel.grid.major.y = element_line(linetype = "dashed"),
            axis.text = element_text(size = 14, color = "black", face = "bold"),
            axis.title = element_text(size = 16, color = "black", face = "bold"),
            plot.title = element_text(size = 20,  color = "black", face = "bold"),
            legend.title = element_text(size = 12,  color = "black", face = "bold"),
            legend.text = element_text(size = 12,  color = "black", face = "bold")
      )+
      scale_y_discrete(
        labels = function(x) str_wrap(x, width = 50)
      )
    ggsave(file.path(outpath, paste0(prefix,"GO_dotplot.png")), p1, width = 16, height = 12, dpi = 300, bg = "white")
    ggsave(file.path(outpath, paste0(prefix,"GO_dotplot.pdf")), p1, width = 16, height = 12)
    cat("\nPlot 1 : ", prefix,"GO_dotplot.png(pdf) has been successfully generated !\n", sep = "")

    p2 <- ggplot(plot_data,
                 aes(x = reorder(Description, GeneRatio),
                     y = GeneRatio,
                     fill = -log10(p.adjust))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      facet_grid(Cluster ~ ., scales = "free_y", space = "free_y") +
      scale_fill_gradient(low = "blue", high = "red", name = "-log10(p.adjust)") +
      labs(x = "", y = "Gene Ratio",
           title = "GO Enrichment Barplot") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold", size = 14),
            panel.grid.major.y = element_line(linetype = "dashed"),
            axis.text = element_text(size = 14, color = "black", face = "bold"),
            axis.title = element_text(size = 16, color = "black", face = "bold"),
            plot.title = element_text(size = 20,  color = "black", face = "bold"),
            legend.title = element_text(size = 12,  color = "black", face = "bold"),
            legend.text = element_text(size = 12,  color = "black", face = "bold"))+
      scale_x_discrete(
        labels = function(x) str_wrap(x, width = 50)  # 自动换行，每行最多50个字符
      )
    ggsave(file.path(outpath, paste0(prefix,"GO_barplot.png")), p2, width = 16, height = 12, dpi = 300, bg = "white")
    ggsave(file.path(outpath, paste0(prefix,"GO_barplot.pdf")), p2, width = 16, height = 12)
    cat("\nPlot 2 : ", prefix,"GO_barplot.png(pdf) has been successfully generated !\n", sep = "")
  }
}

# 4.KEGG分析
# (此处进行联网分析，也可以下载gmt文件进行本地分析)
if (!is.null(genefile)){
  gene_entrez <- bitr(genes,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = orgdb)
  if (!is.null(bg_genes)){
    if (!MSigDB){
      bg_entrez <- bitr(bg_genes,
                        fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = orgdb)
      kegg_result <- enrichKEGG(
        gene          = gene_entrez$ENTREZID,
        universe      = bg_entrez$ENTREZID,
        organism      = organism,
        keyType       = "kegg",
        pAdjustMethod = "BH",
        pvalueCutoff  = p_threshold,
        qvalueCutoff  = q_threshold,
      )
    }else{
      url <- switch(species,
                    human = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c2.all.v2025.1.Hs.symbols.gmt",
                    mouse = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/m2.all.v2025.1.Mm.symbols.gmt",
                    stop("Unsupported species"))
      destfile <- tempfile(fileext = ".gmt")
      download.file(url, destfile, method = "auto", mode = "wb")
      gmt <- read.gmt(destfile)
      kegg_result <- enricher(
        gene          = genes,
        universe      = bg_genes,
        TERM2GENE     = gmt,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_threshold,
        qvalueCutoff  = q_threshold
      )
    }
  }else{
    if (!MSigDB){
    kegg_result <- enrichKEGG(
      gene          = gene_entrez$ENTREZID,
      organism      = organism,
      keyType       = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff  = p_threshold,
      qvalueCutoff  = q_threshold
    )}else{
      url <- switch(species,
                    human = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c2.all.v2025.1.Hs.symbols.gmt",
                    mouse = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/m2.all.v2025.1.Mm.symbols.gmt",
                    stop("Unsupported species"))
      destfile <- tempfile(fileext = ".gmt")
      download.file(url, destfile, method = "auto", mode = "wb")
      gmt <- read.gmt(destfile)
      kegg_result <- enricher(
        gene          = genes,
        TERM2GENE     = gmt,
        pAdjustMethod = "BH",
        pvalueCutoff  = p_threshold,
        qvalueCutoff  = q_threshold
      )
    }
  }
  if(!MSigDB){
    kegg_result <- DOSE::setReadable(kegg_result, OrgDb = orgdb, keyType = "ENTREZID")}
  kegg_result_df <- as.data.frame(kegg_result@result)
  if (nrow(kegg_result_df) == 0) {
    warning("No significant pathways are enriched in KEGG analysis !")
  } else{
    write.csv(kegg_result_df, file.path(outpath, paste0(prefix, "KEGG_result.csv")), row.names = FALSE)
    cat("\nFile 2 : ", prefix, "KEGG_result.csv saved successfully !\n", sep = "")
    #可视化
      plot_data <- kegg_result_df %>%
      separate(GeneRatio, into = c("nGenes", "TotalGenes"), sep = "/") %>%
      mutate(
        GeneRatio = as.numeric(nGenes) / as.numeric(TotalGenes),
        Description_wrapped = stringr::str_wrap(gsub("_", " ", as.character(Description)), width = 50)
      ) %>%
      arrange(desc(GeneRatio)) %>%
      head(10)
  p3 <- ggplot(plot_data,
               aes(x = GeneRatio,
                   y = reorder(Description_wrapped, GeneRatio),  # 按GeneRatio排序
                   size = Count,                        # 点大小=基因数
                   color = -log10(p.adjust))) +         # 颜色=显著性
    geom_point(alpha = 0.8) +                    # 半透明效果
    scale_size_continuous(
      range = c(4, 10),                           # 点大小范围
      name = "Gene Count"
    ) +
    scale_color_gradient(
      low = "blue",
      high = "red",
      name = "-log10(p.adjust)"
    ) +
    labs(
      x = "Gene Ratio",
      y = "",
      title = "KEGG Pathway Enrichment Dotplot",
      caption = paste("Top", nrow(plot_data), "significant pathways")
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      panel.grid.major.y = element_line(linetype = "dashed"),
      axis.text = element_text(size = 14, color = "black", face = "bold"),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      plot.title = element_text(size = 20,  color = "black", face = "bold"),
      plot.caption = element_text(size = 10,  color = "black", face = "bold"),
      legend.title = element_text(size = 12,  color = "black", face = "bold"),
      legend.text = element_text(size = 12,  color = "black", face = "bold")
    )
  ggsave(file.path(outpath, paste0(prefix,"KEGG_dotplot.png")), p3, width = 16, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix,"KEGG_dotplot.pdf")), p3, width = 16, height = 12)
  cat("\nPlot 3 : ", prefix,"KEGG_dotplot.png(pdf) has been successfully generated !\n", sep = "")

  p4 <- ggplot(plot_data,
               aes(x = GeneRatio,
                   y = reorder(Description_wrapped, GeneRatio),
                   fill = -log10(p.adjust))) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_gradient(
      low = "blue",
      high = "red",
      name = "-log10(p.adjust)",
      limits = c(0, max(-log10(plot_data$p.adjust)))  # 固定颜色范围
    ) +
    labs(
      x = "",
      y = "Gene Ratio",
      title = "KEGG Pathways Enrichment (Barplot)",
      caption = paste("Top", nrow(plot_data), "most significant pathways")
    ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.text = element_text(size = 14, color = "black", face = "bold"),
      axis.title = element_text(size = 16, color = "black", face = "bold"),
      plot.title = element_text(size = 20,  color = "black", face = "bold"),
      plot.caption = element_text(size = 10,  color = "black", face = "bold"),
      legend.title = element_text(size = 12,  color = "black", face = "bold"),
      legend.text = element_text(size = 12,  color = "black", face = "bold"),
      panel.grid.major.y = element_blank()  # 隐藏横向网格线
    )
  ggsave(file.path(outpath, paste0(prefix,"KEGG_barplot.png")), p4, width = 16, height = 12, dpi = 300, bg = "white")
  ggsave(file.path(outpath, paste0(prefix,"KEGG_barplot.pdf")), p4, width = 16, height = 12)
  cat("\nPlot 4 : ", prefix,"KEGG_barplot.png(pdf) has been successfully generated !\n", sep = "")
}
}

if(!is.null(gsea_file) && !is.null(gsea_me)){
  gse_input <- function(df){
    geneList <- df[[gsea_me]]
    names(geneList) <- rownames(df)
    geneList <- sort(geneList, decreasing = TRUE)
    return(geneList)
  }
  df <- read.csv(gsea_file, row.names = 1)
  gene_list <- gse_input(df)
  
  # # 基于GO数据库的GSEA分析
  GSEA_GO <- gseGO(
    gene_list,
    ont = "ALL",
    # 可选"BP"、"MF"、"CC"三大类或"ALL"
    OrgDb = orgdb,
    # 使用人的OrgDb
    keyType = "SYMBOL",
    # 基因id类型
    pvalueCutoff = p_threshold,
    pAdjustMethod = "BH"
  )
  GSEA_GO_df <- as.data.frame(GSEA_GO)
  if (nrow(GSEA_GO_df) > 0) {
    GSEA_GO_df_filtered <- GSEA_GO_df %>%
      filter(abs(NES) > nes_threshold & pvalue < p_threshold & p.adjust < q_threshold)
    if (nrow(GSEA_GO_df_filtered) > 0){
      write.csv(GSEA_GO_df_filtered, file.path(outpath, paste0(prefix, "GSEA_GO_filtered_result.csv")), row.names = FALSE, quote = FALSE)
      cat("\nFile 1 : ", prefix, "GSEA_GO_filtered_result.csv saved successfully !\n", sep = "" )
      #可视化
      if(!require(RColorBrewer,quietly = TRUE)){install.packages("RColorBrewer")}
      set3_max <- brewer.pal.info["Set3", "maxcolors"]
      mycolors <- brewer.pal(set3_max, "Set3")
      GSEA_GO_plot_df <- GSEA_GO_df_filtered %>%
        arrange(desc(NES)) %>%
        {if(nrow(.) <= 10) . else bind_rows(head(., 5), tail(., 5))}

      p5 <- gseaNb(
        object = GSEA_GO,
        geneSetID = GSEA_GO_plot_df$ID,
        subPlot = 3,
        newGsea = F,
        curveCol = mycolors,
        rmHt = F,
        addPval = F,
        pvalX = 0.9,
        pvalY = 0.9,
        legend.position = "top",
        termWidth = 30
      )
      ggsave(file.path(outpath, paste0(prefix,"GSEA_GO.png")), p5, width = 16, height = 8, dpi = 300, bg = "white")
      ggsave(file.path(outpath, paste0(prefix,"GSEA_GO.pdf")), p5, width = 16, height = 8)
      cat("\nPlot 1 : ", prefix,"GSEA_GO.png(pdf) has been successfully generated !\n", sep = "")
    }else{
      warning('No significant GSEA result based on GO database was generated!')
      write.csv(GSEA_GO_df, file.path(outpath, paste0(prefix, "GSEA_GO_raw_result.csv")), row.names = FALSE, quote = FALSE)
      cat("\nFile 1 : ", prefix, "GSEA_GO_raw_result.csv saved successfully !\n", sep = "" )
      }
  }else{warning('No GSEA result based on GO database was generated!')}
  
  # 基于KEGG数据库的GSEA分析
  geneList <- gene_list
  if (MSigDB){
    url <- switch(species,
                  human = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c2.all.v2025.1.Hs.symbols.gmt",
                  mouse = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Mm/m2.all.v2025.1.Mm.symbols.gmt",
                  stop("Unsupported species"))
    destfile <- tempfile(fileext = ".gmt")
    download.file(url, destfile, method = "auto", mode = "wb")
    gmt <- read.gmt(destfile)
    GSEA_KEGG <- GSEA(
      geneList = geneList,  
      TERM2GENE = gmt,
      pvalueCutoff = p_threshold
    )
  }else{
    geneList_entrez <- bitr(
      names(geneList),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = orgdb
    )
    geneList <- geneList[geneList_entrez$SYMBOL]
    names(geneList) <- geneList_entrez$ENTREZID
    GSEA_KEGG <- gseKEGG(
      geneList = geneList,
      organism = organism,
      keyType = "ncbi-geneid",
      # 输入Entrez ID
      pvalueCutoff = p_threshold
    )
  }
  GSEA_KEGG_df <- as.data.frame(GSEA_KEGG)
  if (nrow(GSEA_KEGG_df) > 0) {
    GSEA_KEGG_df_filtered <- GSEA_KEGG_df %>% 
      filter(abs(NES) > nes_threshold & pvalue < p_threshold & p.adjust < q_threshold)
    if (nrow(GSEA_KEGG_df_filtered) > 0){
      write.csv(GSEA_KEGG_df_filtered, file.path(outpath, paste0(prefix, "GSEA_KEGG_filtered_result.csv")), row.names = FALSE, quote = FALSE)
      cat("\nFile 2 : ", prefix, "GSEA_KEGG_filtered_result.csv saved successfully !\n", sep = "" )
      #可视化
      if(!require(RColorBrewer,quietly = TRUE)){install.packages("RColorBrewer")}
      set3_max <- brewer.pal.info["Set3", "maxcolors"]
      mycolors <- brewer.pal(set3_max, "Set3")
      GSEA_KEGG_plot_df <- GSEA_KEGG_df_filtered %>%
        arrange(desc(NES)) %>%
        {if(nrow(.) <= 10) . else bind_rows(head(., 5), tail(., 5))}
      
      p6 <- gseaNb(
        object = GSEA_KEGG,
        geneSetID = GSEA_KEGG_plot_df$ID,
        subPlot = 3,
        newGsea = F,
        curveCol = mycolors,
        rmHt = F,
        addPval = F,
        pvalX = 0.9,
        pvalY = 0.9,
        legend.position = "top",
        termWidth = 30
      )
      ggsave(file.path(outpath, paste0(prefix,"GSEA_KEGG.png")), p6, width = 16, height = 8, dpi = 300, bg = "white")
      ggsave(file.path(outpath, paste0(prefix,"GSEA_KEGG.pdf")), p6, width = 16, height = 8)
      cat("\nPlot 2 : ", prefix,"GSEA_KEGG.png(pdf) has been successfully generated !\n", sep = "")
    }else{
      warning('No significant GSEA result based on KEGG database was generated!')
      write.csv(GSEA_KEGG_df, file.path(outpath, paste0(prefix, "GSEA_KEGG_raw_result.csv")), row.names = FALSE, quote = FALSE)
      cat("\nFile 1 : ", prefix, "GSEA_KEGG_raw_result.csv saved successfully !\n", sep = "" )
    }
  }else{warning('No GSEA result based on KEGG database was generated!')}
}



