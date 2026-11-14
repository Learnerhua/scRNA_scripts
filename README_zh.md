# scRNA-seq 分析脚本库

一套**持续更新**的通用单细胞转录组（scRNA-seq）分析脚本库。它包含一条按步骤编号的 **Seurat v5 主分析流程**——每样本 QC 与合并、去双细胞、归一化与批次整合、聚类与差异表达、细胞类型注释——以及若干**独立脚本**——数据准备与格式转换（h5ad 计数清洗、h5ad → Seurat rds、h5ad ⇄ rds）和下游分析（目标基因跨细胞类型表达可视化、样本 × 细胞类型统计、细胞组成成分差异分析等），并附带桑基图所需的配套 R 辅助脚本。所有脚本都以命令行参数接收输入（10x 目录、`.h5`、`rds`、`h5ad`、表格等路径），并把结果写成文件（数据对象、表格、图片等），因此每一步既可单独运行，也可串成完整流水线；**脚本会持续补充**。

[English](README.md) · 完整手册（中文）：[`USAGE.md`](USAGE.md)

## 脚本一览

| 脚本 | 类型 | 主要功能 | 输入 → 输出 | 示例图 |
|---|---|---|---|---|
| `scRNA_spRNA_h5ad_preprocess.py` | 独立（数据准备） | 确保表达矩阵是原始整数 counts（先用 `X`，否则回退 `raw.X`，两者都不是则中止），按指定 `var` 列重设基因名并合并同名基因；同时写出 `X` 与 `layers["counts"]` | `.h5ad`（目录 / 单文件）→ 含 `X` + `layers["counts"]` 的 `.h5ad` | — |
| `scRNA_h5ad2rds.r` | 独立（格式转换） | 批量 `.h5ad` → Seurat `.rds`；两条后端（默认纯 R 的 `anndataR`，或经 reticulate 的 `sceasy`）。输入须带 `layers["counts"]` | `.h5ad`（目录 / 单文件）→ `.rds` | — |
| `scRNA_h5ad_seurat_convert.r` | 独立（格式转换） | 单文件双向转换，按扩展名判方向：`.h5ad` → `.rds`、`.rds` → `.h5ad`（先转 factor 为 character；assay / 层 / 压缩可选） | `.h5ad` ⇄ `.rds` | — |
| `scRNA_mtx2rds.r` | 独立（数据准备） | 把 **Broad Single Cell Portal** 下载的数据组装成 Seurat 对象：稠密矩阵（`-E`）或 `.mtx` + 细胞/基因文件（`-M/-C/-G`），元数据必填（`-m`，细胞名列默认 `NAME`，首行 `TYPE` 自动丢弃）。基因名 `_`→`-` 并 `make.unique`；非整数输入只保留 `data` 层 | 矩阵 + 元数据 → `seurat.rds` | — |
| `scRNA_subcelltype_clustering.r` | 独立（细分亚型） | 对已注释对象里的每个细胞类型**单独**细分：逐类型重跑降维（归一化 → HVG → PCA → Harmony）+ 多分辨率 clustree，再按你记录的每类型分辨率定亚型、跑亚型 DE（marker / top10 / dotplot），最后把 `sub_celltype` 标签回写到原对象 | 注释 rds + 分辨率表 → 亚型 rds + 各亚型 DE | `Images/subcelltype_clustering/` |
| `scRNA_pseudobulk_DE.r` | 独立（差异统计） | **在每种细胞类型内部，找出两组之间差异表达的基因。** 按「细胞类型 × 样本」把原始 counts 求和，用 DESeq2 在**样本层面**比较 `-G1` vs `-G2`（每个细胞类型单独检验），结果可读作「在这类细胞里，这些基因有差异」；单元细胞数不足或每组样本数不足的类型自动跳过。输出 DE 长表 + 单元溯源表，并画火山图 | 已注释 rds（需 counts 层）→ DE csv + 火山图 | `Images/pseudobulk_DE/` |
| `scRNA_QC_step1.r` | 主流程 `step1` | 每样本独立 QC（阈值取自该样本自身的 nFeature_RNA / nCount_RNA 分位数，加 percent.mt 硬上限，可选 percent.HB 红细胞过滤），QC 通过后合并；基因级 `min.cells` 过滤在**合并之后**统一执行 | 10x 目录 / `.h5` 列表（CSV）→ 合并 Seurat 对象（仅 counts）+ QC 小提琴图/散点图 | `Images/01_QC/` |
| `scRNA_dedoublets_step2.r` | 主流程 `step2` | 每样本 **scDblFinder** 去双细胞；可选 **DecontX** 去环境 RNA（同样按样本独立跑） | 合并 rds → 仅含单细胞的 rds（`seu_clean.rds`） | `Images/02_deContamination/` |
| `scRNA_Nor_harmony_step3.r` | 主流程 `step3` | NormalizeData（LogNormalize）→ FindVariableFeatures（vst）→ ScaleData → RunPCA → **Harmony** 批次校正 | 干净 rds → 归一化 + 整合后的 rds（`seu_harmony.rds`） | `Images/03_Nor_Harmony/` |
| `scRNA_clustering_step4.r` | 主流程 `step4` | 两阶段：多分辨率 clustree 探索（不传 `-Res`）→ 聚类 + UMAP + `FindAllMarkers` 差异表达（传 `-Res`）；支持 Louvain / Leiden | 整合 rds → 聚类 rds + `clusters_DE/` marker 表与图 | `Images/04_clustering/`、`Images/04_clustering_louvain/` |
| `scRNA_annotate_step5.r` | 主流程 `step5` | 细胞类型注释，两种互斥模式：`-MB` 用 marker 列表出探索性 dotplot（不写 rds）；`-MA` 做 cluster → cell type 映射（写 `seu_annotated.rds`，并出 celltype UMAP / dotplot） | 聚类 rds + marker TSV → 注释 rds + 图 | `Images/05_annotation/` |
| `scRNA_gene_expr.r` | 独立脚本 | 目标基因在各细胞类型中的平均表达（CSV）+ 5 类图：Ridge / Violin / Feature / Dot / Heatmap | 已注释 rds → 表达量 CSV + 图 | `Images/gene_expr/` |
| `scRNA_gene_DE.r` | 独立脚本 | 单个目标基因：样本 × 细胞类型平均表达表、按 celltype z-score 的热图、每细胞类型一张柱状图、以及样本级箱线图（Wilcoxon 检验 + 显著性横线，**每个点是一个样本而非细胞**） | 已注释 rds → CSV + 热图 / 柱状图 / 箱线图 | `Images/gene_DE/` |
| `scRNA_composition_sccoda.py` | 独立脚本 | 用 **scCODA**（Dirichlet-Multinomial 贝叶斯模型）做细胞组成成分差异分析；`--deterministic` 可让 MCMC 逐字节可复现 | `.h5ad` / 样本 × 细胞类型计数表 / 每行一细胞的 obs 表 → 效应表 + 桑基图 / 箱线图 / 离散度 / 效应森林图 | `Images/composition_sccoda/` |
| `scRNA_spRNA_cellPro.r` | 独立（组成统计） | **每个样本里各种细胞类型各占多少？** 出样本级/组级的计数与占比表、100% 堆叠柱状图、组→细胞类型桑基图，并逐细胞类型比较两组占比（Wilcoxon，**每个点是一个样本**）；可再按亚型重复该比较 | 注释 rds → 占比表 + 柱状图 + 桑基图 + 占比比较图 | `Images/cellPro/` |
| `scRNA_SingleR.r` | 独立（自动注释） | **用参考数据集自动推断细胞类型名称**：把查询对象的每个聚类与参考数据比对（`-S human`/`-S mouse` 用 celldex 内置参考，或 `-RDS`/`-RD` 用自己的参考），输出 cluster→cell type 映射并写回对象，同时画得分热图与按推断标签着色的 UMAP | 已聚类 rds（含 `seurat_clusters`）+ 参考数据 → 映射 csv + 注释 rds + 图 | `Images/SingleR/` |
| `scRNA_spRNA_enrichment.r` | 独立（通路富集） | **把差异基因翻译成通路层面的结论**：逐细胞类型做 GO/KEGG 富集（以该类型自身检出的基因为背景）与 GSEA（全部基因按 log2FC 排序），输出各通路的 p/q 值、NES、dotplot 与 barplot。本脚本只准备输入文件，富集计算由外部命令 `enrichment_analysis` 执行 | 注释 rds + DE 结果 csv → 富集表 + 图 | `Images/enrichment/` |
| `scRNA_GSVA_enrichment.r` | 独立（通路活性） | **不依赖差异分析结果，比较两组之间的通路活性差异**：按「样本 × 细胞类型」聚合 counts，归一化为 log2(CPM+1)，用 GSVA 给每条通路打分（GO:BP、KEGG MEDICUS、REACTOME 三套基因集），再用 limma 比较两组，并画逐细胞类型的通路热图（对照组在左） | 注释 rds → 通路分数表 + 通路差异表 + 热图 | `Images/GSVA_enrichment/` |
| `scRNA_genesetScore.r` | 独立（基因集打分） | **给自己指定的基因集打分**：读取每行一个基因符号的文本文件，用四种方式打分——细胞类型层面的 GSVA 与 ssGSEA、单细胞层面的 AddModuleScore 与 AUCell——并画出细胞类型柱状图与单细胞 UMAP 特征图。适合回答「这个基因集在哪些细胞类型里更活跃」 | 注释 rds + 基因集文件 → 分数表 + 4 类图 | `Images/genesetScore/` |
| `enrichment_analysis.r` | 辅助脚本（**上一脚本必需**） | 用 clusterProfiler 执行实际的 GO/KEGG 富集与 GSEA：输入目标基因列表（`-G`）与可选背景基因集（`-E`），或输入全基因排序表做 GSEA（`-GSEA` + `-GM`），输出富集表、dotplot 与 barplot。也可单独用于 bulk RNA | 基因列表 / 排序表 → 富集表 + 图 | — |
| `scRNA_ROC.r` | 独立（判别分析） | **一组目标基因能把两组分开到什么程度？** 逐基因（加 `-GS` 时是整个基因集的合成分数）在样本层面算 AUC，给每个基因画一条 ROC 曲线，并输出按 AUC 排序的汇总表。这是**判别分析，不做显著性检验**——要 p 值请配合 `scRNA_gene_DE.r` 或 `scRNA_pseudobulk_DE.r` | 注释 rds + 分组 csv + 基因列表 → AUC 表 + ROC 曲线 | `Images/ROC/` |
| `scRNA_cellchat_single.r` | 独立（细胞通讯） | **哪些细胞类型之间在通讯？通过什么配受体？** 对单个注释对象跑 CellChat 流程：建立互作网络、列出显著的配受体对与信号通路、画总体与逐发送方的圈图、概率最高的 9 条通路的弦图，并把细胞群按「发出的信号」聚成通讯模式。需要跑两次——先不带 `-NP` 拿到 `selectK` 图，再按图上的模式数带上 `-NP` 重跑 | 注释 rds → 通讯表 + 网络 / 通路 / 模式图 | `Images/cellchat/single/` |
| `scRNA_spRNA_Monocle3.r` | 独立（拟时序） | **细胞之间的推进顺序是什么？** 用 Monocle3 学习轨迹图，从你指定的根细胞群沿拟时序给细胞排序，找出沿这条轴变化的基因（Moran's I 检验），并把它们归入基因模块。支持单细胞对象与空间对象（metadata 需含 `x`、`y` 坐标） | 注释 rds（或空间对象）→ 拟时序图 + Moran's I 表 + 基因模块 + 模块热图 | `Images/Monocle3/` |
| `scRNA_pyscenic_prepare.r` | 独立（调控网络） | **准备 pySCENIC 输入。** 从 Seurat 对象取出一层表达矩阵，剔除检出细胞过少的基因，可选抽样细胞，写成带 `Gene` / `CellID` 属性的 loom 文件（`pyscenic` CLI 要求的形式）| Seurat rds → loom | |
| `scRNA_pyscenic.py` | 独立（调控网络） | **推断转录因子调控网络。** 连续调用 `pyscenic` CLI 三次：`grn` 做共表达推断（GRNBoost2）、`ctx` 用 motif 富集把共表达裁剪成「直接调控」的靶基因（regulon）、`aucell` 给每个细胞的每个 regulon 打活性分。输出已存在的步骤会自动跳过，ctx 失败后可单独重跑 | loom + 参考文件 → adj.csv + regulons.csv + SCENIC.loom | |
| `scRNA_pyscenic_RSS.r` | 独立（调控网络） | **找出每类细胞的标志转录因子。** 从 AUCell 矩阵算 RSS 特异性分数，画特异性总图与逐 regulon 的各类图，并把调控活性矩阵写回 Seurat 对象——这样以后想画别的 TF 不必重跑 pySCENIC | Seurat rds + SCENIC loom → rss.csv + 图 + seu_SCENIC.rds | `Images/pyscenic/` |
| `scanpy/scRNA_pipeline_scanpy.py` | Scanpy 路线 | **一条命令走完质控、去双细胞、归一化、整合、聚类与 marker 发现。** 对应 Seurat 主流程 step1–4，但用 Python 写、处理 h5ad。QC 是两段式：第一次跑只画 QC 图后退出，供你定阈值；第二次带上阈值才继续。去双细胞用 Scrublet，去批次用 Harmony，聚类用 Leiden（`flavor='igraph'`）| h5ad → preprocessed.h5ad | `Images/scanpy/pipeline/` |
| `scanpy/scRNA_annotate_scanpy.py` | Scanpy 路线 | **给每个 cluster 指定细胞类型。** 阶段一吃「细胞类型 → marker」表并画 dotplot，供你判断 cluster 身份；阶段二吃「cluster → 细胞类型」表，写进 `adata.obs['celltype']` 后重画 UMAP 与 dotplot | preprocessed.h5ad → annotated.h5ad | `Images/scanpy/annotate/` |
| `scanpy/scRNA_DE_scanpy.py` | Scanpy 路线 | **找出每个 cluster 的 marker 基因。** 跑 `rank_genes_groups`（Wilcoxon / t-test / logreg），出 marker 概览图与逐 cluster 的 dotplot，并写出完整 marker 表。这是**细胞级**检验——需要样本级组间比较请用 `scRNA_pseudobulk_DE.r` | h5ad → marker 表 + 图 + DE_analyzed.h5ad | `Images/scanpy/DE/` |
| `scanpy/scRNA_pipeline_scanpy_GPU.py` | Scanpy 路线 | **同一套 13 步流程的 GPU 版。** 流程与输出和 CPU 版一致，只是把每个 `scanpy` 调用换成 `rapids-singlecell`；需要 CUDA 显卡与 `cupy`。额外 `-raw` 参数让保存的 h5ad 保留全部基因（默认只留高变基因）| h5ad → preprocessed.h5ad | |
| `scRNA_scTenifoldKnk.r` | 独立（虚拟敲除） | **不做真实敲除实验，预测拿掉某个基因后哪些基因会受影响。** 用 counts 矩阵重建调控网络，把目标基因的连接置零（虚拟敲除），比较野生型与敲除网络，输出受影响基因表及其富集通路 | Seurat rds（需 counts 层与已算好的 HVG）+ 一个基因名 → diffRegulation 表 + 结果 rds + 4 张图 | `Images/scTenifoldKnk/` |
| `scRNA_hdWGCNA.r` | 独立（共表达网络） | **找出某个细胞类型的基因共表达模块，以及每个模块的核心 hub 基因。** 在 metacell 上建 WGCNA 网络（metacell 把 k 个邻近细胞平均成一个伪细胞，避免 dropout 主导网络），把基因划分成模块、按 kME 排序、用 UCell 给每个细胞打模块活性分，并画出软阈值诊断图、基因聚类树与模块颜色条、模块相关性、细胞类型 × 模块 dotplot | Seurat rds（需 harmony 降维）→ 模块表 + hub 基因 + MEs/hMEs + 10 张图 + seu_hdWGCNA.rds | `Images/hdWGCNA/` |
| `scRNA_velocity.py` | 独立（轨迹） | **细胞正朝哪个方向变化？** 在含 spliced/unspliced 的 AnnData 上跑 scVelo RNA velocity，把速度流线叠在 UMAP 上——每个细胞不再只是一个点，而是带方向的箭头。`dynamical` 模式还会恢复每个基因的转录动力学、给出 latent time 和驱动转变的基因；`stochastic` 更快，只依赖表达量 | velocity 输入 h5ad → 流线图 + latent time 图 + 结果 h5ad | `Images/velocity/` |
| `scRNA_cellchat_compare.r` | 独立（细胞通讯） | **两组的通讯网络有什么不同？** 合并两个 CellChat 对象（由上一脚本产出），对比互作数量、强度、通路组成，以及组间显著差异的配受体对 | 两个 CellChat rds → 合并对象 + 对比图 | `Images/cellchat/compare/` |
| `scRNA_DoubleFinder.r` | 独立（替代方案） | **去双细胞的替代实现**：与 step2 做同一件事，但改用 DoubletFinder，且在整个合并对象上统一建模（不逐样本）；期望双细胞率必须手填 `-DR`（不能自动估计），输出保留全部层。用于用独立算法交叉验证 step2 的结论——绝不要在同一份数据上两个都跑 | 合并 rds → 干净 rds | `Images/DoubleFinder/` |
| `sankey_composition.r` | 辅助脚本（**上一脚本必需**） | 用 `ggalluvial` 渲染组成桑基（冲积）图。由 `scRNA_composition_sccoda.py` 自动调用；**必须与它放在同一目录** | 组成比例 CSV → 桑基图 PNG/PDF | ↑ |

## 流程总览

```text
可选的上游准备 —— 从 .h5ad 起步时（而非 10x 输出）

  原始 .h5ad ──▶ scRNA_spRNA_h5ad_preprocess.py ──▶ scRNA_h5ad2rds.r ──▶ Seurat rds
                    （整数 counts + layers["counts"]）                      │
                                                                          ▼ 作为 step1 的输入

主流程 —— 脚本名带 stepN 前缀，通过 rds 文件依次串联

  10x 目录 / h5 列表
        │  step1 scRNA_QC_step1.r           每样本 QC + 合并
        ▼
        │  step2 scRNA_dedoublets_step2.r   去双细胞（scDblFinder）[可选 DecontX]
        ▼
        │  step3 scRNA_Nor_harmony_step3.r  归一化 + HVG + PCA + Harmony
        ▼
        │  step4 scRNA_clustering_step4.r   clustree → 聚类 + UMAP + DE
        ▼
        │  step5 scRNA_annotate_step5.r     cluster → cell type
  seu_annotated.rds ──┬──▶ scRNA_gene_expr.r   目标基因跨细胞类型表达
                      └──▶ scRNA_gene_DE.r     样本 × 细胞类型表达

独立脚本 —— 不参与主流程编号
  scRNA_composition_sccoda.py   .h5ad / 计数表 → 组成成分差异分析（+ sankey_composition.r）
```

## 特点

- **每样本自适应 QC** —— 阈值来自各样本自身的分位数，测序深度差异大的样本间仍可比。
- **归一化只在 step3 做一次** —— 先删基因后归一化，避免分母不一致。
- **只重画图的缓存机制** —— 关键产物已存在时，多数脚本跳过耗时计算、仅重画图（改图尺寸几秒出图；删掉对应 rds 即强制重跑）。
- **统计以样本为单位而非细胞** —— `scRNA_gene_DE.r` 的箱线图每个点是一个样本，避免伪重复导致 p 值被细胞数撑大。检验默认双侧（`--alt` 可切单侧），结果表同时给 `p_value` 与 `p_adj`（对同一基因的各细胞类型做 BH 校正）。注意双侧秩和检验有硬下限：每组 3 个样本时最小 p 为 0.1，此时不可能达到显著。
- **scCODA 结果可复现** —— 加 `--deterministic` 后重复运行结果逐字节一致。
- **图与参考实现保持一致** —— 组成桑基图交给 R 的 `ggalluvial` 渲染，不在 Python 里重新实现。

## 环境依赖

**R 侧**（主流程各步与独立 R 脚本）

- 核心：`Seurat`（v5）、`tidyverse`（dplyr / ggplot2 / tidyr 等）、`patchwork`、`argparse`
- 去双细胞：`scDblFinder`、`BiocParallel`；可选 DecontX：`celda`
- 批次整合与探索：`harmony`、`clustree`
- 可选 Leiden 聚类：`reticulate`（脚本内部调用 `py_require("leidenalg")`）
- `scRNA_gene_DE.r` 额外需要 `reshape2`
- 桑基图（`scRNA_composition_sccoda.py` 配套 `sankey_composition.r`）：`ggalluvial`、`RColorBrewer`、`scales`、`ggplot2`、`dplyr`、`tidyr`
- 格式转换（`scRNA_h5ad2rds.r`、`scRNA_h5ad_seurat_convert.r`）：`anndataR`（默认后端）；可选 `sceasy`（经 reticulate，其 Python 环境需 `numpy` + `anndata`）
- 矩阵读入（`scRNA_mtx2rds.r`）：`Matrix`（`readMM`）、`data.table`（`fread`）
- 亚型聚类（`scRNA_subcelltype_clustering.r`）：`clustree`、`harmony`、`reticulate`（Leiden 需 `leidenalg`）、`patchwork`
- pseudobulk 差异表达（`scRNA_pseudobulk_DE.r`）：`DESeq2`（Bioconductor）、`Seurat`、`Matrix`、`tidyverse`、`ggrepel`、`RColorBrewer`、`scales`；细胞类型 > 66 时另需 `viridis`
- 细胞比例分析（`scRNA_spRNA_cellPro.r`）：`ggalluvial`（桑基图）、`ggpubr`（显著性标注）、`patchwork`（亚型图拼接）、`tidyverse`、`Seurat`、`RColorBrewer`、`scales`
- 判别分析（`scRNA_ROC.r`）：`pROC`（AUC / ROC）、`tidyverse`、`Seurat`
- 细胞通讯（`scRNA_cellchat_single.r` / `scRNA_cellchat_compare.r`）：`CellChat`、`NMF`（通讯模式聚类）、`ggalluvial`、`patchwork`、`ComplexHeatmap`（outgoing 通路热图）、`Seurat`、`tidyverse`
- 拟时序与轨迹（`scRNA_spRNA_Monocle3.r`）：`monocle3`、`SeuratWrappers`（`as.cell_data_set`）、`pheatmap`（模块热图）、`igraph`、`leidenbase`、`Rcpp`、`Seurat`、`tidyverse`
- pySCENIC（`scRNA_pyscenic_prepare.r` + `scRNA_pyscenic_RSS.r`）：`SCopeLoomR`（读写 loom）、`SCENIC`（`calcRSS` / `plotRSS` / `getAUC`）、`gridExtra`、`Seurat`、`tidyverse`；Python 侧：`pyscenic` 0.12.x（含 `arboreto` / `ctxcore` / `primp`）、`loompy`、`dask`、`pandas`、`numpy`（需独立的 Python 3.7 环境）
- Scanpy 路线（`scanpy/` 子目录）：`scanpy`、`anndata`、`pandas`、`numpy`、`matplotlib`、`harmonypy`（Harmony）、`igraph`（Leiden 用 `flavor='igraph'`，不需要 `leidenalg`）、`scrublet`（去双细胞步骤）；GPU 版另需 `cupy` + `rapids-singlecell` + CUDA 显卡
- 虚拟敲除（`scRNA_scTenifoldKnk.r`）：`scTenifoldKnk`（GitHub `cailab-tamu/scTenifoldKnk`）、`Seurat`、`tidyverse`、`patchwork`、`igraph`、`ggrepel`、`enrichR`（`plotKO` 会联网访问 Enrichr）
- 共表达网络（`scRNA_hdWGCNA.r`）：`hdWGCNA`（GitHub `smorabit/hdWGCNA`）、`WGCNA`、`UCell`、`ggradar`、`Seurat`、`harmonypy`（metacell 步骤依赖 harmony）、`tidyverse`、`patchwork`
- RNA velocity（`scRNA_velocityPrepare.r` + `scRNA_velocity.py`）：`anndataR`（AnnData 转换，R 侧）、`scanpy`、`scvelo`、`anndata`、`numpy`、`pandas`、`h5py`、`matplotlib`（Python 侧）
- 自动注释（`scRNA_SingleR.r`）：`SingleR`、`SingleCellExperiment`、`scater`、`celldex`（内置 human/mouse 参考）、`Seurat`、`tidyverse`
- 通路富集（`scRNA_spRNA_enrichment.r` + 配套 `enrichment_analysis.r`）：主脚本只准备输入文件，富集由配套脚本执行，需 `clusterProfiler`、`enrichplot`、`GseaVis` 与 `org.Hs.eg.db`/`org.Mm.eg.db`
- 通路活性打分（`scRNA_GSVA_enrichment.r`）：`GSVA`、`msigdbr`（GO:BP / KEGG MEDICUS / REACTOME 基因集）、`limma`（组间比较）、`ComplexHeatmap`（热图）、`Seurat`、`tidyverse`
- 自定义基因集打分（`scRNA_genesetScore.r`）：`GSVA`（GSVA + ssGSEA）、`AUCell`、`Seurat`（AddModuleScore）、`tidyverse`

**Python 侧**（`scRNA_spRNA_h5ad_preprocess.py`、`scRNA_composition_sccoda.py`）

- `scCODA`、`tensorflow` / `tensorflow-probability`
- `anndata`（注意与 zarr / numpy 的版本匹配：`anndata < 0.12` 拒 zarr ≥3，`anndata ≥ 0.13` 要求 numpy ≥2.1）
- h5ad 预处理：`scanpy`、`scipy`、`numpy`、`pandas`
- 组成分析：`pandas`、`matplotlib`、`seaborn`
- 桑基图额外需要带 `ggalluvial` 的 `Rscript`；缺失时仅跳过该图，其余产物不受影响

**解释器** —— 脚本 shebang 为 `#!/usr/bin/env Rscript` 与 `#!/usr/bin/env python3`，因此装了上述依赖的环境需在 `PATH` 上（如 `export PATH=/path/to/env/bin:$PATH`），或显式调用解释器：`Rscript script.r …`、`/path/to/env/bin/python script.py …`。

## 快速开始

```bash
# 可选：从 .h5ad 起步
./scRNA_spRNA_h5ad_preprocess.py -I h5ad_dir -V gene_symbols -P 00. -O 00_preprocess
./scRNA_h5ad2rds.r               -I 00_preprocess -P 00. -O 00_rds

# 可选：Broad Single Cell Portal 导出数据（稠密矩阵 + 元数据）
./scRNA_mtx2rds.r                -E exprMatrix.tsv.gz -m metadata.tsv -CC NAME -P 00. -O 00_rds

# 可选：对已注释对象逐细胞类型做亚型聚类（阶段 1 → clustree）
./scRNA_subcelltype_clustering.r -I 05_annotation/05.seu_annotated.rds -AC celltype -P 06. -O 06_subtype

# 可选：每种细胞类型内部，两组之间哪些基因有差异（pseudobulk，DESeq2）
./scRNA_pseudobulk_DE.r -I 05_annotation/05.seu_annotated.rds -CC celltype -SC orig.ident -GC group -G1 treat -G2 control -P 07. -O 07_pseudobulk

# 可选：每个样本的细胞类型构成 + 两组占比比较
./scRNA_spRNA_cellPro.r -I 05_annotation/05.seu_annotated.rds -SC orig.ident -CC celltype -GC group -GO control,treat -EC Unknown -P 08. -O 08_cellPro

# 可选：目标基因对两组的判别能力（AUC / ROC 曲线）
./scRNA_ROC.r -I 05_annotation/05.seu_annotated.rds -SC orig.ident -G group.csv -F genelist.txt -C control,treat -P 09. -O 09_ROC

# 可选：细胞通讯（跑两次：先不带 -NP，再带上）
./scRNA_cellchat_single.r -I 05_annotation/05.seu_annotated.rds -AC celltype -S mouse -TP triMean -TM 0.1 -NP 3 -T 20 -P 10. -O 10_cellChat

# 可选：拟时序 / 轨迹（用 -C 和 -G 指定根细胞群）
./scRNA_spRNA_Monocle3.r -I 05_annotation/05.seu_annotated.rds -DT SC -C celltype -G <root_celltype> -R1 1e-5 -R2 0.5 -T 20 -P 11. -O 11_monocle3

# 可选：RNA velocity（先准备输入，再跑 scVelo）
./scRNA_velocityPrepare.r -R 05_annotation/05.seu_annotated.rds -F full_seurat_with_velocity.rds -P 12. -O 12_velocity
./scRNA_velocity.py -I 12_velocity/12.scRNA_velocity_input.h5ad -AC celltype -M dynamical -T 40 -P 12. -O 12_velocity/result

# 可选：转录因子调控网络（三个脚本，按顺序执行）
./scRNA_pyscenic_prepare.r -I 05_annotation/05.seu_annotated.rds -A RNA -L data -MG 50 -P 01. -O 12_pyscenic
./scRNA_pyscenic.py -I 12_pyscenic/01.expr_norm.loom -S mouse -RD ref_data -T 12 -P 01. -O 12_pyscenic
./scRNA_pyscenic_RSS.r -I 05_annotation/05.seu_annotated.rds -L 12_pyscenic/01.SCENIC.loom -AC celltype -NT 2 -P 02. -O 12_pyscenic/RSS

# 可选：虚拟基因敲除
./scRNA_scTenifoldKnk.r -I 05_annotation/05.seu_annotated.rds -TG GeneX -T 8 -p_val 0.05 -P 01. -O 14_scTenifoldKnk

# 可选：单个细胞类型的共表达模块
./scRNA_hdWGCNA.r -I 05_annotation/05.seu_annotated.rds -AC celltype -SC orig.ident -TC 'CellTypeA' -A RNA -SP 3 -T 8 -P 01. -O 13_hdWGCNA

# 可选：Scanpy 路线（Python，处理 h5ad 而非 rds）
./scanpy/scRNA_pipeline_scanpy.py -I raw.h5ad -BK sample -O 01_scanpy -P 01.          # 第一次：只画 QC 图
./scanpy/scRNA_pipeline_scanpy.py -I raw.h5ad -BK sample -O 01_scanpy -P 01. -minG 200 -maxG 2500 -minU 500 -maxU 10000 -minC 3 -maxMT 5 -maxHB 5
./scanpy/scRNA_annotate_scanpy.py -I 01_scanpy/01.preprocessed.h5ad -MB markers.tsv -P 02. -O 02_annotation
./scanpy/scRNA_DE_scanpy.py -I 02_annotation/02.annotated.h5ad -CK celltype -P 03. -O 03_DE

# 可选：用参考数据集自动推断细胞类型（SingleR）
./scRNA_SingleR.r -I 04_clustering/04.seu_res0.6_clustered.rds -S mouse -P 05. -O 05_SingleR

# 可选：对 DE 结果做 GO/KEGG 富集 + GSEA（逐细胞类型）
./scRNA_spRNA_enrichment.r -I 05_annotation/05.seu_annotated.rds -DER 06_pseudobulk/06.KO_vs_WT_DE_raw.csv -CC celltype -GC group -G control,treat -A RNA -S mouse -IC Unknown -P 07. -O 07_enrichment

# 可选：比较两组之间的通路活性差异（GSVA，不需要 DE 结果）
./scRNA_GSVA_enrichment.r -I 05_annotation/05.seu_annotated.rds -AC celltype -SC orig.ident -GC group -G control,treat -S mouse -P 08. -O 08_GSVA

# 可选：给自己的基因集打分（GSVA / ssGSEA / AddModuleScore / AUCell）
./scRNA_genesetScore.r -I 05_annotation/05.seu_annotated.rds -G geneSet.txt -AC celltype -A RNA -P 09. -O 09_genesetScore

# 主流程
./scRNA_QC_step1.r          -I input_list.csv -DT 10X -P 01. -O 01_QC
./scRNA_dedoublets_step2.r  -I 01_QC/01.seu_afterQC.rds  -P 02. -O 02_deContamination
./scRNA_Nor_harmony_step3.r -I 02_deContamination/02.seu_clean.rds -P 03. -O 03_Nor_Harmony

# step4 分两阶段：先探索，再定分辨率
./scRNA_clustering_step4.r  -I 03_Nor_Harmony/03.seu_harmony.rds -P 04. -O 04_clustering
./scRNA_clustering_step4.r  -I 03_Nor_Harmony/03.seu_harmony.rds -Res 0.7 -P 04. -O 04_clustering

# step5：先用 marker 列表探索，再写最终注释
./scRNA_annotate_step5.r    -I 04_clustering/04.seu_res0.7_clustered.rds -MB markers.tsv    -P 05. -O 05_annotation
./scRNA_annotate_step5.r    -I 04_clustering/04.seu_res0.7_clustered.rds -MA annotation.tsv -P 05. -O 05_annotation

# 独立脚本
./scRNA_gene_expr.r         -I 05_annotation/05.seu_annotated.rds -G GENE1,GENE2,GENE3 -P 06. -O 06_gene_expr
./scRNA_gene_DE.r           -I 05_annotation/05.seu_annotated.rds -G GENE1 -CT celltypeA,celltypeB -P 07. -O 07_gene_DE
./scRNA_composition_sccoda.py -I input.h5ad -S sample -CT celltype -GC group --deterministic -P 08. -O 08_composition

# 单文件格式互转
./scRNA_h5ad_seurat_convert.r -I 05_annotation/05.seu_annotated.rds -O 09_export -P export.
```

每个脚本都能自解释：执行 `<脚本> -h` 查看完整参数。

### step1 的输入列表

```csv
path,sample,group
/path/to/sample1/filtered_feature_bc_matrix,sample1,control
/path/to/sample2/filtered_feature_bc_matrix,sample2,treat
```

`-DT 10X` 表示每个样本给一个 10x 目录；`-DT h5` 表示每个样本一个 `.h5` 文件。`sample` 必须唯一，`group` 会写入对象的 metadata。

## 使用约定

- 文件名带 `stepN` 前缀的是主流程脚本，按顺序串联使用；其余为独立脚本，可作用于任何符合条件的输入。
- 各步之间以文件衔接，中间产物可单独复用；每一步的输出前缀（如 `01.`、`02.`）与其步骤号对应。
- 已存在的关键产物会被当作缓存（仅重画图）；删掉对应 rds 即可强制重跑。
- Harmony 的批次变量必须用样本/供体列，**不要**用生物学分组列（step3 会对 `-BC group` 直接报错拦截，避免把分组差异校正掉）。
- step5 不重算 DE，直接复用 step4 产出的 `clusters_DE/` 表。
- 从 **Broad Single Cell Portal** 下载数据时，用 `scRNA_mtx2rds.r` 直接建对象（细胞名列默认 `NAME`，首行 `TYPE` 自动跳过）。
- `scRNA_subcelltype_clustering.r` 用于把已注释对象再细分：阶段 1 为每个细胞类型出 clustree，阶段 2 按你记录的每类型分辨率（CSV：`celltype`、`resolution`）定亚型并跑亚型 DE。
- `scRNA_pseudobulk_DE.r` 回答「**每种细胞类型内部，两组之间哪些基因有差异**」：按「细胞类型 × 样本」聚合原始 counts 后跑 DESeq2，所以 `-SC` 必须填**生物学重复**。细胞级检验把每个细胞当重复，p 值会虚高；pseudobulk 的结论才能外推到新样本。
- `scRNA_scTenifoldKnk.r` **一次只敲一个基因**（`-TG` 只取一个名字，多个基因换 `-P` 前缀跑多次）。需要原始 counts、对象已算过 `VariableFeatures`，且细胞数要够（每个子网络抽 500 个）。
- `scRNA_hdWGCNA.r` **一次只分析一个细胞类型**：`-TC` 只取一个值，且输入对象必须已带 `harmony` 降维，并含 `-AC` / `-SC` 指定的注释列与样本列。
- Scanpy 路线放在 `scanpy/` 里，是与 Seurat 主线并行的另一条流程：收发 `.h5ad` 而非 `.rds`，并把 Seurat 的 step1–4 合进一个脚本。**两种对象格式不要混用**——需要跨路线时用 `scRNA_h5ad2rds.r` / `scRNA_h5ad_seurat_convert.r` 转换。
- pySCENIC 是三脚本一条链，**必须按顺序执行**：`scRNA_pyscenic_prepare.r`（Seurat → loom）→ `scRNA_pyscenic.py`（grn / ctx / aucell）→ `scRNA_pyscenic_RSS.r`（RSS、出图、活性写回对象）。三者跨两个环境，靠各自 shebang 切换，无需手动 activate。`SCENIC.loom` **只能用 R 的 `SCopeLoomR` 读**，Python 的 `loompy` 读不了（多列属性不符合 loom v3 规范）。
- RNA velocity 分两步：`scRNA_velocityPrepare.r` 把「含 spliced/unspliced 的完整对象」与「分析用的子集对象」对齐，写出带三个 assay 的 h5ad；`scRNA_velocity.py` 再做实际计算。velocity 推断的是表达状态的变化方向，**本身不等于谱系关系**。
- 拟时序分析是有方向的：根细胞群由 `-C` + `-G` 指定，选错根整条时间轴方向就错了——先用轨迹图看懂结构再定根。Monocle3 **不建模细胞群的占比变化**，那要用 scCODA。
- 细胞通讯分两步：`scRNA_cellchat_single.r` 对单个对象建模（**第一次先不带 `-NP`** 拿到 `selectK` 图，读出模式数后再带上 `-NP` 重跑），产出 `cellchat.rds`；`scRNA_cellchat_compare.r` 拿两个这样的对象做对比。两个对象的细胞身份必须完全一致。
- `scRNA_ROC.r` 是判别分析：统计单位是样本，AUC 说的是「这个基因能把两组分开到什么程度」，**不做显著性检验**。要 p 值请用 `scRNA_gene_DE.r`（单基因、样本级 Wilcoxon）或 `scRNA_pseudobulk_DE.r`（逐细胞类型、DESeq2）。
- `scRNA_spRNA_cellPro.r` 回答「**每个样本里各细胞类型占多少、两组之间比例有没有差别**」：统计单位是样本（`-SC`），不是细胞。它是快速描述 + 简单 Wilcoxon 检验；需要带可信区间的正式成分检验时用 `scRNA_composition_sccoda.py`。
- `scRNA_SingleR.r` 用参考数据集自动给出细胞类型名称，输入需已含 `seurat_clusters`（即 step4 的产物）。它的标签属于推断结果，需与 marker 交叉核对，是 step5 marker 注释的补充而非替代。
- 去双细胞有两个版本——`scRNA_dedoublets_step2.r`（scDblFinder，逐样本，双细胞率自动估计，输出只保留 counts）与 `scRNA_DoubleFinder.r`（DoubletFinder，整体处理，双细胞率手填）。二选一，不要在同一份数据上都跑；后者用于用独立算法交叉验证。
- `scRNA_spRNA_enrichment.r` 把 DE 表翻译成通路层面的结论：显著基因做 GO/KEGG 富集，全部基因按 `-GM` 排序做 GSEA。DE 输入需含 `gene`、`cluster`、`p_val_adj`、`avg_log2FC`，所以 `scRNA_pseudobulk_DE.r` 的产物可以直接喂进来。富集在独立后台进程中执行，脚本跑完不代表结果已就绪。
- `scRNA_GSVA_enrichment.r` 回答同样的问题但不需要 DE 表：它给每个细胞类型算通路活性再比较两组。`-G` 的顺序有含义——第一个是对照组（热图里在左），第二个是处理组。热图中的通路不一定显著：脚本会从 `adj.P.Val < 0.05` 依次退到 `P.Value < 0.05`、再到单纯按 `|logFC|` 排序，引用某条通路前请回查差异表里的 adj.P.Val。
- `scRNA_genesetScore.r` 打的是你自己提供的基因集，不是通路库。同一个基因集会给四种分数（GSVA、ssGSEA、AddModuleScore、AUCell）、两个层级（细胞类型、单细胞）。脚本会打印基因集里有多少基因在对象中找不到——先看这个数，基因集与对象不匹配时分数没有意义。
- 从 `.h5ad` 起步时先跑 `scRNA_spRNA_h5ad_preprocess.py`：`scRNA_h5ad2rds.r` 需要 `layers["counts"]` 层；若转换出的 Seurat 对象只有名为 `X` 的层，主流程 step1（按 `layer = "counts"` 取计数）无法使用。

## 示例可视化输出

`Images/` 目录收录了各绘图脚本的示例输出，子目录与产出图的脚本一一对应（数据准备 / 格式转换脚本不产图），详见上表“示例图”一列。示例图均由脚本自身产出，用于展示输出形式与排版。

**示例图只保留 PNG**：脚本运行时同一张图会同时输出 `.png` 与 `.pdf` 两个文件，但示例集中**只保留 PNG**，同名 PDF 已全部删除（既省体积、也便于网页直接浏览）。因此各章「输出」小节里写成 `.png/.pdf` 的文件，在 `Images/` 中只能找到 `.png`；需要 PDF 时用相应脚本重新跑一次即可。

## 文档

- [`USAGE.md`](USAGE.md) —— 完整中文手册：逐脚本参数说明、输出清单、缓存行为与已知坑。
- 各脚本的 `-h` 输出 —— 与代码同步，最权威。

## 构建方式与致谢

本脚本库在 **Hermes 智能体**（[Nous Research](https://hermes-agent.nousresearch.com/docs)）以及 **DeepSeek**、**MiniMax** 等大语言模型的辅助下构建；所有脚本在纳入本库前均由作者实际运行并核验。

## 后续规划

本库持续更新，**后续还会补充更多脚本**（更多分析环节与下游工具）。欢迎提出建议或需求。

---

Copyright © 2026 OYJH. All Rights Reserved.
