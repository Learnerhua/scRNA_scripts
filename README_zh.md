# scRNA-seq 分析脚本库

一套**持续更新**的通用单细胞转录组（scRNA-seq）分析脚本库。它包含一条按步骤编号的 **Seurat v5 主分析流程**——每样本 QC 与合并、去双细胞、归一化与批次整合、聚类与差异表达、细胞类型注释——以及若干**独立下游脚本**（目标基因跨细胞类型表达可视化、样本 × 细胞类型统计、细胞组成成分差异分析等），并附带桑基图所需的配套 R 辅助脚本。所有脚本都以命令行参数接收输入（10x 目录、`.h5`、`rds`、`h5ad`、表格等路径），并把结果写成文件（数据对象、表格、图片等），因此每一步既可单独运行，也可串成完整流水线；**脚本会持续补充**。

[English](README.md) · 完整手册（中文）：[`USAGE.md`](USAGE.md)

## 脚本一览

| 脚本 | 类型 | 主要功能 | 输入 → 输出 | 示例图 |
|---|---|---|---|---|
| `scRNA_QC_step1.r` | 主流程 `step1` | 每样本独立 QC（阈值取自该样本自身的 nFeature_RNA / nCount_RNA 分位数，加 percent.mt 硬上限，可选 percent.HB 红细胞过滤），QC 通过后合并；基因级 `min.cells` 过滤在**合并之后**统一执行 | 10x 目录 / `.h5` 列表（CSV）→ 合并 Seurat 对象（仅 counts）+ QC 小提琴图/散点图 | `Images/01_QC/` |
| `scRNA_dedoublets_step2.r` | 主流程 `step2` | 每样本 **scDblFinder** 去双细胞；可选 **DecontX** 去环境 RNA（同样按样本独立跑） | 合并 rds → 仅含单细胞的 rds（`seu_clean.rds`） | `Images/02_deContamination/` |
| `scRNA_Nor_harmony_step3.r` | 主流程 `step3` | NormalizeData（LogNormalize）→ FindVariableFeatures（vst）→ ScaleData → RunPCA → **Harmony** 批次校正 | 干净 rds → 归一化 + 整合后的 rds（`seu_harmony.rds`） | `Images/03_Nor_Harmony/` |
| `scRNA_clustering_step4.r` | 主流程 `step4` | 两阶段：多分辨率 clustree 探索（不传 `-Res`）→ 聚类 + UMAP + `FindAllMarkers` 差异表达（传 `-Res`）；支持 Louvain / Leiden | 整合 rds → 聚类 rds + `clusters_DE/` marker 表与图 | `Images/04_clustering/`、`Images/04_clustering_louvain/` |
| `scRNA_annotate_step5.r` | 主流程 `step5` | 细胞类型注释，两种互斥模式：`-MB` 用 marker 列表出探索性 dotplot（不写 rds）；`-MA` 做 cluster → cell type 映射（写 `seu_annotated.rds`，并出 celltype UMAP / dotplot） | 聚类 rds + marker TSV → 注释 rds + 图 | `Images/05_annotation/` |
| `scRNA_gene_expr.r` | 独立脚本 | 目标基因在各细胞类型中的平均表达（CSV）+ 5 类图：Ridge / Violin / Feature / Dot / Heatmap | 已注释 rds → 表达量 CSV + 图 | `Images/gene_expr/` |
| `scRNA_gene_DE.r` | 独立脚本 | 单个目标基因：样本 × 细胞类型平均表达表、按 celltype z-score 的热图、每细胞类型一张柱状图、以及样本级箱线图（Wilcoxon 检验 + 显著性横线，**每个点是一个样本而非细胞**） | 已注释 rds → CSV + 热图 / 柱状图 / 箱线图 | `Images/gene_DE/` |
| `scRNA_composition_sccoda.py` | 独立脚本 | 用 **scCODA**（Dirichlet-Multinomial 贝叶斯模型）做细胞组成成分差异分析；`--deterministic` 可让 MCMC 逐字节可复现 | `.h5ad` / 样本 × 细胞类型计数表 / 每行一细胞的 obs 表 → 效应表 + 桑基图 / 箱线图 / 离散度 / 效应森林图 | `Images/composition_sccoda/` |
| `sankey_composition.r` | 辅助脚本（**上一脚本必需**） | 用 `ggalluvial` 渲染组成桑基（冲积）图。由 `scRNA_composition_sccoda.py` 自动调用；**必须与它放在同一目录** | 组成比例 CSV → 桑基图 PNG/PDF | ↑ |

## 流程总览

```text
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
- **统计以样本为单位而非细胞** —— `scRNA_gene_DE.r` 的箱线图每个点是一个样本，避免伪重复导致 p 值被细胞数撑大。
- **scCODA 结果可复现** —— 加 `--deterministic` 后重复运行结果逐字节一致。
- **图与参考实现保持一致** —— 组成桑基图交给 R 的 `ggalluvial` 渲染，不在 Python 里重新实现。

## 环境依赖

**R 侧**（主流程各步与独立 R 脚本）

- 核心：`Seurat`（v5）、`tidyverse`（dplyr / ggplot2 / tidyr 等）、`patchwork`、`argparse`
- 去双细胞：`scDblFinder`、`BiocParallel`；可选 DecontX：`celda`
- 批次整合与探索：`harmony`、`clustree`
- 可选 Leiden 聚类：`reticulate`（脚本内部调用 `py_require("leidenalg")`）
- `scRNA_gene_DE.r` 额外需要 `reshape2`
- 桑基图（`sankey_composition.r`）：`ggalluvial`、`RColorBrewer`、`scales`、`ggplot2`、`dplyr`、`tidyr`

**Python 侧**（`scRNA_composition_sccoda.py`）

- `scCODA`、`tensorflow` / `tensorflow-probability`
- `anndata`（注意与 zarr / numpy 的版本匹配：`anndata < 0.12` 拒 zarr ≥3，`anndata ≥ 0.13` 要求 numpy ≥2.1）、`numpy`、`pandas`、`matplotlib`、`seaborn`
- 桑基图额外需要带 `ggalluvial` 的 `Rscript`；缺失时仅跳过该图，其余产物不受影响

**解释器** —— 脚本 shebang 为 `#!/usr/bin/env Rscript` 与 `#!/usr/bin/env python3`，因此装了上述依赖的环境需在 `PATH` 上（如 `export PATH=/path/to/env/bin:$PATH`），或显式调用解释器：`Rscript script.r …`、`/path/to/env/bin/python script.py …`。

## 快速开始

```bash
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

## 示例可视化输出

`Images/` 目录收录了每个脚本的示例输出（每张图 PNG + PDF 同名成对），子目录与脚本一一对应，详见上表“示例图”一列。示例图均由脚本自身产出，用于展示输出形式与排版。

## 文档

- [`USAGE.md`](USAGE.md) —— 完整中文手册：逐脚本参数说明、输出清单、缓存行为与已知坑。
- 各脚本的 `-h` 输出 —— 与代码同步，最权威。

## 构建方式与致谢

本脚本库在 **Hermes 智能体**（[Nous Research](https://hermes-agent.nousresearch.com/docs)）以及 **DeepSeek**、**MiniMax** 等大语言模型的辅助下构建；所有脚本在纳入本库前均由作者实际运行并核验。

## 后续规划

本库持续更新，**后续还会补充更多脚本**（更多分析环节与下游工具）。欢迎提出建议或需求。

---

Copyright © 2026 OYJH. All Rights Reserved.
