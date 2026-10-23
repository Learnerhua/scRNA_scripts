# scRNA 通用分析脚本库 —— 使用说明

本目录是一套**持续更新**的通用单细胞转录组（scRNA-seq）分析脚本库。其中**主分析流程脚本按步骤编号**（脚本名带 `stepN` 前缀），其余为**独立脚本**，不在流程编号内、输入各自独立，可单独使用。文档每个脚本一章，覆盖从原始计数矩阵到细胞注释、表达可视化与细胞组成差异分析的完整链路；所有脚本都以命令行参数接收输入（10x 目录、`.h5`、`.rds`、`.h5ad`、表格等路径），并把结果写成文件（数据对象、表格、图片等），因此每一步既可单独运行，也可串成完整流水线；脚本会持续补充。各脚本的示例可视化输出集中存放在本目录的 `Images/` 中（见下文“可视化输出示例”）。

脚本自带 shebang（R 脚本 `#!/usr/bin/env Rscript`，Python 脚本 `#!/usr/bin/env python3`）与自解释参数（`-h` 查看帮助）。只要 PATH 上有对应解释器，即可 `./<script>` 直接执行；也可用 `Rscript <script>.r …` / `python3 <script>.py …` 显式指定解释器。环境与依赖清单见文末附录。

## 流程概览

**主流程**（脚本名带 `stepN`，按顺序串联；前一步的 rds 即后一步的输入）：

| 步骤 | 脚本 | 用途 | 输入 → 输出 |
|---|---|---|---|
| step1 | `scRNA_QC_step1.r` | 每样本独立 QC + 合并 | 10x 目录 / h5 列表 → 合并 rds |
| step2 | `scRNA_dedoublets_step2.r` | 去双细胞（可选 DecontX 去环境 RNA） | 合并 rds → 干净 rds |
| step3 | `scRNA_Nor_harmony_step3.r` | 归一化 + Harmony 批次整合 | 干净 rds → 整合 rds |
| step4 | `scRNA_clustering_step4.r` | 多分辨率探索 + 聚类 + 差异表达 | 整合 rds → 聚类 rds + markers |
| step5 | `scRNA_annotate_step5.r` | 细胞类型注释（探索 / 标注两种模式） | 聚类 rds + marker 表 → 注释 rds |

**独立脚本**（不参与主流程编号，可单独使用）：

| 脚本 | 用途 | 输入 → 输出 |
|---|---|---|
| `scRNA_spRNA_h5ad_preprocess.py` | h5ad 计数矩阵预处理：判定并取出原始整数 counts（`X`，否则回退 `raw.X`）、按 `-V` 重设基因名并合并同名基因 | h5ad（目录 / 单文件）→ 含 `X` + `layers["counts"]` 的 h5ad |
| `scRNA_h5ad2rds.r` | 批量 h5ad → Seurat rds（`anndataR` / `sceasy` 两条路径） | h5ad（目录 / 单文件）→ rds |
| `scRNA_h5ad_seurat_convert.r` | 单文件双向格式转换（按扩展名自动判方向） | h5ad ⇄ rds |
| `scRNA_mtx2rds.r` | 把 Broad Single Cell Portal 下载的矩阵 + 元数据组装成 Seurat 对象（稠密矩阵或稀疏三元组两种输入） | 矩阵文件（+ metadata）→ `seurat.rds` |
| `scRNA_gene_expr.r` | 目标基因跨细胞类型表达可视化 | 已注释 rds → 均值 CSV + 5 类图 |
| `scRNA_gene_DE.r` | 样本 × 细胞类型表达矩阵 + 热图 / 柱状图 / 箱线图 | 已注释 rds → 表 + 图 |
| `scRNA_composition_sccoda.py` | 细胞组成成分差异分析（scCODA） | h5ad / 计数表 → 效应表 + 图 |

> `sankey_composition.r` 不是独立调用的分析脚本，而是 `scRNA_composition_sccoda.py` 的**必需配套脚本** —— 后者用 R + `ggalluvial` 渲染细胞组成桑基图，由 py 自动调用（无需手动运行），因此该文件**必须与 `scRNA_composition_sccoda.py` 放在同一目录**。缺 R 或 `ggalluvial` 时仅桑基图跳过，其余产物不受影响。

使用约定：

- 主流程各步以文件衔接，中间产物可单独复用；建议每步输出到独立目录（如 `01_QC/`、`02_deContamination/`），输出前缀（`01.` `02.` …）与步骤号对应，便于追溯。
- 数据准备/格式转换类脚本（`scRNA_spRNA_h5ad_preprocess.py` → `scRNA_h5ad2rds.r`）是从 h5ad 起步时的前置链：先清洗计数矩阵，再转成 Seurat rds，之后即可接主流程 step1；`scRNA_h5ad_seurat_convert.r` 负责单文件双向互转（含 rds → h5ad 导出）。
- 从 Broad Single Cell Portal（singlecell.broadinstitute.org）下载的数据可用 `scRNA_mtx2rds.r` 直接建对象：稠密矩阵走 `-E`，稀疏三元组走 `-M/-C/-G`，元数据里细胞名列用 `-CC` 指定（SCP 导出默认列名就是 `NAME`）。
- 独立脚本不依赖主流程的步骤顺序：`scRNA_gene_expr.r`、`scRNA_gene_DE.r` 需要一个含 celltype 列的已注释 Seurat 对象（典型来源是 step5 的 `seu_annotated.rds`）；`scRNA_composition_sccoda.py` 直接读 h5ad 或样本 × 细胞类型计数表，与主流程无关。
- 脚本会识别已存在的关键产物作为缓存（通常只重画图、不重算），强制重跑时删除对应 rds 即可。

## 可视化输出示例（`Images/`）

本目录下的 `Images/` 存放各脚本的示例可视化输出（每张图 PNG + PDF 同名成对），子目录与产出图的脚本一一对应（数据准备 / 格式转换类脚本不产图，因此没有对应子目录）：

| 子目录 | 对应脚本 | 内容 |
|---|---|---|
| `Images/01_QC/` | `scRNA_QC_step1.r` | QC 前后小提琴图 / 散点图 |
| `Images/02_deContamination/` | `scRNA_dedoublets_step2.r` | doublet UMAP、doublet/QC 指标小提琴图 |
| `Images/03_Nor_Harmony/` | `scRNA_Nor_harmony_step3.r` | PCA dimplot、PC elbow、HVG、Harmony dimplot |
| `Images/04_clustering/`、`Images/04_clustering_louvain/` | `scRNA_clustering_step4.r` | clustree、cluster / sample / group UMAP、各 cluster dotplot（两个目录为同一分析在 Leiden / Louvain 下的示例） |
| `Images/05_annotation/` | `scRNA_annotate_step5.r` | celltype UMAP、celltype × marker dotplot、`-S` 模式下每 cell type 一张 dotplot |
| `Images/gene_expr/` | `scRNA_gene_expr.r` | Ridge / Violin / Feature / Dot / Heatmap |
| `Images/gene_DE/` | `scRNA_gene_DE.r` | 热图、柱状图、箱线图 |
| `Images/composition_sccoda/` | `scRNA_composition_sccoda.py` | 组成桑基、相对丰度箱线图、离散度、效应森林图 |

各章「输出」小节末尾给出对应的示例目录。示例图仅用于展示输出形式与排版，其文件名与各章输出的 `<前缀><文件名>` 规则一致。

---

## 1 scRNA_QC_step1 —— 每样本独立 QC + 合并

### 功能

读取多个单细胞样本（10x 目录或 h5），对**每个样本独立执行 QC**（阈值取自该样本自身的 nFeature_RNA / nCount_RNA 分位数 + percent.mt / percent.HB 硬上限），QC 通过后合并为一个 Seurat 对象。基因级 min.cells 过滤在合并后统一执行，避免低频基因因单样本表达不足被过早删除。按样本自适应阈值，适配样本间测序深度差异大的数据集。归一化在后续步骤统一执行（先删基因后归一化，避免归一化后再删基因导致的分母不一致）。

### 输入

一个 CSV 文件（`-I`），三列，无表头要求（需含列名）：

| 列名 | 说明 |
|---|---|
| path | 每个样本的数据路径。`-DT 10X` 时为 10x 目录（含 matrix.mtx.gz / features.tsv.gz / barcodes.tsv.gz 的 `filtered_feature_bc_matrix/`）；`-DT h5` 时为单个 .h5 文件 |
| sample | 样本名（唯一，不能重复） |
| group | 分组名（如 control/treat），存入对象 metadata 的 group 列 |

示例 `input_list.csv`：

```csv
path,sample,group
/path/to/sample1/filtered_feature_bc_matrix,sample1,groupA
/path/to/sample2/filtered_feature_bc_matrix,sample2,groupB
```

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | input_list.csv 路径 |
| `-DT, --data_type` | 10X | 输入格式：`10X`（目录）或 `h5`（文件） |
| `-GC, --gene_column` | 2 | features.tsv 中基因名列号（默认 2） |
| `-qL, --quantile_low` | 0.01 | 每样本 nFeature/nCount 下分位数（如 0.01 = 1%） |
| `-qH, --quantile_high` | 0.99 | 每样本 nFeature/nCount 上分位数 |
| `-MaxM, --max_mt` | 10 | percent.mt 硬上限（%） |
| `-S, --species` | 无 | 启用红细胞过滤，取值 human / mouse / mammal。需与 -MaxHB 成对使用 |
| `-MaxHB, --max_HBs` | 无 | percent.HB（血红蛋白基因占比）上限（%），需与 -S 成对使用 |
| `-MC, --min_cells` | 3 | 合并后基因至少在多少细胞中检出才保留 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

QC 规则（每个样本独立计算）：
- 保留区间：nFeature_RNA ∈ [qL, qH] 分位数，且 nCount_RNA ∈ [qL, qH] 分位数
- 硬上限：percent.mt ≤ -MaxM；（启用时）percent.HB ≤ -MaxHB
- 合并后：基因需在 ≥ -MC 个细胞中检出

### 用法示例

```bash
# 基础：每样本分位数 QC + 合并，无红细胞过滤
./scRNA_QC_step1.r -I input_list.csv -DT 10X -GC 2 \
            -qL 0.01 -qH 0.99 -MaxM 10 -MC 3 \
            -P 01. -O 01_QC

# 人源样本 + 红细胞过滤（percent.HB ≤ 5%）
./scRNA_QC_step1.r -I input_list.csv -DT 10X -GC 2 \
            -qL 0.01 -qH 0.99 -MaxM 10 \
            -S human -MaxHB 5 -MC 3 \
            -P 01. -O 01_QC
```

### 输出

输出目录（-O）下生成：

| 文件 | 说明 |
|---|---|
| `<前缀>before_QC_vlnplot.png/.pdf` | QC 前小提琴图（2×2：nFeature_RNA / nCount_RNA / percent.mt / [percent.HB]） |
| `<前缀>before_QC_scatter.png/.pdf` | QC 前散点图（nCount×nFeature、nCount×percent.mt、nFeature×percent.mt、[nCount×percent.HB]） |
| `<前缀>after_QC_vlnplot.png/.pdf` | QC 后小提琴图 |
| `<前缀>after_QC_scatter.png/.pdf` | QC 后散点图 |
| `<前缀>QC_summary.txt` | 统计信息：参数、每样本阈值表、QC 前后细胞数、合并后细胞/基因数 |
| `<前缀>seu_afterQC.rds` | **主输出**：QC 后合并的 Seurat 对象（仅 counts 层，归一化在后续脚本执行，供下一步分析输入） |

示例输出：`Images/01_QC/`。

终端同时打印：每样本 `n=入 → 出 (ng=[lo,hi], tc=[lo,hi], mt<=X[, hb<=Y])` 统计、Per-sample QC thresholds 表、QC 前后每样本细胞数、合并后总数。

### 对象 metadata

输出 rds 的 `meta.data` 含：`orig.ident`（样本名）、`sample`、`group`、`percent.mt`、`percent.HB`（-S 时）。

### 注意事项

- 样本数过少（<10 细胞）时无法估计分位数，该样本全部保留（仅应用 mt/HB 上限），终端会提示 "too few cells, kept all"
- 单个样本输入时跳过合并步骤（不加细胞前缀）
- percent.mt 按 `^(?i)mt-` 匹配（兼容人 MT- / 小鼠 mt- 命名）
- 红细胞基因集：human = HBA1/HBA2/HBB/ALAS2/SLC4A1/GYPA；mouse = Hba-a1 等；mammal = ALAS2/FECH/SLC4A1/ANK1/SPTA1/SPTB/BPGM。基因实际不存在时警告并将 percent.HB 置 0

---

## 2 scRNA_dedoublets_step2 —— 去双细胞（per-sample）

### 功能

读取 step1 的 QC 后合并对象（仅 counts），内部执行轻量 Normalize → Scale → PCA → UMAP → FindCluster 作为双细胞建模的辅助，**这些归一化结果不保存**。然后用 **scDblFinder per-sample 建模**（通过 samples 参数）预测每个细胞为 Doublet/Singlet。可选启用 **DecontX per-sample 跑**（每样本独立去环境 RNA 污染），同样按每样本自己的 90% 分位定 L-prob / H-prob 阈值。

最终过滤 Singlet（-D 时还要求 L-prob），用 `DietSeurat` 只保留 counts 层，输出供下一步归一化使用。

**断点重画机制**：脚本会先检查输出目录中 `seu_doublets.rds` 与 `seu_clean.rds` 是否存在：
- 两个 rds 均不存在 → 跑全流程（轻量预处理 + DecontX + scDblFinder + 保存两个 rds + 画图）
- 仅 `seu_doublets.rds` 存在（上次中断）→ 读取缓存，跳过耗时计算，画图并重建 `seu_clean.rds`
- 两个 rds 均存在 → 读取缓存，**只重画图**（配合修改 -W/-H），不重复保存 rds；过滤统计仍会写入 summary

重画图只需几十秒，无需重跑 ~1 小时的 scDblFinder。若需强制重跑分析，删除输出目录中的 `seu_doublets.rds` 即可。

### 输入

`-I` 指向 step1 输出的 `.rds` 文件（counts-only 的合并 Seurat 对象，meta.data 含 `sample` 列）。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | step1 输出的 rds 路径 |
| `-D, --decontx` | false | 启用 DecontX per-sample 去环境 RNA（不传则关闭） |
| `-PC` | 30 | PCA 维度数 |
| `-R` | 0.1 | FindClusters 分辨率（双细胞建模用，0.1 即可） |
| `-CA, --cluster_algorithm` | Louvain | FindClusters 算法；Leiden 需 `-P_env` 指向装了 leidenalg 的 Python |
| `-P_env` | 无 | Leiden 算法用的 Python 环境路径 |
| `-DR, --double_rate` | auto | 期望双细胞率，传给 scDblFinder。不传则按细胞数自动估 |
| `-SC, --sample_column` | sample | meta.data 中样本标识列名（per-sample 建模依据） |
| `-T` | 10 | scDblFinder 线程数 |
| `-W, --width` | 16 | 图宽（英寸） |
| `-H, --height` | 8 | 图高（英寸） |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 默认（无 DecontX）
./scRNA_dedoublets_step2.r -I 01_QC/01.seu_afterQC.rds \
    -PC 30 -R 0.1 -T 20 -P 02. -O 02_deContamination

# 开启 DecontX（per-sample 去环境 RNA）
./scRNA_dedoublets_step2.r -I 01_QC/01.seu_afterQC.rds \
    -D -PC 30 -R 0.1 -T 20 -P 02. -O 02_deContamination
```

### 输出

输出目录 `-O` 下生成：

| 文件 | 说明 |
|---|---|
| `<前缀>seu_doublets.rds` | 含 scDblFinder 标记（Double_score、Is_Double）、可选 DecontX 标记（decontX_class）——未过滤的中间对象（仅 counts 层） |
| `<前缀>seu_clean.rds` | **主输出**：过滤 Singlet（-D 时还要求 L-prob）后的 Seurat 对象，仅 counts 层，供 step3 归一化输入 |
| `<前缀>decontX_umap.png/.pdf` | 仅 -D 时生成：DecontX contamination 与 L/H-prob 分类 |
| `<前缀>decontX_violin.png/.pdf` | 仅 -D 时生成：L/H-prob 组的 nCount/nFeature/percent.mt 小提琴图 |
| `<前缀>doublet_umap.png/.pdf` | Doublet score 与 Singlet/Doublet 分类 |
| `<前缀>doublet_violin.png/.pdf` | Singlet/Doublet 组的 QC 指标小提琴图 |
| `<前缀>Dedoublets_summary.txt` | 统计信息：参数、每样本 Doublet/Singlet 计数、过滤前后每样本细胞数 |

示例输出：`Images/02_deContamination/`（示例为未启用 `-D` 时的输出）。

### 注意事项

- step1 输出的 rds 必须含 `sample` 列（默认即有，可通过 `-SC` 修改）
- DecontX 必须 **per-sample 跑**（混样本会引入偏差），本脚本按样本拆分执行；H-prob 阈值也按每样本自己的 90% 分位定
- 脚本输出的 rds **仅保留 counts 层**（PCA/UMAP/聚类等中间结果丢弃），归一化留给 step3，避免重复计算
- `seu_clean.rds` 是下游输入，`seu_doublets.rds` 留作可追溯性
- Leiden 算法需要额外配置 Python 环境（`-P_env`），简单分析用默认 Louvain 即可
- scDblFinder 不传 `-DR` 会按细胞数自动估双细胞率（约 0.8%/千细胞）；知道真实多细胞率时可手动传

---

## 3 scRNA_Nor_harmony_step3 —— 归一化与批次整合（Harmony）

### 功能

读取 step2 的 counts-only 对象，按顺序执行：**NormalizeData（LogNormalize）→ FindVariableFeatures（vst）→ ScaleData → RunPCA → RunHarmony** 批次校正。输出含 `pca` 与 `harmony` 两个降维的结果对象，供 step4 聚类使用。

step1/2 刻意不归一化（避免删基因前归一化导致分母不一致），归一化在本步统一完成。

### 输入

`-I` 指向 step2 输出的 `seu_clean.rds`（counts-only，含 `sample` 列）。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | step2 输出的 rds 路径 |
| `-BC, --batch_col` | sample | Harmony 批次变量列。**必须用样本列（sample），不要用 group**——否则会把疾病分组的生物学差异当作批次抹掉 |
| `-NF, --n_features` | 2000 | vst 高变基因数 |
| `-HT, --harmony_theta` | 2 | Harmony 多样性聚类惩罚。越大混匀越强；0 关闭惩罚 |
| `-PC` | 30 | PCA 维度数 |
| `-T` | 10 | Harmony 线程数 |
| `-J, --doJackStraw` | false | 额外跑 JackStraw 分析 PCA 显著性 |
| `-W, --width` | 12 | 图宽（英寸） |
| `-H, --height` | 8 | 图高（英寸） |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_Nor_harmony_step3.r -I 02_deContamination/02.seu_clean.rds \
    -BC sample -NF 2000 -HT 2 -PC 30 -T 20 \
    -P 03. -O 03_Nor_Harmony
```

同 step2，若 `seu_harmony.rds` 已存在则跳过计算只重画图（改 -W/-H 时用）；强制重跑需删除该 rds。

### 输出

输出目录 `-O` 下生成：

| 文件 | 说明 |
|---|---|
| `<前缀>seu_harmony.rds` | **主输出**：归一化 + 校正后对象（counts/data/scale.data + pca/harmony reductions + HVG 信息），供 step4 聚类输入 |
| `<前缀>pca_dimplot.png/.pdf` | PCA 按样本着色（校正前状态，看批次效应） |
| `<前缀>pca_elbowplot.png/.pdf` | PC 方差贡献（选 PC 数参考） |
| `<前缀>HVGs.png/.pdf` | 高变基因分布（标注 top10） |
| `<前缀>harmony_dimplot.png/.pdf` | Harmony 按样本着色（校正后状态，样本应混匀） |
| `<前缀>Nor_harmony_summary.txt` | 统计信息：参数、每批次细胞数 |

示例输出：`Images/03_Nor_Harmony/`。

### 注意事项

- **批次变量是 sample 不是 group**：本脚本对 `-BC group` 直接报错拦截。用 group 做 Harmony 会把分组间的生物学差异校正掉，下游比较将失效
- 归一化必须先于删基因或与删基因顺序一致——step1 已先删基因（min.cells 过滤）再输出 counts，本步归一化顺序正确
- HVG 合并后统一选择（vst）会受大样本轻微影响；Harmony 的批次校正可部分补偿，如样本量极不平衡可考虑 per-sample HVG
- 后续 step4 聚类将基于 `harmony` reduction 而非 `pca`

---

## 4 scRNA_clustering_step4 —— 聚类与差异表达

### 功能

读取 step3 输出的 harmony 对象，分两阶段产出：
1. **探索模式（不传 -Res）**：FindNeighbors + 多分辨率 FindClusters（`seq(0, 1, by=0.1)`，Louvain；`seq(0.1, 1, 0.1)`，Leiden）→ 生成 clustree.png（用于人工选 Res）
2. **聚类模式（传 -Res）**：在选定的 Res 下聚类 + RunUMAP + FindAllMarkers DE → 生成 `seu_res<Res>_clustered.rds` 与 dotplot 集合

支持 `-CA Louvain`（默认）/ `-CA Leiden -P_env <python>` 两种聚类算法。

### 双产物缓存 + 完整性检查

脚本按两个独立产物判断分支，并对最终 rds 做完整性检查（大小 + readRDS 试读）：

| 状态 | 行为 | 提示 |
|---|---|---|
| `seu_res<Res>_clustered.rds` 完整存在 | 跳过所有计算，只重画图 | "File is complete. Skipping... re-plotting only." |
| rds 存在但**不完整**（<1MB 或 readRDS 失败）| 警告 + 自动全流程重跑（覆盖半成品）| "[Warning] Found ... but it is too small / readRDS failed - will re-run" |
| 无 -Res + clustree.png 不存在 | FindNeighbors + 多分辨率 + 画 clustree + 退出 | "No -Res provided: exploration mode only" |
| 无 -Res + clustree.png 存在 | 跳过 clustree，提示看图选 Res | "clustree.png already exists, skipping multi-resolution exploration" |
| 有 -Res + rds 不存在 | 跳过 clustree（如已有）→ 单分辨率聚类 + UMAP + 存 rds | 常规流程提示 |
| 缓存重画 + markers csv 已存在 | 跳过 DE 重算 | "[Info] Cached run and cluster markers csv already exists - skipping DE" |

**强制重跑**：删除 `seu_res<Res>_clustered.rds` 即可。

### 输入

`-I` 指向 step3 的 `seu_harmony.rds`（含 harmony reduction）。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | step3 输出 rds 路径 |
| `-PC` | 30 | PCA 维度数（聚类时基于哪个 reduction 取前 N 维） |
| `-Red` | harmony | 聚类所用的 reduction（pca 或 harmony） |
| `-Res` | 无 | 聚类分辨率。**不传 = 仅探索模式**，传了才生成最终 rds |
| `-CA` | Louvain | 聚类算法：Louvain（R 原生）或 Leiden（需 `-P_env`） |
| `-P_env` | 无 | Leiden 用的 Python 环境路径 |
| `-A` | RNA | assay 名 |
| `-skipDE` | false | 跳过 FindAllMarkers（不跑 DE，省时间） |
| `-T` | 10 | 线程数（用于 FindAllMarkers 的 future 并行） |
| `-W, --width` | 10 | UMAP 图宽（英寸） |
| `-H, --height` | 8 | UMAP 图高（英寸） |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 阶段 1: 多分辨率探索（不传 -Res）
./scRNA_clustering_step4.r -I 03_Nor_Harmony/03.seu_harmony.rds \
    -Red harmony -PC 30 -T 20 -P 04. -O 04_clustering \
    -CA Leiden -P_env /your_path/env/python

# 阶段 2: 选定 Res（如 0.7）后聚类 + UMAP + DE
./scRNA_clustering_step4.r -I 03_Nor_Harmony/03.seu_harmony.rds \
    -Red harmony -PC 30 -Res 0.7 -T 20 -P 04. -O 04_clustering \
    -CA Leiden -P_env /your_path/env/python

# 重画图：rds 完整存在时只重画，不重算
./scRNA_clustering_step4.r -I 03_Nor_Harmony/03.seu_harmony.rds \
    -Red harmony -PC 30 -Res 0.7 -W 12 -H 8 -P 04. -O 04_clustering \
    -CA Leiden -P_env /your_path/env/python
```

### 输出

输出目录 `-O` 下生成：

| 文件 | 说明 |
|---|---|
| `clustree.png/.pdf` | 探索图（多分辨率聚类树，仅探索模式 / 阶段 1 生成）|
| `<前缀>seu_res<Res>_clustered.rds` | **主输出**：含 UMAP、cluster 标记，可直接用于注释脚本 |
| `<前缀>res<Res>_cluster_umap.png/.pdf` | UMAP 按 cluster 着色 |
| `<前缀>res<Res>_sample_umap.png/.pdf` | UMAP 按 sample 着色（看整合效果） |
| `<前缀>res<Res>_group_umap.png/.pdf` | UMAP 按 group 着色 |
| `<前缀>clusters_DE/<前缀>res<Res>_cluster_markers.csv` | 每个 cluster 的差异表达结果（全部 marker） |
| `<前缀>clusters_DE/<前缀>res<Res>_cluster_top10_markers.csv` | 每个 cluster top 10 marker |
| `<前缀>clusters_DE/<前缀>res<Res>_cluster_<N>_dotplot.png/.pdf` | 每个 cluster 的 top 10 marker dotplot |
| `<前缀>Clustering_summary.txt` | 统计信息：参数、cluster 细胞数、每样本细胞数 |

示例输出：`Images/04_clustering/`（Leiden）、`Images/04_clustering_louvain/`（Louvain），均为 res 0.1。

### 注意事项

- **聚类算法选 Leiden 还是 Louvain**：本流程统一选 Leiden（leidenalg + igraph 已配好），保证连通社区
- **clustree 是必看图**：先探索、再聚类，跳过阶段 1 直接传 Res 是盲选
- **同一 Res 重复跑会读缓存**：rds 完整时几秒出图（仅 ggsave）；强制重跑需删 rds
- **DE 是耗时大头**：全数据集（数万细胞 × 十几个 cluster）跑 `FindAllMarkers` 通常需十几分钟；后续注释脚本不再重跑 DE，直接用 `clusters_DE/` 里的 csv
- **`-T` 只对 FindAllMarkers 生效**：FindNeighbors/FindClusters/RunUMAP 自身单线程（Seurat 限制）

---

## 5 scRNA_annotate_step5 —— 细胞类型注释

### 功能

读取 step4 的聚类对象（`seu_res<Res>_clustered.rds`），用两种模式：
1. **探索模式（-MB）**：用户提供 cell type marker 列表（2 列 TSV），脚本画 cluster × marker dotplot 用于人工判断每个 cluster 的身份。**不写 rds。**
2. **标注模式（-MA）**：用户提供 cluster → cell type 映射表（3 列 TSV：Cluster / Markers / CellType），脚本画 celltype UMAP + celltype dotplot，并写 `seu_annotated.rds`（含 celltype meta 列）。

**两种模式互斥**，每次只能传一种。**不重算 DE**——复用 step4 的 `clusters_DE/`。

### 输入

`-I` 指向 step4 的 `seu_res<Res>_clustered.rds`。
`-MB` 或 `-MA` 之一指向 marker 列表/标注表 TSV。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | step4 输出 rds 路径 |
| `-MB, --markers_before` | 无 | marker 列表 TSV（探索模式）。格式：V1=celltype, V2=markers（逗号分隔）|
| `-MA, --markers_after` | 无 | 标注表 TSV（标注模式）。格式：Cluster, Markers, CellType |
| `-S, --split` | false | -MB 模式下每 cell type 单独画 dotplot（marker 多时用）|
| `-CO, --celltype_order` | 无 | dotplot 中 cell type 顺序（逗号/空格分隔，如 `B_cell, T_cell, NK`）|
| `-W, --width` | 24 | dotplot 图宽（英寸，UMAP 硬编码 12×12）|
| `-H, --height` | 14 | dotplot 图高（英寸，UMAP 硬编码 12×12）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### rds 缓存 + 完整性检查

脚本对 `seu_annotated.rds` 做完整性检查（大小 + readRDS 试读），完整存在时：

| 模式 | 行为 |
|---|---|
| -MB | 仍跑：只画探索 dotplot，不写 rds（无副作用）|
| -MA | **仍读 TSV**（按用户选项 B 设计，便于改映射后重跑）：重新映射 cluster → celltype、画 UMAP/dotplot，**覆盖写 rds** |

不完整（< 1MB 或 readRDS 失败）→ 警告 + 全流程重跑。

### 用法示例

```bash
# 阶段 1：探索（看每个 cluster 的候选身份）
./scRNA_annotate_step5.r -I 04_clustering/04.seu_res0.1_clustered.rds \
    -MB Markers.txt -P 05. -O 05_annotation

# 阶段 2：基于探索结果人工决策后，最终标注
./scRNA_annotate_step5.r -I 04_clustering/04.seu_res0.1_clustered.rds \
    -MA 05_annotation/annotation.txt -P 05. -O 05_annotation

# 重画图：rds 完整时只重画（不重映射）
./scRNA_annotate_step5.r -I 04_clustering/04.seu_res0.1_clustered.rds \
    -MA 05_annotation/annotation.txt -W 24 -H 14 -P 05. -O 05_annotation
```

### 输出

输出目录 `-O` 下生成：

| 文件 | 说明 |
|---|---|
| `<前缀>cluster_markers.png/.pdf` | 仅 -MB 模式：cluster × marker dotplot（一张图，按 cell type 分面）|
| `<前缀>markers_dotplot/<前缀><celltype>_cluster_markers.png/.pdf` | 仅 -MB + -S 模式：每 cell type 一张 dotplot |
| `<前缀>umap_celltype.png/.pdf` | 仅 -MA 模式：UMAP 按 cell type 着色 |
| `<前缀>dotplot_celltype.png/.pdf` | 仅 -MA 模式：celltype × marker dotplot |
| `<前缀>seu_annotated.rds` | 仅 -MA 模式：含 celltype meta 列的主输出（下游分析输入）|
| `<前缀>Annotation_summary.txt` | 统计信息：参数、每 cell type 细胞数（仅 -MA）|

示例输出：`Images/05_annotation/`（含 `-MB -S` 模式下每 cell type 一张 dotplot）。

### 注意事项

- **两模式互斥**：传 -MB 和 -MA 会报错
- **未映射 cluster 会变 NA**：若 -MA 文件里没列出某 cluster，脚本会警告"cluster X is NOT in -MA file and were set to NA"，需补齐映射
- **-S 仅用于 -MB**：与 -MA 一起传会报错
- **rds 缓存行为**：完整 rds + -MB 仅重画；完整 rds + -MA 重映射并覆盖 rds（选项 B）；不完整 rds → 警告 + 重跑
- **不重算 DE**：脚本完全复用 step4 的 clusters_DE/，自己不再跑 FindAllMarkers
- **celltype_order 透传**：用 `-CO` 指定 dotplot 中 cell type 的固定显示顺序，便于多次跑结果一致
- **UMAP 硬编码 12×12**：与 dotplot 分开，因 14+ celltype 图例在右侧占大，UMAP 不需要 -W/-H 默认宽度；dotplot 仍用 -W/-H
- **-MA 文件格式**：3 列带表头（Cluster / Markers / CellType），Markers 列内**必须统一逗号**（分号会被当作基因名的一部分）；Markers 值可用 `"..."` 包裹
- **marker 列表之外的细胞类型**：经典 marker 可能提示出输入列表以外的细胞类型，应在结果中说明其判定依据

---

## 6 scRNA_gene_expr（独立脚本）—— 目标基因跨细胞类型表达可视化

### 功能

读取已注释的 Seurat 对象（典型来源是 step5 的 `seu_annotated.rds`），计算每个目标基因在每种细胞类型中的平均表达量（基于 data 层，LogNormalize 后），导出为 CSV，并基于 Seurat 内置函数生成 5 种可视化图（Ridge / Violin / Feature / Dot / Heatmap）。

### 输入

- `-I` 指向已注释的 Seurat rds（如 step5 的 `seu_annotated.rds`；必须含 `data` 层、含 celltype meta 列）
- `-G` 目标基因列表：逗号分隔字符串（如 `GENE1,GENE2,GENE3`），或 TSV/TXT 文件路径（一行一个基因）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | 注释后的 Seurat rds 路径 |
| `-G, --genes` | （必填） | 基因列表：逗号字符串 或 文件路径 |
| `-CT, --celltype_col` | celltype | 用于分组的 meta 列名 |
| `-CTO, --celltype_order` | 按列默认排序 | celltype 在图中的固定顺序（逗号/空格分隔，如 `T,B,Macro`） |
| `-RI, --ridge_cols` | 0 | RidgePlot 列数（0 = 默认 4 列；基因少时可用 2-3） |
| `-DS, --downsample` | 100 | 热图每组（celltype）最多取多少细胞 |
| `-W, --width` | 12 | 图宽（英寸） |
| `-H, --height` | 10 | 图高（英寸） |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 方式 1：命令行直接给基因列表
./scRNA_gene_expr.r -I 05_annotation/05.seu_annotated.rds \
    -G GENE1,GENE2,GENE3 \
    -CTO celltypeA,celltypeB,celltypeC \
    -W 16 -H 12 -P 06. -O 06_gene_expression

# 方式 2：从文件读基因列表
./scRNA_gene_expr.r -I 05_annotation/05.seu_annotated.rds \
    -G genes.txt -P 06. -O 06_gene_expression
```

### 输出

输出目录 `-O` 下生成：

| 文件 | 说明 |
|---|---|
| `<前缀>avg_expression.csv` | **主输出**：gene × celltype 平均表达矩阵（基于 data 层） |
| `<前缀>ridgeplot.png/.pdf` | RidgePlot — 表达分布山脊图 |
| `<前缀>violinplot.png/.pdf` | VlnPlot — 小提琴图 |
| `<前缀>featureplot.png/.pdf` | FeaturePlot — UMAP 上表达分布（需 rds 含 umap）|
| `<前缀>dotplot.png/.pdf` | DotPlot — 表达比例 × 平均表达气泡图 |
| `<前缀>heatmap.png/.pdf` | DoHeatmap — 热图（用 scale.data，跨基因 z-score）|
| `<前缀>GeneExpression_summary.txt` | 参数 + 基因存在性 + 缺失列表 |

示例输出：`Images/gene_expr/`。

### 注意事项

- **依赖 data 层**：rds 必须已 LogNormalize（seu_annotated.rds 默认满足）；无 data 层则表达值为零，前4 张图无意义
- **缺失基因**：输入基因中不在 Seurat 对象里的会被警告 +不参与绘图（不会让脚本崩溃）
- **FeaturePlot 需 UMAP**：输入 rds 必须有 `umap` reduction（由聚类步骤产生）；无则自动跳过并打印 "[Info] No 'umap' reduction..."
- **数据层差异**：CSV + Ridge/Vln/Feature/Dot 用 **data 层**（LogNormalize，跨 celltype 可比）；Heatmap 用 **scale.data**（z-score，仅看每个 celltype 内的相对高低）
- **Heatmap 为什么要先 ScaleData**：DoHeatmap 默认读 scale.data 层，而 scale.data 通常只含高变基因（HVG）——非 HVG 的目标基因会被静默省略（表现为"只显示几个基因"）。脚本在画热图前对目标基因单独 `ScaleData(features=目标基因)`，确保**全部输入基因都出现**
- **Heatmap 顶部色块 + 45° 名字**：`group.bar = TRUE` 显示 celltype 颜色块，`label = TRUE`（Seurat 默认 `angle = 45`）显示细胞类型名字，斜向上 45° 排列
- **Heatmap 按组采样**：`subset(seu, downsample = -DS)`（默认每组最多 100 细胞）避免 5 万+ 细胞渲染过慢；因已采样，celltype 用默认顺序（不再按细胞数排序）。`-DS` 调大更接近真实分布，调小渲染更快
- **FeaturePlot 子图 1:1**：每个子图 `aspect.ratio = 1`；列数按基因数自动选，使整体宽高比最接近 `-W/-H`。7 个基因 → 3×3（比例约 0.9）；9-12 个基因 → 4×3（约 1.2）
- **全图尺寸随基因数自适应**（脚本自动算，不需手调）：
  - **RidgePlot**：4 列默认；高度按 `行数 × celltype × 0.22 + 1.2` 算，但 `≤ -W × 2`（避免长条）
  - **ViolinPlot**：自动选 1-8 列使整体比例接近 `-W/-H`；1:1 子图布局，高度按 `行数 × unit_w` 算，同样有 `≤ -W × 2` 上限
  - **DotPlot**：宽度按 `基因 × 0.55 + 2` 扩展，高度按 `celltype × 0.32 + 1.5` 算；基因数 ≤ ~15 时与 `-W -H` 一致
  - **Heatmap**：宽度按 `基因 × 0.45 + 1.5`，高度按 `celltype × 0.45 + 2`，同样有 `≤ -W × 1.5` 上下限
  - 所有图都有 `-W × 1.5` 或 `-W × 2` 高度上限保护，避免 30+ 基因变超长纸带
- **RidgePlot 列数**：`-RI` 控制（0 = 默认 4 列），图例已隐藏（x 轴即细胞类型）。注意 Seurat 的参数名是 **`ncol`**——`cols` 是颜色向量（误传会报 "Insufficient values in manual scale"）
- **celltype 顺序**：默认沿用对象内的水平顺序；`-CTO` 可指定固定顺序（如 `-CTO T,B,Macro`）

---

## 7 scRNA_gene_DE（独立脚本）—— 样本 × 细胞类型表达矩阵 + 热图 + 柱状图 + 箱线图

### 功能

读取已注释的 Seurat 对象（如 step5 的 `seu_annotated.rds`），对**单个目标基因**计算其在每种"样本 × 细胞类型"组合下的平均表达量，并一次性产出数值表与三类图：

1. **样本 × 细胞类型均值表**（CSV）—— 缺失组合填 0，供下游统计软件直接读
2. **热图** —— 该表的按-celltype z-score 版本，看样本间相对高低
3. **柱状图**（每个指定 celltype 一张）—— 绝对平均表达量，看样本绝对值分布
4. **箱线图**（每个指定 celltype 一张）—— 按分组（如 对照 vs 处理）比较，**每个点是一个样本**（不是细胞），附 Wilcoxon p 值与显著性横线

### 输入

- `-I` 已注释的 Seurat rds（须含样本列、celltype 列、分组列与 `data` 层）
- `-G` 单个基因符号
- `-CT` 要画柱状图/箱线图的细胞类型（逗号分隔）；留空则只输出 CSV + 热图

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | 已注释的 Seurat rds 路径 |
| `-G, --gene` | （必填） | **单个**基因符号（含逗号/空格会被拦截） |
| `-S, --sample_col` | sample | 样本 meta 列名 |
| `-CTC, --celltype_col` | celltype | 细胞类型 meta 列名 |
| `-GC, --group_col` | group | 分组 meta 列名（仅供箱线图） |
| `-L, --layer` | data | 用于取值的 Seurat 层（默认 LogNormalize 后） |
| `-CT, --celltypes` | 空 | 要画柱状图 + 箱线图的细胞类型（逗号分隔） |
| `-SO, --sample_order` | 空 | **柱状图**的样本顺序（左→右）；留空用表内顺序。不影响热图/箱线图 |
| `-GO, --group_order` | 空 | **箱线图**的分组顺序（左→右）；留空按字典序。首个为参照组 |
| `--alt` | auto | 箱线图检验方向：`auto`（等效 `less`，即右侧组 > 左侧组）/ `two.sided` / `greater` / `less` |
| `--signif_bar` | auto | 是否画箱线图显著性横线；`none` 关闭 |
| `--bar_gap` | 0.5 | 箱线图横线距最高数据点的间距（y 轴 padding 的倍数），越大越远 |
| `--col_order` | **hclust** | **热图列顺序**：`hclust`（默认）= 按样本表达谱聚类（列 z-score 后 Euclidean 距离 + complete linkage），让表达模式相近的细胞类型相邻；`csv` = 保持表内顺序 |
| `--box_coef` | 1.5 | 箱线图须线长度（IQR 的倍数）。默认 1.5 = 经典 Tukey 规则；超出围栏的样本画成**空心黑圆**（且从散点层剔除，不重复画）。传 `1e9` 可让须线覆盖真实极值、不判离群点。**注意 ggplot 的 `coef=0` 不等于 base R 的 `range=0`**（详见注意事项） |
| `-W / -H` | 14 / 10 | **热图**宽 / 高（英寸） |
| `-BW / -BH` | 8 / 6 | **柱状图**宽 / 高（英寸） |
| `-QW / -QH` | 5 / 6 | **箱线图**宽 / 高（英寸） |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `--no-plot` | FALSE | 跳过全部可视化（仅 CSV + summary） |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 完整输出：CSV + 热图 + 2 张柱状图 + 2 张箱线图 + 统计表
./scRNA_gene_DE.r -I 05_annotation/05.seu_annotated.rds -G <gene> \
    -CT celltypeA,celltypeB \
    -SO sample1,sample2,sample3 \
    -GO groupB,groupA \
    -O 07_geneDE

# 只要数值表 + 热图
./scRNA_gene_DE.r -I 05_annotation/05.seu_annotated.rds -G <gene> -O 07_geneDE

# 只要数值表
./scRNA_gene_DE.r -I 05_annotation/05.seu_annotated.rds -G <gene> --no-plot -O 07_geneDE
```

### 输出

输出目录 `-O` 下生成：

| 文件 | 说明 |
|---|---|
| `<前缀><基因名>_sample_celltype_avg.csv` | **数值表** — N_sample × N_celltype 均值矩阵（缺失组合填 0） |
| `<前缀><基因名>_heatmap.png/.pdf` | **热图** — sample × celltype，按 celltype 列 z-score，NAVY → WHITE → FIREBRICK，范围 [-3, 3]，图例标题为 **Expression**（注意：图例文字为 Expression，但色阶实际是 z-score 值）|
| `<前缀><基因名>_<celltype>_barplot.png/.pdf` | **柱状图** — 该 celltype 各样本的绝对平均表达量（每 `-CT` 一项一张） |
| `<前缀><基因名>_<celltype>_boxplot.png/.pdf` | **箱线图** — 该 celltype 两组的样本级比较 + Wilcoxon p + 显著性横线 |
| `<前缀><基因名>_boxplot_stats.csv` | **检验结果表** — 每个 celltype 一行：分组细胞数、中位数、均值、检验名、p 值 |
| `<前缀>GeneDE_summary.txt` | 参数与统计汇总 |

示例输出：`Images/gene_DE/`。

### 注意事项

- **单基因输入**：`-G` 必须是单个基因符号；要批量分析请循环调用脚本
- **缺失组合填 0**：某 sample × celltype 无细胞时对应值为 0。注意这与"真实零表达"不可区分——箱线图中若某样本该 celltype 无细胞，其 0 会参与检验，解读时须结合每组细胞数
- **箱线图的点 = 样本，不是细胞**：用细胞做重复会产生伪重复（p 值被细胞数撑到 1e-16）。每组样本数少时要留意检验的分辨率下限——如 4 vs 3 的两侧 Wilcoxon 最小可能 p 为 2/C(7,4) = 0.0571，无论差异多大都无法 < 0.05
- **箱线图显著性横线**：横线锚定在"最高数据点之上"（`--bar_gap` 控制距离），否则 y 轴上限会裁掉最高样本；两端小竖线朝下。若日志出现 `Removed N row ... (stat_boxplot())` 说明有样本被裁掉
- **箱线图用经典 Tukey 规则（`--box_coef 1.5`）**：须线终点是"围栏（Q3+1.5×IQR）内**最远的那个数据点**"，**不是围栏值本身**。超出围栏的样本画成**空心黑圆**并**从散点层剔除**（避免同一样本画两次）。
- **箱线图须线可能为 0（小样本陷阱）**：ggplot2 的 `geom_boxplot` 用 **type-7 分位数**算 Q1/Q3。当 n 很小（如 4）且存在远点时，Q3 会被拉到高于所有非离群值，导致上须线长度为负 → 画不出来（该样本改由空心圆标出）。例如某组 4 个值为 `10, 20, 30, 100`：type-7 分位数给出 Q1 = 17.5、Q3 = 47.5、IQR = 30、围栏上界 = 92.5；最大值 100 超出围栏被判为离群点，而围栏内最大的非离群值（30）**小于** Q3（47.5）→ 上须线长度为负，无法绘制。
  - **注意**：ggplot 的 `coef=0` 与 base R 的 `range=0` **语义不同**——ggplot 的 `coef=0` 令围栏等于 Q1/Q3，所有超出箱体的点都被判为离群点，须线退化成 0 长度（实测确认），**不要用 0**。
  - base R 的 `boxplot()` 用 Tukey hinge（`fivenum`）而非 type-7 分位数，**同一组数据**会给出更长的须线且不判离群点——**两个实现的结论不一致是正常的**。
  - 若希望须线覆盖真实极值、不判离群点，传 `--box_coef 1e9`。
- **`outlier.*` 参数在该 ggplot2 版本下不生效**：`geom_boxplot(outlier.shape=1, outlier.size=2.8)` 实测既不改变 `ggplot_build` 里的 `shape`（仍为 19）也不渲染任何标记。脚本改为 `outlier.shape = NA` + 自己用 `geom_point(shape=1, size=2.8, colour="black", stroke=0.9)` 画离群点，保证一定出现。
- **热图列聚类（默认开启，`--col_order hclust`）**：先对“样本 × 细胞类型”矩阵按列 z-score，再对细胞类型做层次聚类（Euclidean + complete），把表达模式相近的细胞类型排到一起。**必须在列 z-score 之后聚类**——否则表达量绝对水平最高的列会主导距离量度。聚类只影响热图列顺序，**CSV 列顺序与柱状图/箱线图均不受影响**；若要保持表内顺序传 `--col_order csv`。
- **热图 z-score 按列（celltype）做**：让跨 celltype 颜色可比；"全列偏蓝"表示该 celltype 在所有样本中表达均低
- **`-SO` 只影响柱状图**，`-GO` 只影响箱线图；热图始终用表内固定顺序
- **列顺序依赖 factor level**：celltype 列顺序取自对象 `celltype` 列的 factor level 顺序（非字母序），改注释后顺序可能变化
- **图可复现**：脚本在绘图前固定随机种子（`set.seed(42)`），jitter 散点位置每次一致；PNG 逐字节可复现，PDF 仅嵌入时间戳不同

## 8 scRNA_composition_sccoda（独立脚本）—— 细胞组成成分差异分析（scCODA）

### 简介

scCODA 是基于 Dirichlet-Multinomial 模型的贝叶斯成分差异分析方法。给定一个**样本 × 细胞类型**的计数表，识别哪些细胞类型在两组之间的**相对丰度**显著不同（同时避免"一类升高导致另一类机械性下降"的成分数据假阳性）。

### 标准调用

```
./scRNA_composition_sccoda.py -I <input.h5ad> \
    -S sample -CT celltype -GC group -R <reference> -GR <ref_level> \
    --deterministic \
    -O <output_dir> -P <prefix>.
```

`<reference>` 取 `automatic`（默认）或某个 cell type 名称；`<ref_level>` 是该分组列的参考水平。

**建议正式分析一律带上 `--deterministic`**——否则同样的命令重复运行会得到略有差异的结果（只有设了种子也不够，原因见下）；带上它可保证结果**逐位可复现**、报告引用的数字永远对得上，代价几乎为零。桑基图（Plot 1）由同目录的 `sankey_composition.r`（**本脚本的必需配套脚本**，须与 py 同目录）渲染，脚本会自动调用，无需手动运行。

### 核心参数速查

| 参数 | 默认 | 说明 |
|---|---|---|
| `-I/--input` | 必填 | 输入文件：AnnData `.h5ad`、样本 × 细胞类型计数表、每行一细胞的 obs 表（见 `-T`）|
| `-T/--input_type` | auto | `h5ad` / `counts` / `obs`；`auto` 按扩展名与列名判定 |
| `-S/--sample_col` | sample | h5ad 或 obs 模式下的样本列名 |
| `-CT/--celltype_col` | celltype | h5ad 或 obs 模式下的细胞类型列名 |
| `-GC/--group_col` | group | 协变量列：分组定义所在列 |
| `-F/--formula` | 空 | patsy 公式；与 `-GC/-GR` 二选一；显式优先级最高 |
| `-GR/--group_ref` | 空 | 参考组水平；给出后自动生成 `C(col, Treatment('X'))` |
| `-EC/--extra_covariates` | 空 | 额外协变量，逗号分隔 |
| `-R/--reference_cell_type` | automatic | 参考细胞类型；automatic 自动选"出现率 ≥ `1 - -RA`（默认 95%）且相对丰度离散度最低" |
| `-RA/--ref_absence_threshold` | 0.05 | automatic 的缺失阈值 |
| `-M/--method` | hmc | 采样器：`hmc` / `nuts` / `hmc_da` |
| `-N/--num_results` | 20000 | 后验采样数 |
| `-B/--num_burnin` | 5000 | burn-in 步数 |
| `-L/--num_leapfrog_steps` | 10 | HMC leapfrog 步数 |
| `-TD/--max_tree_depth` | 10 | NUTS 最大树深 |
| `-SS/--step_size` | 0.01 | HMC/NUTS 步长 |
| `-SD/--seed` | 42 | 随机种子。**注意**：仅设种子不足以让结果逐位可复现（见 `--deterministic`）|
| `--deterministic` | 关 | **让 MCMC 逐位可复现**。scCODA 有三处会导致每次运行结果漂移：① 构造 `CompositionalAnalysis` 时用 `tf.random.normal` 抽 MCMC 初始值；② 调用 `tfp.mcmc.sample_chain` 时**没传 seed**；③ TF 的 CPU 浮点算子默认非确定性。本开关会：在 TF 导入前设 `TF_DETERMINISTIC_OPS=1`、**在模型构造前**设 `tf.random.set_seed`、并把 `-SD` 注入 `tfp.mcmc.sample_chain`。开启后重复运行可得到**逐字节一致**的结果，且几乎不增加耗时 |
| `-NT/--num_threads` | 0 | CPU 线程数；0 = TF 自行决定。**必须在 TF 导入前**通过环境变量生效（脚本已自动处理）。小样本 × 少细胞类型属于**延迟受限**（耗时主要来自 Python / TF 调用开销），用较小的线程数或默认值通常就够，加大线程数收益有限 |
| `--fdr` | 0.05 | 目标错误发现率（spike-and-slab spike 概率的估计 FDR）|
| `--fdr_sweep` | 0.05,0.1,0.2,0.4 | 同时给出多 FDR 的可信效应表，便于敏感性分析 |
| `--no-plot` | 关 | 跳过全部图（组成桑基、离散度、箱线、效应森林），只出数 + 采样 |
| `-RS/--rscript` | 自动查找 | `Rscript` 路径（组成桑基图用）。查找顺序：本参数 → `PATH` → 与本 Python 解释器同目录 |
| `-ET/--end_celltype` | 无 | 桑基图中**固定排在最后**的细胞类型（可空）。**排序规则本身是全局丰度降序**（对齐参考实现）。注：参考脚本里此参数叫 `-EC`，本脚本 `-EC` 已被 `--extra_covariates` 占用，故改名 `-ET` |
| `-P/--prefix` | `<输入前缀>.` | 输出文件名前缀 |
| `-O/--outpath` | 必填 | 输出目录（默认拒绝与输入同目录）|
| `-W/--width` `-H/--height` | 8 / 6 | 图尺寸（英寸）|
| `-dpi` | 300 | 图分辨率 |
| `--verbose_sampler` | 关 | 打开 HMC 步级进度；25,000 步会刷出几万行日志，默认关 |

### 输出

| 文件 | 说明 |
|---|---|
| `<前缀>composition_sankey.png/.pdf` | 各组细胞组成桑基（冲积）图。**对照组固定显示在上、处理组在下**；细胞类型按**全局丰度降序**排列（可用 `-ET` 把某类型钉到末尾）。由**必需配套脚本** `sankey_composition.r`（须与本 py 同目录）用 R + `ggalluvial` 渲染（需 Rscript 可用；可用 `-RS` 指定路径）|
| `<前缀>composition_boxplot.png/.pdf` | 各细胞类型相对丰度箱线图（按组分面）|
| `<前缀>composition_dispersion.png/.pdf` | 各细胞类型在样本间的相对丰度离散度——参考类型选择依据 |
| `<前缀>composition_effects.png/.pdf` | 效应量森林图：log2FC + 95% HDI，红色 = 在 `--fdr` 下可信 |
| `<前缀>composition_intercepts.csv` | 截距表 |
| `<前缀>composition_effects.csv` | 效应表（log2FC、HDI、最终纳入 spike 的概率）|
| `<前缀>composition_credible_effects.csv` | 在 `--fdr` 下可信变化的细胞类型（布尔列）|
| `<前缀>composition_fdr_sweep.csv` | 多 FDR 下的可信效应对比 |
| `<前缀>composition_result.pkl` | 完整采样结果（ARVI InferenceData + 模型），用 `pickle.load(open(p,'rb'))` 读回 |
| `<前缀>Composition_sccoda_summary.txt` | 全过程汇总（配置、FDR sweep、耗时）|

示例输出：`Images/composition_sccoda/`。

### 注意事项

- **输入文件三种形式**：`h5ad`（仅读 `adata.obs` 用 `backed='r'`，不加载表达矩阵，1 GB 文件也快）/ `counts`（每行一样本）/ `obs`（每行一细胞）。**R 的 .rds 不能直接被 Python 读**（S4 对象），需先转换为 h5ad（例如 `sceasy::convertFormat`、`zellkonverter::SceToAnnData`、`SeuratDisk::SaveH5Seurat` + `Convert`）。
- **每个样本只对应一个分组**：`sample × group` 必须一一对应，否则脚本会显式报错。这是有意约束——scCODA 以"样本"为最小单位。
- **零计数自动加伪计数 0.5**：脚本会打印 `Zero counts encountered in data! Added a pseudocount of 0.5.`，这是 scCODA 的内置 Dirichlet-multinomial 处理。**报告中需说明**：这些 0 是"该样本未检出该细胞类型"（不是"表达为零"）。
- **参考细胞类型 (`-R`)**：scCODA 的数学结构要求选定一个细胞类型作为"参考"——其效应被自动设为 0、无法被检出。**这是成分数据 N-1 自由度的代数约束**，详见下方"参考细胞类型选择"小节。默认 `automatic`，对大样本（≥ 8 例）够用；小样本下需显式指定。
- **采样耗时与问题规模强相关**：MCMC 迭代时间与问题规模近似线性关系。**小样本 + 少 cell type** 属于**延迟受限**：主要时间花在 Python/TF 调用开销上，并行收益有限，此时**增大 `-N`**（如 50,000 / 100,000）以提高后验估计的稳定性是更好的选择。对于大样本（≥30），HMC 是主要计算成本，线程数才有意义。
- **spike-and-slab 的 `diverging` 警告**：README 明确说明 spike-and-slab 先验会导致部分样本被标记为 `diverging`——这是先验特性，**不当作采样失败**。
- **报告数字必须实测**：可信效应表只显示 True/False；具体 log2FC 与 HDI 看 `composition_effects.csv`，并与 `composition_dispersion.png`（各 celltype 在样本间的相对丰度）结合——单纯看 credible 表容易把"微小但可信"的变化与"大但有 HDI 跨越 0"的变化混在一起。
- **"结果可复现"需要显式开启 `--deterministic`**：只设 `-SD` 是**不够的**——相同参数重复运行，各细胞类型的 Final Parameter 与 inclusion probability 仍会漂移。三个根因：① scCODA 构造 `CompositionalAnalysis` 时用 `tf.random.normal(...)` 抽 MCMC 初始值，**若种子设得太晚（晚于模型构造）这一步就是随机的**；② scCODA 调用 `tfp.mcmc.sample_chain` 时**没有传 seed**；③ TF 的 CPU 浮点算子默认不保证逐位确定。`--deterministic` 同时修掉这三点，开启后重复运行结果**逐字节一致**。
- **临界细胞仍会因采样量不足而摆动**：即使逐位复现，Monte Carlo 误差也会让贴近判定阈值的细胞类型（inclusion probability 长期在阈值附近）在采样量不足时结果不稳。**处理办法：扩大 `-N`** 让它收敛后再判定；对收敛后仍贴近阈值的细胞应标注"临界/不稳健"。
- **不要用超大线程数做最终结果**：本脚本默认 `-NT 0`（交给 TF 决定）。超大线程数收益极小（见上），最终交付建议用默认值或 `-NT 1`（配合 `--deterministic`）

### 参考细胞类型选择（`-R` 参数方法论）

scCODA 在数学上要求选定一个"参考"细胞类型——其效应强制设为 0。这是成分数据 N-1 自由度的代数约束：每行（样本）的细胞数总和固定，因此 n 个细胞类型中只有 n-1 个独立自由度。选定参考后，其余细胞类型的相对丰度变化可以无歧义估计，但**参考类型本身的变化在该分析中无法回答**（被吸收到 intercept）。

#### 选参考类型的标准

| 标准 | 说明 |
|---|---|
| **1. 必须：参考类型不能随分组变化** | 选完参考 = 假设它不变。若参考实际上在变，所有其他 cell type 的 log-ratio 估计会被它的变化污染。 |
| **2. 推荐：生物学上"稳态"的类型** | 通用淋巴细胞（如 T、B）通常较稳；髓系细胞（如中性粒细胞、巨噬细胞）易随炎症变化；上皮细胞易随结构损伤变化。具体选哪个取决于研究背景。 |
| **3. 推荐：丰度高、跨样本稳定出现** | 避免选有零值的——scCODA 内部会对 0 加伪计数 0.5，污染 log-ratio。 |
| **4. 避免：本身就疑似在变的** | 计算每组的相对丰度，若两组差 2× 以上则该类型不适合作参考。 |
| **5. 实务：跑多次取交集** | 换 2–3 个候选参考各跑一次，**所有参考都 True 的结论最稳健**。这是 scCODA 论文推荐的方法。 |

#### `automatic` 模式怎么选、什么时候够用

`automatic` 的选取逻辑：**所有出现率 ≥ 95% 样本的 cell type 中，相对丰度离散度最低的那一个**。

| 场景 | `automatic` 是否适用 |
|---|---|
| 样本数 ≥ 8，且没有"明显在变"的 cell type | ✓ 可信 |
| 样本数少（"95% 出现"这一门槛在小样本下会退化为更低的有效样本数）| ⚠️ 可能选到低丰度、低统计把握度的 cell type |
| 你有先验知识知道哪个 cell type 不变 | ✓ 仍可显式指定；**默认就是 automatic**——只是显式覆盖更稳 |
| 想让模型自己"无偏"地选 | ✓ 用 automatic 即可 |

`automatic` 的**常见风险**：当满足"出现率 ≥ 95%" 门槛的 cell type 很少时，它倾向于选**丰度极低但 absolute 波动也极小**的 cell type。该 cell type 作参考时统计把握度差，最终 log-ratio 估计对参考的具体选择会变得敏感。

#### 推荐的工作流

```
Step 1: 不传 -R，用默认 automatic 跑一次 → 记录选了谁
Step 2: 显式传 -R <celltype>（选丰度高 + 跨样本稳 + 生物学上预期不变的类型）
Step 3: 对比两次结果
        交集（两次都 True） = 最稳健结论
        并集               = 可能更宽松
Step 4: 报告里说明：用什么参考、为什么、两次一致性如何
```

#### 复用建议

- 样本数 < 10 时**必须显式指定** `-R`，不要用默认
- 样本数 ≥ 30 时可以直接用 `automatic`
- 不论哪种情况，**至少跑 2 个参考取交集**作为最终报告

### 已知坑

- **桑基图（Plot 1）需要 R + `ggalluvial`**：组成桑基图由**必需配套脚本** `sankey_composition.r`（须与本 py 同目录）用 R 的 `ggalluvial` 渲染（不重新实现，以保证与原参考实现一致——纯 Python 复刻的流带轮廓明显不如原包）。运行环境需有 `Rscript` 与 `ggalluvial` 包；脚本按"`-RS` 显式路径 → PATH → 与本 Python 解释器同目录"的顺序自动查找。**若 R 或 `ggalluvial` 缺失，该图会跳过并在日志中说明（`Plot 1 skipped (...)`），其余产物不受影响。**
- **h5ad 读取需 anndata 与 zarr 兼容**：`anndata < 0.12` 硬拒 zarr ≥3；`anndata ≥ 0.13` 要求 numpy ≥2.1。旧版环境中常见故障。
- **效应图文件名固定为** `<prefix>composition_effects.png/.pdf`（不再把协变量水平拼进文件名；多协变量水平时后写的会覆盖前一个，如需区分请分别指定 `-P`）。
- **导入 TensorFlow 较慢**：scCODA 依赖 tensorflow / tensorflow-probability，首次导入约 5–10 s（之后命中缓存）。

---

## 9 scRNA_spRNA_h5ad_preprocess.py（独立脚本）—— h5ad 计数矩阵预处理

### 功能

把 h5ad 整理成后续分析可直接使用的形态，做三件事：

1. **确定哪一份矩阵是原始整数 counts**：先看 `adata.X`（随机抽 1000 个细胞，按 1e-6 容差判断取值是否全为整数）；不是整数则回退 `adata.raw.X`；两者都不是整数则报错中止。
2. **重设基因名并合并同名基因**：用 `-V` 指定的 `adata.var` 列作基因名，同名基因按细胞**求和**合并（sparse COO + `sum_duplicates`，整数提升为 int64 防溢出，不落稠密矩阵）。
3. **产出标准结构**：把计数同时写入 `X` 与 `adata.layers["counts"]`，供下游（尤其 `scRNA_h5ad2rds.r`）直接取用。

它是 `scRNA_h5ad2rds.r` 的前置：后者要求 h5ad 内已有 `layers["counts"]`。

### 输入

`-I` 目录（处理其中所有 `.h5ad`，扩展名大小写不敏感）或 `-F` 单个 h5ad 文件，二者互斥且必须给其一。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_dir` | 无 | h5ad 所在目录（与 `-F` 互斥）|
| `-F, --file` | 无 | 单个 h5ad 文件（与 `-I` 互斥）|
| `-V, --var_index` | 无 | `adata.var` 中用作基因名的列；不传则沿用原 index 并跳过同名基因合并（仅警告）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 单文件
./scRNA_spRNA_h5ad_preprocess.py -F input.h5ad -V gene_symbols -P 00. -O 00_preprocess

# 整个目录
./scRNA_spRNA_h5ad_preprocess.py -I h5ad_dir -V gene_symbols -P 00. -O 00_preprocess
```

### 输出

`<前缀><原名>.h5ad`（gzip 压缩）：含 `X`、`layers["counts"]`、原 `obs`；var 只保留基因名 index。

终端打印：每个文件的 counts 判定结果（`[counts check]` 行的 is_integer / zero_fraction / min / max）、选用了哪个矩阵、去重前后基因数。

### 注意事项

- **先跑本脚本再转 rds**：`scRNA_h5ad2rds.r` 会删掉名为 `X` 的 layer，因此输入 h5ad 必须已有 `layers["counts"]`——本脚本正是产出该层的那一步。
- 计数判定是**抽样**（1000 细胞）且未固定随机种子：混合型数据（部分整数、部分非整数）可能判定不稳，看终端 `[counts check]` 行确认。
- 不传 `-V` 时保留原 var index 且**不合并同名基因**；若原 index 存在重复基因名，后续建 Seurat 对象会出问题。
- 同名基因合并是**求和**（不是取平均），与常见 10x 处理一致。
- `-I` 模式下文件顺序取自 `os.listdir()`（未排序），不影响结果。

---

## 10 scRNA_h5ad2rds.r（独立脚本）—— 批量 h5ad → Seurat rds

### 功能

把一批（或单个）h5ad 转成 Seurat `.rds`，两条后端：

- `-M anndataR`（默认，纯 R）：`anndataR::read_h5ad(path, as = "Seurat")`，随后删掉名为 `X` 的 layer，再 `saveRDS`。
- `-M sceasy`：`sceasy::convertFormat(from = "anndata", to = "seurat")`，经 reticulate 调用 `-P_env` 指定的 Python 环境；脚本会检查该环境已安装 `numpy` 与 `anndata`。

### 输入

`-I` 目录（取其中 `.h5ad`，大小写不敏感）或 `-F` 单文件，互斥必填其一。**h5ad 需含 `layers["counts"]`**（见注意事项）。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无 | h5ad 目录（与 `-F` 互斥）|
| `-F, --file` | 无 | 单个 h5ad（与 `-I` 互斥）|
| `-M, --method` | anndataR | 转换后端：`anndataR` 或 `sceasy` |
| `-T, --threads` | 12 | 线程数（仅 sceasy 分支注册 BiocParallel；见注意事项）|
| `-P_env` | 无 | `-M sceasy` 时的 Python 解释器路径（该分支必填）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 目录批量（默认 anndataR）
./scRNA_h5ad2rds.r -I 00_preprocess -O 00_rds -P 00.

# 单文件
./scRNA_h5ad2rds.r -F 00_preprocess/sample1.h5ad -O 00_rds -P 00.

# sceasy 后端（需指定 Python 环境）
./scRNA_h5ad2rds.r -I 00_preprocess -M sceasy -P_env /path/to/env/bin/python -O 00_rds
```

### 输出

`<前缀><原名>.rds`（`.h5ad` 后缀替换为 `.rds`），逐文件打印保存结果。

### 注意事项

- **输入必须带 `layers["counts"]`**：只有 `X` 的 h5ad 会在删除 `X` 层时触发 `invalid class "Assay5" object: 'default' must be between 0 and the number of layers present` 并中止。先跑 `scRNA_spRNA_h5ad_preprocess.py` 即可满足。
- anndataR 分支直接操作 `@layers`（S4 内部结构），升级 anndataR 版本后需重测。
- `-T` 目前无实际加速：sceasy 分支虽注册了 `MulticoreParam` / `SnowParam`，但循环内逐个 `convertFormat`；anndataR 分支完全不用该参数。
- sceasy 路径的 Python 环境需预装 `numpy` + `anndata`；脚本只检查、不安装。
- `-I` 模式下文件顺序取自 `list.files()`（未排序），不影响结果。

---

## 11 scRNA_h5ad_seurat_convert.r（独立脚本）—— 单文件双向格式转换

### 功能

按输入扩展名自动判定方向，单文件互转：

- 输入 `.h5ad` → `anndataR::read_h5ad(as = "Seurat")` → 存 `.rds`
- 输入 `.rds` → `readRDS` → 把 `meta.data` 的 **factor 列转 character**（否则写 h5ad 会报错）→ `as_AnnData(assay_name, x_mapping)` → `write_h5ad(compression)`
- 其他扩展名 → `Unsupported file format` 报错

### 输入

`-I` 单个 `.rds` 或 `.h5ad`（必填）。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | 输入文件（`.rds` 或 `.h5ad`）|
| `-O, --outpath` | 当前目录 | 输出目录 |
| `-P, --prefix` | 输入名（去扩展名） | 输出文件名前缀 |
| `-A, --assay_name` | RNA | 仅 rds → h5ad：要导出的 assay |
| `-X, --x_mapping` | counts | 仅 rds → h5ad：映射到 AnnData `X` 的层 |
| `-C, --compression` | gzip | 仅 rds → h5ad：`none` / `gzip` / `lzf` |

### 用法示例

```bash
# rds → h5ad（导出给 scanpy / cellxgene 等）
./scRNA_h5ad_seurat_convert.r -I 05_annotation/05.seu_annotated.rds -O 09_export -P export.

# h5ad → rds
./scRNA_h5ad_seurat_convert.r -I external.h5ad -O 00_rds -P 00. -A RNA -X counts
```

### 输出

`<前缀><输入名>.rds` 或 `<前缀><输入名>.h5ad`；结束打印对象信息与总耗时。

### 注意事项

- **h5ad → rds 不重命名层**：只有 `X` 的 h5ad 转出来的 Seurat 对象，层名就是 `X`（不是 `counts`）。本库 step1 用 `GetAssayData(layer = "counts")` 取计数，这种对象会找不到层——接主流程前请确认源 h5ad 带 `layers["counts"]`（或先跑 `scRNA_spRNA_h5ad_preprocess.py`）。
- rds → h5ad 会把 `meta.data` 的 factor 转成 character，导出的 h5ad 层名来自 `-X`（默认 `counts`）。
- 与 `scRNA_h5ad2rds.r -F` 功能重叠，但本脚本**双向**且不删层，行为更保守。
- `as_AnnData()` / `write_h5ad()` 未加包名前缀，依赖 `library(anndataR)` 已附加。


---

## 12 scRNA_mtx2rds.r（独立脚本）—— Broad SCP 矩阵转 Seurat 对象

### 功能

把 **Broad Single Cell Portal**（singlecell.broadinstitute.org）下载的表达矩阵和元数据直接组装成 Seurat 对象，省去手工读矩阵、对齐细胞名、构造对象的过程。两种输入模式（互斥，必须给其中一种）：

- **稠密矩阵**（`-E`）：SCP 导出的 `*_exprMatrix.tsv` / `.tsv.gz` / `.csv`，第一列是基因、表头是细胞名。用 `fread` 读，失败回退 `read.csv`。
- **稀疏三元组**（`-M` + `-C` + `-G`）：`.mtx` 矩阵 + 每行一个细胞名的文件 + 每行一个基因名的文件，脚本会校验行列数是否与矩阵一致。

元数据 `-m` 必填（TSV / CSV 自动识别分隔符），细胞名列由 `-CC` 指定，**默认 `NAME`**——正是 SCP metadata 的列名；若元数据首行首格含 `TYPE`，该行会被自动丢弃（SCP 导出常见）。

处理流程：读矩阵 → 基因名 `_` 统一换成 `-` 并 `make.unique()` 去重 → 读元数据并以其细胞名列作行名 → 取「矩阵列 ∩ 元数据行」的细胞 → NA 填 `"Unknown"` → 判定数值类型并建对象。

### 输入

- 表达矩阵：`-E` 单个稠密矩阵文件，或 `-M`/`-C`/`-G` 三个文件一起给
- 元数据：`-m`（必填）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-E, --expr_file` | 无 | 稠密表达矩阵（首列基因、表头细胞名）|
| `-M, --mtx_file` | 无 | 稀疏矩阵文件（`.mtx`）|
| `-C, --cell_file` | 无 | 每行一个细胞名，行数须等于矩阵列数 |
| `-G, --gene_file` | 无 | 每行一个基因名，行数须等于矩阵行数 |
| `-m, --meta_file` | （必填） | 元数据文件（TSV / CSV 自动识别）|
| `-CC, --cell_column` | NAME | 元数据中细胞名所在列 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# SCP 稠密矩阵 + 元数据
./scRNA_mtx2rds.r -E exprMatrix.tsv.gz -m metadata.tsv -CC NAME -P 00. -O 00_rds

# 稀疏三元组
./scRNA_mtx2rds.r -M matrix.mtx -C barcodes.tsv -G genes.tsv -m metadata.tsv -P 00. -O 00_rds
```

### 输出

`<前缀>seurat.rds`。终端打印矩阵维度、元数据维度、数值类型判定结果，以及对象概览（维度、metadata 列、前几个 feature 名）。

### 注意事项

- **数值类型决定层结构**：全整数 → `counts` + `data`（`data = log1p(counts)`）；出现非整数（SCP 有时直接给归一化矩阵）→ 造一个全 0 的 counts 再**删除 counts 层**，只保留 `data`。后者不能接主流程 step1 / step2（它们要读 `counts`），只能作为已归一化对象做降维/可视化。
- **`-E` 模式不接受重复基因名**：设 `rownames` 时会直接报 `duplicate 'row.names' are not allowed`（`make.unique()` 在其后执行，救不回来）；`-M` 模式能容忍重复，重复基因会被改名为 `G-A` / `G-A.1` 两个独立 feature，**不会合并求和**。
- **NA 处理是整表替换**：`meta[is.na(meta)] <- "Unknown"`，某个数值列只要有 NA，整列会被转成 character。
- **`-CC` 列名写错时报错不直观**：会落到 `CreateSeuratObject` 的 `Rownames must be supplied`。
- **细胞取交集**：矩阵与元数据细胞不完全重合时，按交集保留，多余细胞静默丢弃。
- 基因名统一把 `_` 换成 `-`（`gsub("_", "-")`），与常见 10x 基因名写法一致。
- 非整数分支里 `seo@assays$RNA@layers$counts <- NULL` 是直接操作 slot，仅适用 Seurat v5。


---

## 附录：环境与依赖

**R 侧**（主流程各步与独立 R 脚本）：

- 核心：`Seurat`（v5）、`tidyverse`（dplyr / ggplot2 / tidyr 等）、`patchwork`
- 参数解析：`argparse`
- 去双细胞：`scDblFinder`、`BiocParallel`；可选 DecontX：`celda`
- 批次整合与探索：`harmony`、`clustree`
- Leiden 聚类（可选）：`reticulate`（脚本内部 `py_require("leidenalg")`）
- `scRNA_gene_DE.r` 额外：`reshape2`
- 桑基图（`scRNA_composition_sccoda.py` 配套 `sankey_composition.r`）：`ggalluvial`、`RColorBrewer`、`scales`、`ggplot2`、`dplyr`、`tidyr`
- 格式转换（`scRNA_h5ad2rds.r`、`scRNA_h5ad_seurat_convert.r`）：`anndataR`（默认路径）；可选 `sceasy`（经 reticulate，其配套 Python 环境需预装 `numpy` + `anndata`）
- 矩阵读入（`scRNA_mtx2rds.r`）：`Matrix`（`readMM`）、`data.table`（`fread`）

**Python 侧**（`scRNA_spRNA_h5ad_preprocess.py`、`scRNA_composition_sccoda.py`）：

- `scCODA`、`tensorflow` / `tensorflow-probability`
- `anndata`（注意与 `zarr` / `numpy` 的版本匹配：`anndata < 0.12` 拒 zarr ≥3，`anndata ≥ 0.13` 要求 numpy ≥2.1）
- h5ad 预处理：`scanpy`、`anndata`、`scipy`、`numpy`、`pandas`
- 组成分析：`numpy`、`pandas`、`matplotlib`、`seaborn`
- 桑基图需 `Rscript` 与 R 的 `ggalluvial`（缺失时该图跳过，其余产物不受影响）

**解释器**：脚本 shebang 用 `/usr/bin/env` 解析 `Rscript` / `python3`，运行前请确保所用环境在 `PATH` 上（例如 `export PATH=/path/to/env/bin:$PATH`）；若不在 PATH 上，直接 `./<script>` 会报解释器未找到，此时用 `Rscript <script>.r …` / `/path/to/env/bin/python <script>.py …` 显式调用即可。

R 脚本启动时会检查缺失的依赖并提示安装命令；Python 侧缺包会在导入处报错。
