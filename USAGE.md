# scRNA 通用分析脚本库 —— 使用说明

本目录是一套**持续更新**的通用单细胞转录组（scRNA-seq）分析脚本库。其中**主分析流程脚本按步骤编号**（脚本名带 `stepN` 前缀），其余为**独立脚本**，不在流程编号内、输入各自独立，可单独使用。文档每个脚本一章，覆盖从原始计数矩阵到细胞注释、表达可视化与细胞组成差异分析的完整链路；所有脚本都以命令行参数接收输入（10x 目录、`.h5`、`.rds`、`.h5ad`、表格等路径），并把结果写成文件（数据对象、表格、图片等），因此每一步既可单独运行，也可串成完整流水线；脚本会持续补充。各脚本的示例可视化输出集中存放在本目录的 `Images/` 中（见下文“可视化输出示例”）。

脚本自带 shebang（R 脚本 `#!/usr/bin/env Rscript`，Python 脚本 `#!/usr/bin/env python3`）与自解释参数（`-h` 查看帮助）。只要 PATH 上有对应解释器，即可 `./<script>` 直接执行；也可用 `Rscript <script>.r …` / `python3 <script>.py …` 显式指定解释器。环境与依赖清单见文末附录。

## 流程概览

**主流程**（脚本名带 `stepN`，按顺序串联；前一步的 rds 即后一步的输入）：

| 步骤 | 脚本 | 用途 | 输入 → 输出 |
|---|---|---|---|
| step1 | `scRNA_QC_step1.r` | **质控并合并样本**：每个样本按自身分位数定阈值剔除低质量细胞，再把所有样本合并为一个对象 | 10x 目录 / h5 列表 → 合并 rds |
| step2 | `scRNA_dedoublets_step2.r` | **去除双细胞与环境污染**：逐样本建模识别双细胞（可选 DecontX 去除环境 RNA）并剔除 | 合并 rds → 干净 rds |
| step3 | `scRNA_Nor_harmony_step3.r` | **归一化并消除批次效应**：LogNormalize 后用 Harmony 把不同样本/批次的细胞整合到一起 | 干净 rds → 整合 rds |
| step4 | `scRNA_clustering_step4.r` | **聚类并找出各群差异基因**：多分辨率探索聚类结构，再在选定分辨率下聚类并做 marker 分析 | 整合 rds → 聚类 rds + markers |
| step5 | `scRNA_annotate_step5.r` | **给细胞群标注细胞类型**：用 marker 列表把每个 cluster 映射到具体细胞类型 | 聚类 rds + marker 表 → 注释 rds |

**独立脚本**（不参与主流程编号，可单独使用）：

| 脚本 | 用途 | 输入 → 输出 |
|---|---|---|
| `scRNA_spRNA_h5ad_preprocess.py` | **清洗 h5ad 计数矩阵**：判定并取出原始整数 counts（`X`，否则回退 `raw.X`）、按 `-V` 重设基因名并合并同名基因，得到可直接用于差异分析的结构 | h5ad（目录 / 单文件）→ 含 `X` + `layers["counts"]` 的 h5ad |
| `scRNA_h5ad2rds.r` | **批量 h5ad 转 Seurat rds**：两条后端可选，`anndataR`（纯 R）/ `sceasy`（Python）| h5ad（目录 / 单文件）→ rds |
| `scRNA_h5ad_seurat_convert.r` | **单文件 h5ad ⇄ rds 双向转换**：按扩展名自动判定转换方向 | h5ad ⇄ rds |
| `scRNA_mtx2rds.r` | **把 Broad Single Cell Portal 下载的数据组装成 Seurat 对象**：读入表达矩阵与元数据（支持稠密矩阵与稀疏三元组两种输入），对齐细胞名后建对象 | 矩阵文件（+ metadata）→ `seurat.rds` |
| `scRNA_subcelltype_clustering.r` | **在细胞类型内部细分亚型**：逐细胞类型重跑降维与聚类出 clustree，按指定分辨率定亚型并做亚型 DE，最后把亚型标签回写到原对象 | 注释 rds + 分辨率表 → 亚型 rds + 各亚型 DE |
| `scRNA_pseudobulk_DE.r` | **在每个细胞类型中，找出两组之间差异表达的基因**（每个细胞类型单独检验，样本层面的 DESeq2），并画火山图 | 注释 rds → DE csv + 火山图 |
| `scRNA_gene_expr.r` | **展示目标基因在不同细胞类型中的表达**：逐基因给出跨细胞类型的平均表达量，并生成 5 类表达分布图 | 已注释 rds → 均值 CSV + 5 类图 |
| `scRNA_gene_DE.r` | **比较单个目标基因在样本与细胞类型间的表达差异**：给出样本 × 细胞类型均值表与热图 / 柱状图 / 箱线图，两组间比较以样本为单位 | 已注释 rds → 表 + 图 |
| `scRNA_composition_sccoda.py` | **检验两组间细胞类型比例构成的差异**：用 Dirichlet-Multinomial 贝叶斯模型做正式成分差异分析，输出带可信区间的效应表 | h5ad / 计数表 → 效应表 + 图 |
| `scRNA_spRNA_cellPro.r` | **每个样本里各种细胞类型各占多少**：出样本级/组级占比表、组成柱状图、桑基图，并比较两组之间各细胞类型占比是否不同 | 注释 rds → 占比表 + 柱状图 + 桑基图 + 占比比较图 |
| `scRNA_SingleR.r` | **用参考数据集自动推断细胞类型**：把查询对象的每个聚类与参考数据比对，输出 cluster→cell type 映射并写回对象 | 已聚类 rds + 参考数据 → 注释结果 csv + rds + 图 |
| `scRNA_spRNA_enrichment.r` | **把差异基因翻译成通路层面的结论**：逐细胞类型做 GO/KEGG 富集与 GSEA，把「哪些基因变了」变成「哪些通路变了」 | 注释 rds + DE 结果 csv → 富集表 + 图 |
| `scRNA_GSVA_enrichment.r` | **比较两组之间的通路活性差异**：不依赖差异分析结果，直接用 GSVA 给每个细胞类型算通路活性分数并做组间检验 | 注释 rds → 通路分数表 + 差异表 + 热图 |
| `scRNA_genesetScore.r` | **给自定义基因集打分**：算 GSVA / ssGSEA / AddModuleScore / AUCell 四种分数，看这个基因集在哪些细胞类型或细胞里更活跃 | 注释 rds + 基因集文件 → 分数表 + 4 类图 |
| `scRNA_DoubleFinder.r` | **替代方案**：同样去双细胞与环境污染，但改用 DoubletFinder 且整体（不逐样本）处理；用于与 step2 交叉验证 | 合并 rds → 干净 rds |
| `scRNA_ROC.r` | **评估目标基因的组间判别能力**：逐基因在样本层面算 AUC 并画 ROC 曲线，给出按 AUC 排序的汇总表 | 注释 rds + 分组 csv + 基因列表 → AUC 表 + ROC 曲线 |
| `scRNA_cellchat_single.r` | **推断细胞类型之间的通讯信号**：用 CellChat 找出谁给谁发什么信号、走哪条通路，并给出网络图与通讯模式 | 注释 rds → 通讯表 + 网络图 + 模式图 |
| `scRNA_cellchat_compare.r` | **比较两组的通讯网络差异**：合并两个 CellChat 对象，逐项对比通讯数量、强度与通路组成 | 两个 CellChat rds → 合并对象 + 对比图 |
| `scRNA_spRNA_Monocle3.r` | **给细胞排出发育/转变的先后顺序**：用 Monocle3 学习轨迹图、指定根细胞、计算拟时序，并找出沿轨迹变化的基因与基因模块 | 注释 rds（或空间对象）→ 拟时序图 + Moran's I 表 + 模块表 + 热图 |
| `scRNA_velocity.py` | **推断细胞的状态变化方向**：用 scVelo 算 RNA velocity，叠在 UMAP 上得到流线图，并给出 latent time 与驱动转变的基因 | velocity 输入 h5ad → 流线图 + latent time 图 + 结果 h5ad |
| `scRNA_velocityPrepare.r` | **配套脚本**：把含 spliced/unspliced 的完整对象与分析子集对齐，生成带 RNA / spliced / unspliced 三层的 velocity 输入对象 | 参考 rds + 完整 rds → velocity 输入 rds + h5ad |
| `scRNA_pyscenic_prepare.r` | **准备 pySCENIC 输入**：从 Seurat 对象过滤低检出基因、可选抽细胞，写成带 Gene/CellID 属性的 loom | Seurat rds → loom |
| `scRNA_pyscenic.py` | **推断转录因子调控网络**：依次调 pyscenic 的 grn（共表达）、ctx（motif 富集得到 regulon）、aucell（逐细胞活性打分） | loom + 参考文件 → adj.csv + regulons.csv + SCENIC.loom |
| `scRNA_pyscenic_RSS.r` | **找出每类细胞的标志转录因子**：算 RSS 特异性分数、出 5 张图，并把调控活性写回 Seurat 对象（5 张图的读法见第 28 章「结果图怎么读」）| Seurat rds + SCENIC loom → rss.csv + 图 + seu_SCENIC.rds |
| `scRNA_scTenifoldKnk.r` | **预测「敲掉某个基因」会连带影响哪些基因**：用单细胞数据重建调控网络，在网络里把目标基因的连接置零做虚拟敲除，比较敲除前后的网络差异，得到受影响基因表与富集通路 | Seurat rds（需 counts 层 + HVG）+ 基因名 → diffRegulation 表 + 4 张图 |
| `scRNA_hdWGCNA.r` | **找出某个细胞类型内部的基因共表达模块与核心 hub 基因**：在 metacell 上建 WGCNA 网络，划分模块、算 kME、给细胞打模块活性分，并出模块诊断与可视化图 | Seurat rds（需 harmony 降维）→ 模块表 + hub 基因 + MEs/hMEs + 10 张图 + seu_hdWGCNA.rds |
| `enrichment_analysis.r` | **配套脚本**：执行实际的 GO/KEGG 富集与 GSEA（第 17 章脚本调用它）| 基因列表 / 排序表 → 富集表 + 图 |
| `sankey_composition.r` | **配套脚本**：用 `ggalluvial` 渲染细胞组成桑基图（第 8 章脚本调用它）| 组成占比表 → 桑基图 |

**Scanpy 路线（`scanpy/` 子目录，输入输出为 `.h5ad`）**

| 脚本 | 用途 | 输入 → 输出 |
|---|---|---|
| `scanpy/scRNA_pipeline_scanpy.py` | **质控、去双细胞、归一化、整合、聚类并找出各群差异基因**：一条命令走完前四步预处理，产出可注释的对象 | h5ad → preprocessed.h5ad |
| `scanpy/scRNA_annotate_scanpy.py` | **给细胞群标注细胞类型**：先按 marker 表画 dotplot 判断 cluster 身份，再按映射表写回 celltype | preprocessed.h5ad → annotated.h5ad |
| `scanpy/scRNA_DE_scanpy.py` | **找出每个细胞群的 marker 基因**：Wilcoxon / t-test / logreg 差异检验 + 概览图与逐群 dotplot | h5ad → DE 表 + 图 + DE_analyzed.h5ad |
| `scanpy/scRNA_pipeline_scanpy_GPU.py` | **GPU 加速版的前四步预处理**：流程与 CPU 版完全相同，改用 rapids-singlecell，需 CUDA 显卡 | h5ad → preprocessed.h5ad |

> Scanpy 路线与 Seurat 路线**产物格式不同**（`.h5ad` vs `.rds`），不能混用；跨路线请用 `scRNA_h5ad2rds.r` / `scRNA_h5ad_seurat_convert.r` 转换。详见「Scanpy 路线」章。

> `enrichment_analysis.r` 不是独立调用的分析脚本，而是 `scRNA_spRNA_enrichment.r` 的**必需配套脚本** —— 后者用 `system()` 调用它做实际的 GO/KEGG 富集与 GSEA（命令名 `enrichment_analysis` 需在 `PATH` 中；本文件提供的是 `.r` 版本，命名为 `enrichment_analysis` 的软链/副本即可被调用）。缺该命令时主脚本仍会跑完并生成所有中间文件，但富集结果为空。
>
> `sankey_composition.r` 不是独立调用的分析脚本，而是 `scRNA_composition_sccoda.py` 的**必需配套脚本** —— 后者用 R + `ggalluvial` 渲染细胞组成桑基图，由 py 自动调用（无需手动运行），因此该文件**必须与 `scRNA_composition_sccoda.py` 放在同一目录**。缺 R 或 `ggalluvial` 时仅桑基图跳过，其余产物不受影响。

使用约定：

- 主流程各步以文件衔接，中间产物可单独复用；建议每步输出到独立目录（如 `01_QC/`、`02_deContamination/`），输出前缀（`01.` `02.` …）与步骤号对应，便于追溯。
- 数据准备/格式转换类脚本（`scRNA_spRNA_h5ad_preprocess.py` → `scRNA_h5ad2rds.r`）是从 h5ad 起步时的前置链：先清洗计数矩阵，再转成 Seurat rds，之后即可接主流程 step1；`scRNA_h5ad_seurat_convert.r` 负责单文件双向互转（含 rds → h5ad 导出）。
- 从 Broad Single Cell Portal（singlecell.broadinstitute.org）下载的数据可用 `scRNA_mtx2rds.r` 直接建对象：稠密矩阵走 `-E`，稀疏三元组走 `-M/-C/-G`，元数据里细胞名列用 `-CC` 指定（SCP 导出默认列名就是 `NAME`）。
- 注释完成后若要对每个细胞类型再细分亚型，用 `scRNA_subcelltype_clustering.r`：它按 `-AC` 指定的细胞类型逐个重跑降维与聚类，先出多分辨率 clustree 供人工选分辨率，再按 `-Res` 指定的「细胞类型 → 分辨率」表定亚型、跑亚型 DE，最后把 `sub_celltype` 标签回写到原对象。
- 想回答「**每种细胞类型里，两组之间的差异基因是哪些**」用 `scRNA_pseudobulk_DE.r`：它按「细胞类型 × 样本」把原始 counts 求和成 pseudobulk，再用 DESeq2 在**样本层面**比较 `-G1` vs `-G2`。重复单位是样本而不是细胞，所以 p 值不会因细胞数多而虚高，结论可外推到新样本。
- 想看「每个样本里各细胞类型占多少、两组之间比例有没有差别」用 `scRNA_spRNA_cellPro.r`：它按样本统计各细胞类型的数量与占比，出计数表/占比表、100% 堆叠柱状图、组成桑基图，并在两组时逐细胞类型做 Wilcoxon 检验（**每个点是一个样本**）。若要的是带贝叶斯置信区间的严格成分差异检验，用 `scRNA_composition_sccoda.py`；本脚本是快速描述 + 简单检验。
- 想快速、无marker 地给细胞类型「起个名字」时用 `scRNA_SingleR.r`：它把查询对象按 `seurat_clusters` 与参考数据集比对（`-S human`/`-S mouse` 用 celldex 内置参考，或 `-RDS`/`-RD` 用自己的参考），输出 cluster→cell type 映射并写回对象。这是**自动推断**，结果需结合 marker 人工确认，不能替代 step5 的 marker 注释。
- DE 结果拿到之后，想知道「这些差异基因背后是哪条通路在变」用 `scRNA_spRNA_enrichment.r`：它逐细胞类型做 GO/KEGG 富集（用该类型自身检出的基因作背景）与 GSEA（用全基因按 log2FC 排序），输出各通路的 p/q 值、NES 与 dotplot/barplot。
- 想比较「两组之间哪些通路的活性不一样」但手上还没有 DE 结果时，用 `scRNA_GSVA_enrichment.r`：它按「样本 × 细胞类型」聚合成伪 bulk，算 log2(CPM+1) 后用 GSVA 打分（GO:BP / KEGG MEDICUS / REACTOME 三套基因集），再用 limma 在样本层面比较两组。它与 `scRNA_spRNA_enrichment.r` 互补：前者看「通路活性」，后者看「差异基因落在哪些通路上」。
- 想看「某个自定义基因集在哪些细胞类型里更活跃」用 `scRNA_genesetScore.r`：它同时给 GSVA / ssGSEA（细胞类型层面）与 AddModuleScore / AUCell（单细胞层面）四种分数，可互相印证，也能看出这个基因集是不是被某一个细胞类型特异性表达。
- 去双细胞有两个脚本，**二选一，不要串联**：`scRNA_dedoublets_step2.r`（主流程推荐，scDblFinder，逐样本建模、双细胞率可自动估计、输出只保留 counts） 与 `scRNA_DoubleFinder.r`（替代方案，DoubletFinder，整体处理、双细胞率需手填 `-DR`）。后者主要用于**交叉验证**——两者是独立算法，结论差异大时值得警惕。
- pySCENIC 是三脚本一条链，**必须按顺序执行**：`scRNA_pyscenic_prepare.r`（Seurat → loom）→ `scRNA_pyscenic.py`（grn/ctx/aucell）→ `scRNA_pyscenic_RSS.r`（RSS + 出图 + 活性写回对象）。三者跨 R 与 Python 两个环境，靠各自的 shebang 切换，不需要手动 activate。`SCENIC.loom` **只能用 R 的 `SCopeLoomR` 读**，Python 的 `loompy` 读不了。
- RNA velocity 分析分两步：先跑 `scRNA_velocityPrepare.r` 把「含 spliced/unspliced 的完整对象」与「分析用的子集对象」对齐，生成同时带 RNA / spliced / unspliced 三个 assay 的 h5ad；再跑 `scRNA_velocity.py` 做实际计算。**velocity 推断的是表达状态的变化方向，不等于真实谱系关系**，需要结合生物学解释。
- 拟时序分析（`scRNA_spRNA_Monocle3.r`）是**有方向**的：根细胞由 `-C` + `-G` 指定，选错根整条时间轴方向就错了——先用轨迹图看懂结构再定根。细胞群之间的**比例**变化它不建模，那要用 scCODA。
- 独立脚本不依赖主流程的步骤顺序：`scRNA_gene_expr.r`、`scRNA_gene_DE.r` 需要一个含 celltype 列的已注释 Seurat 对象（典型来源是 step5 的 `seu_annotated.rds`）；`scRNA_composition_sccoda.py` 直接读 h5ad 或样本 × 细胞类型计数表，与主流程无关。
- 脚本会识别已存在的关键产物作为缓存（通常只重画图、不重算），强制重跑时删除对应 rds 即可。

## 可视化输出示例（`Images/`）

本目录下的 `Images/` 存放各脚本的示例可视化输出，子目录与产出图的脚本一一对应（数据准备 / 格式转换类脚本不产图，因此没有对应子目录）。

**示例图只保留 PNG**：脚本运行时同一张图会同时输出 `.png` 与 `.pdf` 两个文件，但 `Images/` 中**只保留 PNG**，同名 PDF 已全部删除（既省体积、也便于网页直接浏览）。因此下表各脚本「输出」小节里写成 `.png/.pdf` 的文件，在 `Images/` 中只能找到 `.png`；需要 PDF 时用相应脚本重新跑一次即可。

| 子目录 | 对应脚本 | 内容 |
|---|---|---|
| `Images/01_QC/` | `scRNA_QC_step1.r` | QC 前后小提琴图 / 散点图 |
| `Images/02_deContamination/` | `scRNA_dedoublets_step2.r` | doublet UMAP、doublet/QC 指标小提琴图 |
| `Images/DoubleFinder/` | `scRNA_DoubleFinder.r` | 替代方案的产物：doublet UMAP、doublet 小提琴图、DecontX 污染 UMAP |
| `Images/03_Nor_Harmony/` | `scRNA_Nor_harmony_step3.r` | PCA dimplot、PC elbow、HVG、Harmony dimplot |
| `Images/04_clustering/`、`Images/04_clustering_louvain/` | `scRNA_clustering_step4.r` | clustree、cluster / sample / group UMAP、各 cluster dotplot（两个目录为同一分析在 Leiden / Louvain 下的示例） |
| `Images/05_annotation/` | `scRNA_annotate_step5.r` | celltype UMAP、celltype × marker dotplot、`-S` 模式下每 cell type 一张 dotplot |
| `Images/gene_expr/` | `scRNA_gene_expr.r` | Ridge / Violin / Feature / Dot / Heatmap |
| `Images/gene_DE/` | `scRNA_gene_DE.r` | 热图、柱状图、箱线图 |
| `Images/composition_sccoda/` | `scRNA_composition_sccoda.py` | 组成桑基、相对丰度箱线图、离散度、效应森林图 |
| `Images/subcelltype_clustering/` | `scRNA_subcelltype_clustering.r` | 每个细胞类型的 clustree、亚型 UMAP，以及各亚型 DE 目录下的 marker dotplot |
| `Images/cellPro/` | `scRNA_spRNA_cellPro.r` | 100% 堆叠柱状图（各样本的细胞类型构成）、组→细胞类型桑基图、逐细胞类型的两组占比比较图 |
| `Images/pseudobulk_DE/` | `scRNA_pseudobulk_DE.r` | 按细胞类型分区的差异表达火山图（底部色块为细胞类型，红上/蓝下，标注每型 top 差异基因）|
| `Images/SingleR/` | `scRNA_SingleR.r` | SingleR 得分热图（各聚类与各参考标签的相似度）、按推断标签着色的 UMAP |
| `Images/enrichment/` | `scRNA_spRNA_enrichment.r` | GO 富集 dotplot / barplot、KEGG 富集 dotplot / barplot、GSEA 的 GO 与 KEGG 图 |
| `Images/GSVA_enrichment/` | `scRNA_GSVA_enrichment.r` | 通路活性热图：行 = 通路（按 GOBP / KEGG / REACTOME 分区），列 = 样本（对照组在左），右侧标注数据库来源 |
| `Images/genesetScore/` | `scRNA_genesetScore.r` | GSVA 分数柱状图、ssGSEA 分数柱状图、AddModuleScore 与 AUCell 的 UMAP 特征图 |
| `Images/ROC/` | `scRNA_ROC.r` | 目标基因的 ROC 曲线图（示例为判别能力强、中等、弱三档）|
| `Images/cellchat/single/` | `scRNA_cellchat_single.r` | 单对象分析的**全部 16 类图**：总体网络圈图、逐发送方圈图、top9 信号通路弦图、selectK 参考图、通讯模式热图/河流图/点图、通路 LR 贡献与表达图、信号角色网络与散点、贡献最大的 LR 对圈图、全部 LR 对气泡图 |
| `Images/scanpy/pipeline/` | `scanpy/scRNA_pipeline_scanpy.py` | 质控前后小提琴图、高变基因、PCA 方差贡献、6 个分辨率的 UMAP 网格、各 cluster marker 图 |
| `Images/scanpy/annotate/` | `scanpy/scRNA_annotate_scanpy.py` | cluster 的 UMAP、按 celltype 着色的 UMAP、细胞类型 × marker 的 dotplot |
| `Images/scanpy/DE/` | `scanpy/scRNA_DE_scanpy.py` | marker 概览图、单个 cluster 的 top marker dotplot |
| `Images/scTenifoldKnk/` | `scRNA_scTenifoldKnk.r` | 虚拟敲除效应网络图（中心为目标基因、节点按富集通路着色）、top20 受影响基因条形图、Z vs −log10(p.adj) 散点图、显著基因占比饼图 |
| `Images/hdWGCNA/` | `scRNA_hdWGCNA.r` | 软阈值诊断图、基因聚类树与模块颜色条、各模块 kME 排序、hME / hub 基因打分投 UMAP、模块雷达图、模块相关性、细胞类型 × 模块 dotplot |
| `Images/pyscenic/` | `scRNA_pyscenic_RSS.r` | RSS 特异性总图（点大小 = RSS、颜色 = Z）、逐细胞类型 RSS 排名、regulon 活性小提琴图 / Feature 图 / 点图 |
| `Images/velocity/` | `scRNA_velocity.py` | RNA velocity 分析：dynamical 流线图、velocity genes 图、latent time 与 velocity pseudotime 的 UMAP 分布、stochastic 流线图 |
| `Images/Monocle3/` | `scRNA_spRNA_Monocle3.r` | 拟时序与轨迹分析：UMAP 轨迹图（按 Monocle3 簇着色，含 principal graph）、拟时序着色 UMAP、沿拟时序变化基因的表达趋势、基因模块热图、空间坐标上的拟时序分布 |
| `Images/cellchat/compare/` | `scRNA_cellchat_compare.r` | 两组对比的**全部 12 类图**：通讯数量/强度对比柱状图、差异网络圈图、2×2 总体圈图、强度对比热图、信号角色 2D 散点、通路功能相似性聚类与排名、组间差异显著通路图、outgoing 通路热图、LR 对比较气泡图（含只显示增强部分）、差异最大通路的 LR 表达对比 |

各章「输出」小节末尾给出对应的示例目录。示例图仅用于展示输出形式与排版，其文件名与各章输出的 `<前缀><文件名>` 规则一致。

---

## 1 scRNA_QC_step1 —— 每样本独立 QC + 合并

### 功能

**回答的问题：哪些细胞该被剔除？把多个样本清理干净并合并成一个可分析的对象。**

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

**回答的问题：哪些细胞是「两个细胞粘在一起」的双细胞？哪些是环境 RNA 污染？**

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

**回答的问题：不同样本/批次的细胞混在一起了吗？把批次效应消掉，让同一类细胞聚到一起。**

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

**回答的问题：这些细胞可以分成几群？每群之间哪些基因有差异？**

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

**回答的问题：每一群细胞是什么细胞类型？**

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
| `-CO, --celltype_order` | 无 | dotplot 中 cell type 顺序（逗号/空格分隔，如 `A, B, C`）|
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

**回答的问题：我想看的这几个基因，在不同细胞类型里分别表达多少、长什么样？**

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
| `-CTO, --celltype_order` | 按列默认排序 | celltype 在图中的固定顺序（逗号/空格分隔，如 `A,B,C`） |
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
- **celltype 顺序**：默认沿用对象内的水平顺序；`-CTO` 可指定固定顺序（如 `-CTO A,B,C`）

---

## 7 scRNA_gene_DE（独立脚本）—— 样本 × 细胞类型表达矩阵 + 热图 + 柱状图 + 箱线图

### 功能

**回答的问题：某个基因在每种细胞类型里、在各个样本中的表达量如何？两组样本之间这个基因有差别吗？**

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
| `--alt` | **two.sided** | 箱线图 Wilcoxon 检验方向：`two.sided`（默认，无方向）/ `auto`（等效 `less`，即右侧组 > 左侧组）/ `greater` / `less` |
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
| `<前缀><基因名>_boxplot_stats.csv` | **检验结果表** — 每个 celltype 一行：分组细胞数、中位数、均值、检验名、`p_value`、`p_adj`（对同一基因的各细胞类型做 BH 校正）|
| `<前缀>GeneDE_summary.txt` | 参数与统计汇总 |

示例输出：`Images/gene_DE/`。

### 注意事项

- **单基因输入**：`-G` 必须是单个基因符号；要批量分析请循环调用脚本
- **缺失组合填 0**：某 sample × celltype 无细胞时对应值为 0。注意这与"真实零表达"不可区分——箱线图中若某样本该 celltype 无细胞，其 0 会参与检验，解读时须结合每组细胞数
- **箱线图的点 = 样本，不是细胞**：用细胞做重复会产生伪重复（p 值被细胞数撑到 1e-16）。
- **默认是双侧检验**（`--alt two.sided`）：两侧 Wilcoxon 在两组完全分离时的最小 p 值为 `2/C(n1+n2, n1)` —— 3 vs 3 时 = 0.1、3 vs 4 时 = 0.0571、4 vs 4 时 = 0.0286、5 vs 5 时 = 0.0079。也就是说**每组只有 3 个样本时，双侧检验无论差异多大都不可能 < 0.05**，此时箱线图只能看差异方向和幅度，不能作显著性结论。单侧检验（`--alt less/greater`）的最小 p 是双侧的一半（3 vs 3 为 0.05），但方向必须在看数据之前就确定，否则属事后选择。
- **`p_adj` 是 BH 校正值，校正范围是同一次运行里各细胞类型的检验**：脚本一次只处理一个基因，所以校正族 = 该基因在 `-CT` 指定的各个细胞类型上的检验。跨基因不做校正——若要比较多个基因，请自行在合并后的表上再校正一次。无法检验的细胞类型（某组样本数 < 2）`p_value` 为 `NA`，`p_adj` 也是 `NA`，但仍计入校正的检验个数。
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

### 功能

**回答的问题：两组之间，细胞类型的「比例构成」有没有变化？（不是表达量的差异，而是「有多少这种细胞」）**

scCODA 是基于 Dirichlet-Multinomial 模型的贝叶斯成分差异分析方法。给定一个**样本 × 细胞类型**的计数表，识别哪些细胞类型在两组之间的**相对丰度**显著不同（同时避免"一类升高导致另一类机械性下降"的成分数据假阳性）。

### 输入

三种输入格式，脚本按扩展名与列名自动判定（也可用 `-T` 强制指定）：

| 格式 | 说明 | 需要的列 |
|---|---|---|
| `.h5ad` | AnnData 对象 | `obs` 里的样本列（`-S`）与细胞类型列（`-CT`）|
| 计数表 | 样本 × 细胞类型，行=样本、列=细胞类型 | 首列为样本名 |
| 每行一细胞 | 长表 | 样本列（`-S`）与细胞类型列（`-CT`）|

其中 h5ad 与长表模式还会读 `-GC` 指定的分组列。

### 用法示例

```
./scRNA_composition_sccoda.py -I <input.h5ad> \
    -S sample -CT celltype -GC group -R <reference> -GR <ref_level> \
    --deterministic \
    -O <output_dir> -P <prefix>.
```

`<reference>` 取 `automatic`（默认）或某个 cell type 名称；`<ref_level>` 是该分组列的参考水平。

**建议正式分析一律带上 `--deterministic`**——否则同样的命令重复运行会得到略有差异的结果（只有设了种子也不够，原因见下）；带上它可保证结果**逐位可复现**、报告引用的数字永远对得上，代价几乎为零。桑基图（Plot 1）由同目录的 `sankey_composition.r`（**本脚本的必需配套脚本**，须与 py 同目录）渲染，脚本会自动调用，无需手动运行。

### 参数

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

- **样本量是硬门槛**：scCODA 的经验要求是每组 ≥ 6–8 个样本。样本少时，自动参考类型选择（要求「出现率 ≥ 95% 的样本」）会退化成「几乎所有样本都要出现该类型」，几乎找不到合格参考，模型本身也不稳定。**样本量不满足时，任何方法都给不出可靠结论**——此时只能做描述性比较（占比表、柱状图、桑基图），不要下显著性结论。
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

**回答的问题：这个 h5ad 里的矩阵是原始计数吗？能不能直接拿来算差异分析？**

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

**回答的问题：怎么把 h5ad 文件变成 Seurat 能直接用的 rds？**

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

**回答的问题：怎么把单个 h5ad / rds 互相转换？**

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

**回答的问题：从 Broad Single Cell Portal 下载的矩阵和元数据，怎么直接变成可用的 Seurat 对象？**

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

## 13 scRNA_subcelltype_clustering.r（独立脚本）—— 细胞类型内亚型聚类

### 功能

**回答的问题：同一个细胞类型内部还能不能再细分？每一群亚型的特征基因是什么？**

在**已完成细胞类型注释**的对象上，对每个细胞类型**单独**再做一轮亚型细分：

1. **逐细胞类型重跑降维与聚类**（探索阶段）：对每个细胞类型取子集 → `NormalizeData` → `FindVariableFeatures(vst, 2000)` → `ScaleData` → `RunPCA` → `RunHarmony(group.by.vars = "orig.ident")` → `FindNeighbors` → `FindClusters`（Louvain 扫 0–1、Leiden 扫 0.1–1，步长 0.1）→ 输出该细胞类型的 **clustree** 供人工挑选分辨率。同时保存 `<celltype>_clustering.rds`（已存在则跳过该细胞类型）。
2. **按指定分辨率定亚型并做 DE**（应用阶段）：读取 `-Res` 指定的「细胞类型 → 分辨率」CSV，对每个细胞类型用该分辨率重新 `FindClusters`，生成 `sub_celltype`（形如 `<celltype>_<cluster id>`），跑 `FindAllMarkers`（`only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25`），导出全部 marker、top10 marker 与每亚型 dotplot。
3. **回写原对象**：把各亚型标签合并写回原始 Seurat 对象的 `sub_celltype` 列，输出整体 UMAP 与 `<前缀>subcelltype.rds`。

### 输入

- `-I` 已注释的 Seurat rds（需含细胞类型列，且带 `umap` reduction）
- `-AC` 该细胞类型列名（必填）
- `-Res` 分辨率表 CSV（**第二阶段必填**），两列：`celltype`、`resolution`

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | （必填） | 注释后的 Seurat rds |
| `-AC, --annotation_column` | （必填） | metadata 中细胞类型列名 |
| `-Res, --resolution_file` | 无 | 分辨率表 CSV（`celltype`, `resolution` 两列）|
| `-PC` | 30 | PCA 维数 |
| `-CA, --cluster_algorithm` | Louvain | `Louvain` 或 `Leiden`（Leiden 需 `-P_env` 或已装 `leidenalg`）|
| `-A, --assay` | RNA | 表达所在 assay |
| `-P_env` | 无 | Leiden 用的 Python 环境路径 |
| `-T, --threads` | 10 | 线程数（Harmony 与 future）|
| `-Red, --reduction` | harmony | 聚类所用降维（`pca` 或 `harmony`）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 阶段 1：逐细胞类型出 clustree，人工挑分辨率
./scRNA_subcelltype_clustering.r -I 05_annotation/05.seu_annotated.rds -AC celltype     -CA Leiden -T 12 -P 06. -O 06_subtype

# 阶段 2：填好分辨率表后定亚型、跑 DE
./scRNA_subcelltype_clustering.r -I 05_annotation/05.seu_annotated.rds -AC celltype     -Res resolution_recording.csv -CA Leiden -T 12 -P 06. -O 06_subtype
```

分辨率表格式：

```csv
celltype,resolution
celltypeA,0.4
celltypeB,0.8
```

### 输出

输出目录 `-O` 下（`<ct>` 为细胞类型名，空格会替换为下划线）：

| 文件 | 说明 |
|---|---|
| `<前缀><ct>_clustering.rds` | 该细胞类型的多分辨率聚类对象（阶段 1 产物，已存在则跳过）|
| `<前缀><ct>_clustree.png/.pdf` | 该细胞类型的多分辨率聚类树 |
| `<前缀><ct>_res<res>_umap.png/.pdf` | 该细胞类型在指定分辨率下的 UMAP |
| `<前缀><ct>_clustering_res<res>_results.rds` | 该细胞类型定亚型后的对象（含 `sub_celltype`）|
| `<ct>_res<res>_DE/<前缀><ct>_markers.csv` | 该细胞类型内各亚型的全部 marker |
| `<ct>_res<res>_DE/<前缀><ct>_top10_markers.csv` | 每亚型 top 10 marker（`p_val_adj<0.05` 且 `pct.1>0.3`）|
| `<ct>_res<res>_DE/<前缀><亚型名>_dotplot.png/.pdf` | 每亚型一张 top10 marker dotplot |
| `<前缀>subcelltype_umap.png/.pdf` | 全对象按 `sub_celltype` 着色的 UMAP |
| `<前缀>subcelltype.rds` | **主输出**：原始对象 + `sub_celltype` 列 |

示例输出：`Images/subcelltype_clustering/`（示例为单个样本的完整一轮：17 个细胞类型的 clustree / 亚型 UMAP + 各亚型 dotplot）。

### 注意事项

- **细胞类型子集按 `-AC` 指定的列取**：阶段 1 用 `seo[, seo@meta.data[[annotation_column]] == ct]` 逐类型切分，所以对象里的细胞类型列叫什么都能用。
- **但阶段 2 的分辨率表，其列名必须叫 `celltype`**：脚本用 `celltype_res$celltype` 读取该表（第 230 行），表的这个列名是写死的；表里另一列必须叫 `resolution`。
- **分辨率表按替换后的名字匹配**：表里的 `celltype` 值要与**空格换成下划线后**的名字一致（如 `T cell` 写作 `T_cell`）。
- **两阶段是分开跑的**：阶段 1 只需 `-AC`；阶段 2 必须给 `-Res`，否则脚本报 `Please provide a valid resolution file using the -Res parameter.` 后中止。阶段 1 的 `<ct>_clustering.rds` 是阶段 2 的输入——若删掉它，阶段 2 会因找不到文件报错。
- **细胞类型名含空格会被替换**：文件名前缀里的空格统一换成 `_`。
- Harmony 固定按 `orig.ident` 分组校正，与 `-Red` 无关；`-Red` 只影响 `FindNeighbors` 用哪个 reduction。
- Leiden 分两支：给了 `-P_env` 就 `use_python()` 并检查 `leidenalg`；没给则 `reticulate::py_require("leidenalg")` 现场装（需可用的 Python + 编译工具）。
- 亚型 DE 的 `FindAllMarkers` 是耗时大头；亚型多时逐个 dotplot 渲染也较慢。
- `subset()` 保留原对象的 `umap` reduction（继承自原对象），因此第 2 阶段的 UMAP 与原注释 UMAP 同一坐标系。


---

## 14 scRNA_pseudobulk_DE.r（独立脚本）—— 按细胞类型找两组差异基因（pseudobulk）

### 功能

**回答的问题：在每一种细胞类型内部，两组（`-G1` vs `-G2`）之间哪些基因的表达存在差异？**

具体做法：把同一个「细胞类型 × 样本」内的所有细胞的原始 counts 相加，得到一行 pseudobulk 计数（每个细胞类型得到一张「样本 × 基因」的表），再用 DESeq2 在**样本层面**比较两组，得到该细胞类型的差异基因列表；**每个细胞类型各做一次独立检验**，最后把所有细胞类型的结果合并成一张表并画成火山图（每个细胞类型一个分区，标注该类型的 top 差异基因）。

因此输出是一张**长表**：`cluster` 列标明每个差异基因来自哪一种细胞类型，可以直接读出「哪类细胞里的哪个基因在两组间不同」。

**统计口径要点**：重复单位是**样本**（生物学重复），不是细胞。细胞层面的 Wilcoxon 检验会把每个细胞当作一个重复，p 值随细胞数增大而虚高；pseudobulk 把同一样本的细胞先合并，p 值反映的才是「换个样本还看不看得到这个差异」，因此结论可外推到新样本。

### 输入

一个已注释的 Seurat rds，**必须含 `counts` 层**（原始整数计数）。若对象是多个对象 `merge()` 得来、counts 被拆成 `counts.1`/`counts.2` 等多层，脚本会自动 `JoinLayers` 合并后再聚合。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| 已注释 Seurat rds（含 counts 层）|
| `-CC, --celltype_column` | 无（必填）| 细胞类型列；逐个水平各做一次检验 |
| `-SC, --sample_column` | 无（必填）| 样本/个体列，即生物学重复 |
| `-GC, --group_column` | 无（必填）| 分组列 |
| `-G1, --group1` | 无（必填）| 作为 ident.1 的组；log2FC > 0 表示在该组更高 |
| `-G2, --group2` | 无（必填）| 作为 ident.2（参照）的组 |
| `-A, --assay` | RNA | 含 counts 层的 assay |
| `-MC, --min_cells` | 20 | 单个「细胞类型 × 样本」单元的最少细胞数，低于则丢弃该单元 |
| `-MS, --min_samples` | 3 | 每组最少样本数，不足则该细胞类型跳过（不检验）|
| `-p_val` | 0.05 | 原始 p 值阈值（仅用于过滤输出）|
| `-p_adj` | 0.05 | 校正 p 值阈值（仅用于过滤输出）|
| `-lgFC` | 1 | \|avg_log2FC\| 阈值（仅用于过滤输出）|
| `-top_n` | 5 | 火山图每细胞类型每侧标注的基因数 |
| `-W, --plot_width` | 18 | 图宽（英寸）|
| `-H, --plot_height` | 12 | 图高（英寸）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_pseudobulk_DE.r -I 05_annotation/05.seu_annotated.rds \
    -CC celltype -SC orig.ident -GC group -G1 treat -G2 control \
    -P 06. -O 06_pseudobulk
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀><G1>_vs_<G2>_DE_raw.csv` | 全部基因 × 全部细胞类型的检验结果，未过滤 |
| `<前缀><G1>_vs_<G2>_DE_filtered.csv` | 按 p / 校正 p / log2FC 阈值过滤后的差异基因 |
| `<前缀><G1>_vs_<G2>_pseudobulk_units.csv` | **溯源表**：每个细胞类型的单元数、两组样本数、单元细胞数的最小/最大值 |
| `<前缀><G1>_vs_<G2>_DE_volcano.png/.pdf` | 火山图（灰阴影柱表示每个细胞类型的 log2FC 范围，红上/蓝下，标注 top 基因，底部细胞类型色块）|

示例输出：`Images/pseudobulk_DE/`。

CSV 列名固定为 `gene,p_val,avg_log2FC,pct.1,pct.2,p_val_adj,cluster`（`cluster` 即细胞类型），与 Seurat `FindAllMarkers` 的输出格式一致，便于下游工具直接读取。

`pct.1` / `pct.2` 在 pseudobulk 语境下的含义是：**该组内有多少比例的样本（pseudobulk 单元）检出该基因**（非零计数），而不是细胞比例。

### 注意事项

- **必须用原始 counts**：脚本从 `counts` 层聚合；若对象只有 `data`（归一化后的值），pseudobulk 求和没有意义，DESeq2 也会拒绝非整数输入。
- **不做 CPM 归一化**：聚合方式是原始 counts 直接相加，样本间测序深度差异由 DESeq2 的 size factor（median-of-ratios）吸收。这样保持整数计数以符合负二项模型，也不引入「总量归一化」的额外假设。
- **重复单位是样本**：`-SC` 必须是**生物学重复**（个体/样本），不能填细胞条形码或分组标签，否则退化成伪重复。
- **`-MC` 是数据质量闸门**：单元细胞数太少时聚合噪声大。典型 10–50 的区间；调高会减少可用单元数、调低会放进噪声单元。经验规律：`-MC` 偏低时个别细胞类型会因某组样本数不足而被跳过；调得过高则细胞类型陆续退出检验。默认值 20 在多数数据上兼顾了保留单元数与控制噪声。
- **`-MS` 是统计可行性闸门**：不足 3 个样本/组时 DESeq2 的离散度估计不稳定，该细胞类型直接跳过并打印 `>>> Skip ...`，不产出 p 值。
- **过滤阈值只影响输出文件与图**，不影响检验本身：`raw` 表保留全部基因的原始 p 值，可自行用别的阈值重筛。
- **BH 校正是在每个细胞类型内做的**（每个细胞类型一次 DESeq2 运行），不是跨所有细胞类型合并校正。若需要全局 FDR，请对 raw 表自行重算。
- 细胞类型名不做任何字符替换，含括号或空格的名称会原样出现在 CSV 与图中。
- 火山图使用 `pseudo_log` 变换的 y 轴（log2FC 而非 -log10 p），与库内其它 DE 脚本风格一致；因此**点的位置只反映效应量与显著性方向的组合，不反映 p 值大小**。


---

## 15 scRNA_spRNA_cellPro.r（独立脚本）—— 细胞类型占比分析（每样本构成 + 组间比较）

### 功能

**回答的问题：每个样本里各种细胞类型各占多少？两组之间这些比例有没有差别？**

脚本做四件事，全部基于已注释对象的 `meta.data`（只用细胞计数，不需要表达矩阵）：

1. **计数与占比表**：按样本统计每个细胞类型的细胞数，以及该细胞类型在这个样本中的占比；再按组汇总一份。
2. **100% 堆叠柱状图**：每个样本一根柱子，按组分面，一眼看出组间构成差异。
3. **组成桑基图**：组 → 细胞类型的流向图，看整体构成比例。
4. **组间占比比较**（需 `-GO` 恰好给两组）：逐细胞类型画箱线图 + 散点 + Wilcoxon 检验与显著性标注；散点**每个点是一个样本**。
5. **亚型层面的占比比较**（可选 `-SCC`）：若对象有亚型列（如 `sub_celltype`），再按每个父细胞类型出一张亚型占比比较图，最后拼成一张大图。

### 输入

已注释的 Seurat rds（只需 `meta.data` 中有样本列、细胞类型列、分组列；亚型比较另需亚型列）。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| 已注释 Seurat rds |
| `-SC, --sample_column` | 无（必填）| 样本身份列（如 `orig.ident`），决定统计单位 |
| `-CC, --celltype_column` | 无（必填）| 细胞类型列（如 `celltype`）|
| `-GC, --group_column` | 无（必填）| 分组列（如 `Group`）|
| `-SCC, --sub_celltype_column` | 无 | 亚型列（如 `sub_celltype`）；给了才做亚型占比比较 |
| `-GO, --group_order` | 无 | 分组显示顺序，逗号或空格分隔（如 `control,treat`）。**只有恰好两组时才做占比检验**，并决定柱状图分面顺序与桑基图上下位置 |
| `-EC, --end_celltype` | 无 | 指定某个细胞类型钉在图的最末端（常用于把 `Unknown` 放到最后）|
| `-BW / -BH` | 10 / 12 | 柱状图宽 / 高（英寸）|
| `-SW / -SH` | 12 / 6 | 桑基图宽 / 高（英寸）|
| `-DW / -DH` | 14 / 8 | 细胞类型占比比较图宽 / 高（英寸）|
| `-SDW / -SDH` | 36 / 28 | 亚型占比比较大图宽 / 高（英寸）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 计数 + 占比表 + 柱状图 + 桑基图 + 两组占比比较
./scRNA_spRNA_cellPro.r -I 06_annotation/06.seu_annotated.rds \
    -SC orig.ident -CC celltype -GC group -GO control,treat -EC Unknown \
    -P 07. -O 07_cellPro

# 追加亚型层面比较（对象需有 sub_celltype 列）
./scRNA_spRNA_cellPro.r -I 06_annotation/06.seu_annotated.rds \
    -SC orig.ident -CC celltype -SCC sub_celltype -GC group -GO control,treat \
    -P 07. -O 07_cellPro
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀><细胞类型列>_count_per_sample.csv` | 每个样本 × 每种细胞类型的细胞数 |
| `<前缀><细胞类型列>_proportion_per_sample.csv` | 每个样本 × 每种细胞类型的占比（行内合计为 1）|
| `<前缀><细胞类型列>_proportion_stat.csv` | 长表：每行 = 一个（样本 × 细胞类型），含细胞数、占比、所属组 |
| `<前缀><细胞类型列>_count_per_group.csv` | 每组 × 每种细胞类型的细胞数 |
| `<前缀><细胞类型列>_proportion_per_group.csv` | 每组 × 每种细胞类型的占比 |
| `<前缀>bar.png/.pdf` | 100% 堆叠柱状图（按组分面）|
| `<前缀>sankey.png/.pdf` | 组 → 细胞类型桑基图 |
| `<前缀><细胞类型列>_DE.png/.pdf` | 逐细胞类型的两组占比比较图（需 `-GO` 恰好两组）|
| `<前缀><亚型列>_proportion_stat.csv` | 亚型层面的占比长表（需 `-SCC`）|
| `<前缀><亚型列>_DE.png/.pdf` | 各父细胞类型下的亚型占比比较大图（需 `-SCC`，3 列拼接）|

### 注意事项

- **统计单位是样本**：占比按样本内部归一化（该类型细胞数 ÷ 该样本总细胞数），散点与箱线图里**每个点是一个样本**，不是细胞。
- **`-GO` 决定是否做检验**：只有 `-GO` 恰好给出两个组时，才输出 `<细胞类型列>_DE` 图并做 Wilcoxon 检验；`-GO` 给三个以上组时，柱状图与桑基图仍会画（按给定顺序分面/上下排列），但不做组间检验。
- **`-GO` 同时决定显示顺序**：柱状图按组分面时以 `-GO` 顺序排列；桑基图用 `rev(group_order)`，即**对照组在上、处理组在下**。
- **细胞类型排序按全局丰度降序**（`-EC` 指定的类型钉到末尾），柱状图、桑基图、比较图都用同一顺序，便于跨图对照。
- **子集过滤**：给了 `-GO` 时，metadata 里不属于这些组的细胞会被**直接剔除**，因此计数表只含 `-GO` 覆盖的组。
- **`-EC` 只影响绘图顺序**，不影响任何数值（占比仍按全部细胞类型计算）。
- **Wilcoxon 检验用的是每样本一个比例值**，因此**检验效能完全由每组的样本数决定**：秩和检验在小样本下能达到的最小 p 值有硬下限——每组 3 个样本时最小双侧 p = 0.1，每组 4 个时 = 0.0286，每组 5 个时 = 0.0079。也就是说**每组只有 3 个样本时，这张图不可能出现 `*` 或 `**` 标注**，此时它只能当「差异排序」的参考，不能作为显著性证据。
- **不显著 ≠ 没有差异**：样本量不足导致的 p 值偏大，不能解读为「两组无差异」。
- 显著性标注用 `p.signif`（`ns` / `*` / `**` …），`hide.ns = FALSE`，所以不显著的也会标 `ns`。
- **箱线图在比例极小的类型上会退化**：某组完全缺失该类型时该组只剩一个点，箱体不可见。
- **所有宽高参数都正常生效**：`-BW`/`-BH`（堆叠柱状图）、`-SW`/`-SH`（桑基图）、`-DW`/`-DH`（细胞类型比较图）、`-SDW`/`-SDH`（亚型比较图），单位均为英寸。
- 运行时会出现 ggplot2 的 `size` → `linewidth` 弃用警告（`element_line(size=...)`、桑基 `geom_flow(size=...)`），属正常提示，不影响输出。

### 与 `scRNA_composition_sccoda.py` 的区别

两者都分析「细胞组成比例」，但定位不同：

| | `scRNA_spRNA_cellPro.r`（本脚本）| `scRNA_composition_sccoda.py` |
|---|---|---|
| 定位 | 快速描述 + 简单 Wilcoxon 检验 | 正式的成分差异分析（Dirichlet-Multinomial 贝叶斯模型）|
| 输出 | 计数/占比表 + 柱状图 + 桑基图 + 箱线图 | 可信效应表 + 离散度图 + 效应森林图 + 桑基 |
| 统计 | 每样本一个比例值 → Wilcoxon | 计数表 → 贝叶斯后验，**给出可信区间，能处理「一类升高导致另一类机械性下降」的成分数据假阳性** |
| 样本量要求 | 极低（但**每组 ≤ 3 时无法给出显著结果**）| 较高（经验上每组 ≥ 6–8 个样本）|
| 依赖 | 纯 R | Python + scCODA + R（桑基图）|

成分数据有三个固有性质，决定了这两个脚本的分工：

1. **占比合计恒为 1**：某类型上升必然挤掉其它类型，n 个类型只有 n-1 个独立自由度。逐个独立做检验会把同一次变化拆成「A 升 + B 降」两个发现；scCODA 用参考细胞类型显式处理这个约束。
2. **比例 ≠ 数量**：两组捕获细胞总数不同时，占比不变但绝对数量可能成倍变化，单看占比会漏掉。
3. **检验效能由样本量决定**：本脚本的 Wilcoxon 每组样本数少时连 `p<0.05` 都达不到。

建议：**先用本脚本快速看构成、排出可疑的细胞类型，再用 scCODA 对这些类型做正式检验**（前提是样本量满足 scCODA 的要求；不满足时任何方法都给不出可靠结论）。


---

## 16 scRNA_SingleR.r（独立脚本）—— 用参考数据集自动推断细胞类型

### 功能

**为细胞群自动推断细胞类型名称**：把查询对象的每个聚类与一份参考数据集做 SingleR 比对，得到「cluster → cell type」的映射，并把推断结果写回 Seurat 对象的 `SingleR` 列。

与 step5（`scRNA_annotate_step5.r`）的区别：step5 是**基于 marker 的人工注释**（你给 cluster→cell type 映射表，脚本负责画图与写对象）；本脚本是**基于参考数据集的自动推断**（参考数据自带标签，脚本去匹配）。两者输出同名信息，但本脚本不需要你预先知道 marker，适合快速摸底或对已有注释做交叉验证。

参考数据三种来源（**互斥，必须给其中一种**）：
- `-S human` / `-S mouse`：用 celldex 内置参考（`HumanPrimaryCellAtlasData` / `MouseRNAseqData`），首次使用会联网下载并缓存
- `-RDS`：自己的 Seurat 对象作参考，用 `-AC` 指定其中哪一列是标签
- `-RD`：自己的 RData（SCE 对象）作参考，用 `-AC` 指定 `colData` 里哪一列是标签

### 输入

- 查询对象：已聚类的 Seurat rds，**必须含 `seurat_clusters` 列**（注释是按 cluster 给出的，不逐细胞）
- 参考数据：`-S` 内置参考，或 `-RDS` / `-RD` 提供的参考文件

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| 查询用的 Seurat rds，须含 `seurat_clusters` |
| `-S, --species` | 无 | 内置参考的物种：`human` 或 `mouse` |
| `-RD, --rdata_file` | 无 | RData 格式的参考（SCE 对象）|
| `-RDS, --rds_file` | 无 | RDS 格式的参考（Seurat 对象）|
| `-AC, --annotation_column` | 无 | 参考数据中存放标签的列名；用 `-RD`/`-RDS` 时必填 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 用 celldex 内置小鼠参考
./scRNA_SingleR.r -I 04_clustering/04.seu_res0.6_clustered.rds -S mouse -P 05. -O 05_SingleR

# 用自己的 Seurat 参考（-AC 指定标签列）
./scRNA_SingleR.r -I 04_clustering/04.seu_res0.6_clustered.rds -RDS ref_annotated.rds -AC celltype -P 05. -O 05_SingleR

# 用自己的 RData 参考
./scRNA_SingleR.r -I 04_clustering/04.seu_res0.6_clustered.rds -RD ref.rdata -AC label.main -P 05. -O 05_SingleR
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>SingleR_result.csv` | cluster → cell type 映射表（两列：`Cluster`、`celltype`）|
| `<前缀>singleR_heatmap.png/.pdf` | SingleR 得分热图（各 cluster 与各参考标签的相似度）|
| `<前缀>singleR_umap.png/.pdf` | UMAP 按推断出的 `SingleR` 分组着色 |
| `<前缀>seo_SingleR.rds` | 查询对象 + `SingleR` 元数据列 |
| `<前缀>SingleR_ref.rdata` | **仅 `-RDS` 模式**：把 RDS 参考转成 RData 存档，后续可直接用 `-RD` 复用，不必重复转换 |

示例输出：`Images/SingleR/`。

### 注意事项

- **三种参考来源互斥**：`-S`、`-RD`、`-RDS` 只能给一个；都不给会报错，给多个也会报错。
- **注释粒度是 cluster 不是 cell**：脚本用 `seo$seurat_clusters` 作为分组依据，所以输出是「每个聚类一个细胞类型」。输入对象若没有这一列会直接失败。
- **必须已有聚类结果**：本脚本不做聚类，期望输入是 step4 的产物（含 `seurat_clusters` 与 UMAP）。
- **UMAP 图要求对象里有 umap reduction**：脚本用 `DimPlot(reduction = "umap")` 出图，输入对象没有 UMAP 会报错。
- **查询对象用的是 counts 层**（`GetAssayData(..., layer = "counts")`），参考对象则优先用 counts（经 `scater::logNormCounts` 归一化）、无 counts 时退回 data 层。
- **`-RDS` 参考模式会做一次格式转换**：参考对象按 counts（经 `logNormCounts` 归一化，无 counts 时退回 data 层）转成 SingleCellExperiment，并**保持稀疏矩阵不稠密化**——大参考库请务必用这一模式，不要事先 `as.data.frame()`。
- **`-RD` 模式下的标签列**：脚本把 `-AC` 指定的列**复制一份**到标准列名 `label.main` 上（SingleR 只读这一列）。**原始 rdata 文件不会被改写**；但若参考文件里本来就有 `label.main` 列，它在本次运行中会被替换成 `-AC` 指向的列——需要保留原标签时，请让参考文件只含一个标签列并直接用它做 `-AC`。
- **`-AC` 缺省时的行为**：用 `-RD`/`-RDS` 却不给 `-AC` 时，脚本会先打印参考的 metadata 再报错退出，便于你查看有哪些列可用。
- **自动推断不等于最终注释**：SingleR 只看表达谱相似度，可能给出与实际不符的名字（例如把相近的免疫亚群混淆）。正式结果建议用 step5 的 marker 注释，或至少与 marker 结果交叉核对。
- **参考标签的粒度决定输出粒度**：参考数据里的标签是细胞类型级还是更细的亚型级，输出就是对应粒度。
- 运行时会出现 ggplot2 `size` → `linewidth` 的弃用警告，属正常提示，不影响输出。


---

## 17 scRNA_spRNA_enrichment.r（独立脚本）—— 差异基因的 GO/KEGG 富集与 GSEA

### 功能

**把差异基因翻译成通路层面的结论**：对每个细胞类型，回答「这个类型里变化的基因主要参与哪些通路？」

脚本做三步：

1. **准备背景基因**：按细胞类型取出该类型中**至少 1 个细胞表达 >1 次**的基因，作为该类型的背景基因集（比「全部基因」作背景更贴合单细胞情境）
2. **准备差异基因**：从 `-DER` 指定的 DE 结果表里按 `p_val_adj < -p` 且 `|avg_log2FC| > -log2FC` 过滤，再按 `cluster` 列分配到各细胞类型
3. **调用富集程序**：对每个细胞类型分别调用外部命令 `enrichment_analysis`，做
   - **GO / KEGG 富集**：以该类型自身检出的基因为背景，检验差异基因在通路上的富集
   - **GSEA**：用该类型**全部基因**按 `-GM` 指定的指标（默认 `avg_log2FC`）排序做预排序富集，不依赖差异基因阈值

### 重要前提：依赖配套脚本 `enrichment_analysis.r`

**本脚本本身不做富集计算**，它只准备输入文件并通过 `system()` 调用 `enrichment_analysis` 这个命令名。该命令由**配套脚本 `enrichment_analysis.r`**（见第 18 章）提供，内部使用 `clusterProfiler`（`enrichGO` / `enrichKEGG` / `gseGO` / `gseKEGG`）、`enrichplot`、`GseaVis` 和物种对应的注释库 `org.Hs.eg.db` / `org.Mm.eg.db`。

- **命令名不带 `.r` 后缀**：脚本里写的是 `enrichment_analysis`，所以需要让它在 `PATH` 中可用——把 `enrichment_analysis.r` 建立一个名为 `enrichment_analysis` 的软链（或副本）即可，例如：
  ```bash
  ln -s "$(pwd)/enrichment_analysis.r" ~/bin/enrichment_analysis
  ```
- 若该命令不在 `PATH` 中，脚本仍会正常跑完并生成所有中间文件，但**富集结果为空**，且每个 `<细胞类型>_run.log` 里会记录「command not found」
- 脚本会把每次调用的完整命令**写成 `<细胞类型>_run.sh` / `<细胞类型>_GSEA_run.sh`** 保存下来，便于手工重跑单个细胞类型

### 输入

- `-I`：已注释的 Seurat rds（用于提取背景基因）
- `-DER`：DE 结果 csv，需含 `gene`、`cluster`、`p_val_adj`、`avg_log2FC` 列（与 `scRNA_pseudobulk_DE.r` 的输出格式兼容）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| 已注释 Seurat rds |
| `-DER, --DE_result` | 无（必填）| DE 结果 csv（raw，未过滤）|
| `-CC, --celltype_column` | 无（必填）| 细胞类型列（与 DE 表的 `cluster` 列对应）|
| `-GC, --group_column` | 无（必填）| 分组列（用于筛出参与比较的组）|
| `-G, --groups` | 无 | 参与比较的组（逗号或空格分隔）；只用这些组的细胞算背景基因 |
| `-A, --assay` | RNA | 取 counts 层的 assay |
| `-p, --p_threshold` | 0.05 | DE 基因过滤：校正 p 值上限 |
| `-log2FC, --log2FC_threshold` | 0 | DE 基因过滤：\|avg_log2FC\| 下限 |
| `-GKP, --GO_KEGG_p_threshold` | 0.05 | GO/KEGG 富集的 p 值阈值 |
| `-GKQ, --GO_KEGG_q_threshold` | 0.2 | GO/KEGG 富集的 q 值（FDR）阈值 |
| `-GM, --GSEA_metric` | avg_log2FC | GSEA 的排序指标（DE 表里的列名）|
| `-GSEA_P` | 0.05 | GSEA 的 p 值阈值 |
| `-GSEA_Q` | 0.25 | GSEA 的 q 值阈值（GSEA 惯例 0.25）|
| `-NE, --GSEA_NES` | 1 | GSEA 保留 \|NES\| 大于该值的通路 |
| `-S, --species` | human | 物种：`human` 或 `mouse`，决定用哪个注释库 |
| `-IC, --ignore_celltype` | 无 | 要排除的细胞类型（逗号或空格分隔），如 `Unknown` |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_spRNA_enrichment.r -I 05_annotation/05.seu_annotated.rds \
    -DER 06_pseudobulk/06.treat_vs_control_DE_raw.csv \
    -CC celltype -GC group -G control,treat -A RNA -S mouse \
    -IC Unknown -P 07. -O 07_enrichment
```

### 输出

输出目录 `-O` 下（`<ct>` = 细胞类型）：

| 文件 | 说明 |
|---|---|
| `<前缀><ct>_bg_genes.csv` | 该细胞类型的背景基因（单列 `Gene`）|
| `<前缀><ct>_DEGs.txt` | 该细胞类型的差异基因（每行一个）|
| `<前缀><ct>_GSEA_input.csv` | 该细胞类型 GSEA 的输入表（行名为基因，含排序指标列）|
| `<前缀><ct>_run.sh` / `<前缀><ct>_GSEA_run.sh` | 对应的外部命令，便于手工重跑 |
| `<前缀><ct>_run.log` / `<前缀><ct>_GSEA_run.log` | 外部命令的运行日志 |
| `result/<前缀><ct>_GO_result.csv` | GO 富集结果表 |
| `result/<前缀><ct>_GO_dotplot.png/.pdf`、`_GO_barplot.png/.pdf` | GO 富集图 |
| `result/<前缀><ct>_KEGG_result.csv` | KEGG 富集结果表 |
| `result/<前缀><ct>_KEGG_dotplot.png/.pdf`、`_KEGG_barplot.png/.pdf` | KEGG 富集图 |
| `result/<前缀><ct>_GSEA_GO_raw_result.csv`、`_GSEA_GO_filtered_result.csv` | GSEA 的 GO 结果（全部 / 按阈值过滤）|
| `result/<前缀><ct>_GSEA_GO.png/.pdf` | GSEA 的 GO 图 |
| `result/<前缀><ct>_GSEA_KEGG_raw_result.csv`、`_GSEA_KEGG_filtered_result.csv` | GSEA 的 KEGG 结果 |
| `result/<前缀><ct>_GSEA_KEGG.png/.pdf` | GSEA 的 KEGG 图 |

示例输出：`Images/enrichment/`。

`result/` 子目录由脚本自动创建，所有富集产物都写在里面。

### 注意事项

- **细胞类型数不足 10 个则跳过**：背景基因少于 10 个、或差异基因少于 10 个的细胞类型会被丢弃（打印 `less than 10, discard it !`），对应的文件不会生成。
- **背景基因的过滤条件是 `rowSums(counts > 1) > 0`**：即该基因至少在一个细胞里表达 >1 次。
- **GO/KEGG 与 GSEA 的输入不同**：GO/KEGG 用「该类型的差异基因 vs 该类型的背景基因」；GSEA 用「该类型全部基因按 log2FC 排序」，**不经过差异基因阈值过滤**。所以 GSEA 能检出那些「整体一致但幅度不大」的通路，与 GO 富集互补。
- **DE 表必须与对象对齐**：`-DER` 的 `cluster` 列要与 `-CC` 指定的列取值一致，否则该类型的差异基因数为 0 而被跳过。
- **物种必须在 DE 结果准备好之前确定**：`-S` 决定注释库（human→`org.Hs.eg.db`，mouse→`org.Mm.eg.db`），选错会导致基因 ID 全部匹配不上、结果为空。
- **富集是后台并行的**：外部命令以 `&` 结尾并发执行，脚本本身会先返回。**脚本显示「All Done」不代表富集已跑完**——需等 `result/` 里的文件出现，或查看 `*_run.log`。
- **KEGG 需要联网**：KEGG 富集会调用 KEGG REST API；离线环境只有 GO 部分可用。
- **`-G` 指定的组之外**的细胞不参与背景基因统计，但不影响差异基因的取用（差异基因来自 DE 表）。


---

## 18 enrichment_analysis.r（配套脚本）—— GO/KEGG 富集与 GSEA 的实际执行者

### 功能

**执行通路富集计算本身**：`scRNA_spRNA_enrichment.r` 只负责准备输入文件，真正的富集由本脚本完成。本脚本面向的是「一组基因 → 这些基因富集在哪些 GO 术语 / KEGG 通路上」这一问题，既可用于单细胞（逐细胞类型调用），也可直接用于 bulk RNA。

两种分析模式：

- **GO / KEGG 富集（`-G`）**：给定目标基因列表，用 `enrichGO` / `enrichKEGG` 检验它们在通路上的富集程度。给 `-E` 时以该文件里的基因为背景（universe），不给则用注释库的默认背景。
- **GSEA（`-GSEA` + `-GM`）**：给定全基因排序表（如按 log2FC 排序），用 `gseGO` / `gseKEGG` 做预排序富集，输出 NES 与 p/q 值。**不依赖差异基因阈值**，能检出「整体一致变化但幅度不大」的通路。

### 输入

两种分析模式对应两种输入（`-G` 与 `-GSEA` 互斥）：

| 模式 | 输入 | 含义 |
|---|---|---|
| GO / KEGG 富集 | `-G` 基因列表文件（每行一个基因符号，至少 10 个）| 要检验这些基因在通路上的富集；`-E` 可另给背景基因文件（行名为基因符号），不给则用注释库默认背景 |
| GSEA | `-GSEA` 排序表（行名为基因符号）+ `-GM` 排序指标列名 | 全基因按某指标排序后做预排序富集，不依赖差异基因阈值 |

### 重要前提

- **`-S` 决定物种**：human → `org.Hs.eg.db`，mouse → `org.Mm.eg.db`。选错会导致基因 ID 全部匹配不上、结果为空。
- **`-M` 可切换 KEGG 数据源**：默认走在线 KEGG REST API；加 `-M` 则改用 MSigDB 内置的 KEGG 基因集（v2025.1），**网络不稳定或离线时推荐加这个参数**。
- **KEGG 在线模式需要联网**。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-G, --genes` | 无 | 目标基因文件（每行一个基因，至少 10 个）。与 `-GSEA` 互斥 |
| `-E, --expr` | 无 | 背景基因文件（行名为基因符号）。给了 `-E` 就必须给 `-G` |
| `-S, --species` | human | `human` 或 `mouse` |
| `-p, --pvalue` | 0.05 | 显著通路的 p 值阈值 |
| `-Q, --qvalue` | 0.2 | 显著通路的 q 值（FDR）阈值；GSEA 建议 0.25 |
| `-GSEA, --GSEA_file` | 无 | GSEA 输入表（行名为基因符号）。必须与 `-GM` 同时给 |
| `-GM, --GSEA_metric` | 无 | GSEA 的排序指标列名 |
| `-NE, --NES_threshold` | 1 | GSEA 保留 \|NES\| 大于该值的通路 |
| `-M, --MSigDB` | FALSE | 改用 MSigDB 的 KEGG 基因集，不走在线 KEGG |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# GO/KEGG 富集（自带背景基因）
enrichment_analysis -G genes.txt -E bg_genes.csv -S mouse -P 01. -p 0.05 -Q 0.2 -O result/

# GSEA（离线时加 -M）
enrichment_analysis -GSEA GSEA_input.csv -GM avg_log2FC -S mouse -NE 1 -p 0.05 -Q 0.25 -M -P 02. -O result/
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>GO_result.csv` | GO 富集结果表（BP / MF / CC 合并）|
| `<前缀>GO_dotplot.png/.pdf`、`<前缀>GO_barplot.png/.pdf` | GO 富集图（按 ontology 分面）|
| `<前缀>KEGG_result.csv` | KEGG 富集结果表 |
| `<前缀>KEGG_dotplot.png/.pdf`、`<前缀>KEGG_barplot.png/.pdf` | KEGG 富集图 |
| `<前缀>GSEA_GO_raw_result.csv` | GSEA 的 GO 全部结果 |
| `<前缀>GSEA_GO_filtered_result.csv` | 按 `-p`/`-Q`/`-NE` 过滤后的 GSEA GO 结果 |
| `<前缀>GSEA_GO.png/.pdf` | GSEA 的 GO 图 |
| `<前缀>GSEA_KEGG_raw_result.csv`、`<前缀>GSEA_KEGG_filtered_result.csv` | GSEA 的 KEGG 结果（全部 / 过滤后）|
| `<前缀>GSEA_KEGG.png/.pdf` | GSEA 的 KEGG 图 |

### 注意事项

- **`-G` 与 `-GSEA` 互斥**，同时给会报错；`-E` 不能单独给。
- **`-GSEA` 与 `-GM` 必须成对出现**，只给其一也会报错。
- **基因标识符用 SYMBOL**：注释库的 `keyType = "SYMBOL"`，输入文件必须是基因符号，不能是 Ensembl ID 或 Entrez ID。
- **目标基因少于 10 个时结果无意义**，脚本不做检查，需自行保证。
- **GO 富集做 BH 校正**（`pAdjustMethod = "BH"`），`universe` 背景只在给了 `-E` 时生效。
- **`-E` 的背景文件格式**：第一行是表头、每行一个基因符号，读取时用行名（`read.csv(row.names = 1)`）。`scRNA_spRNA_enrichment.r` 生成的 `<细胞类型>_bg_genes.csv` 是「表头 `Gene` + 一列基因」的格式，可直接使用。
- **KEGG 的基因 ID 需先转成 ENTREZID**，`readable = TRUE` 会把结果转回符号名。
- **输出目录会自动创建**，不存在时递归建。
- 运行时会打印 ggplot2 `size` → `linewidth` 的弃用警告，属正常提示。


---

## 19 scRNA_GSVA_enrichment.r（独立脚本）—— GSVA 通路活性打分与组间比较

### 功能

**不依赖差异分析结果，直接回答「哪些通路的活性在两组之间不一样」**：对每个细胞类型，用 GSVA 把表达谱转换成通路活性分数，再比较两组之间每条通路的分数差异。

与 `scRNA_spRNA_enrichment.r` 的分工：那一个要先有 DE 结果、用「差异基因 vs 背景基因」做通路富集；本脚本**只用表达谱**，把每条通路在每个样本上的活性算出来再比较，因此不受「差异基因阈值」限制，也能检出「整体一致变化但单个基因幅度不大」的通路。

四步流程：

1. **伪 bulk 聚合**：用 `AggregateExpression` 按「样本 × 细胞类型」把 counts 求和
2. **归一化**：log2(CPM + 1)
3. **GSVA 打分**：从 msigdbr 取三套基因集——GO:BP（`C5`）、KEGG MEDICUS（`C2`）、REACTOME（`C2`），对每个细胞类型单独算 GSVA 分数
4. **组间比较**：用 limma（`~ 0 + group`）对每条通路做两组比较，输出 logFC / P.Value / adj.P.Val，并按细胞类型画通路热图

### 输入

- 已注释的 Seurat rds，需含细胞类型列、样本列、分组列三列
- 样本列必须是**生物学重复**（每个样本一个独立个体），pseudobulk 聚合与组间比较都以此为单位

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| Seurat rds |
| `-AC, --annotation_column` | 无（必填）| 细胞类型列 |
| `-SC, --sample_column` | 无（必填）| 样本列（生物学重复）|
| `-GC, --group_column` | 无（必填）| 分组列 |
| `-G, --groups` | 无 | 要比较的两组，逗号分隔。**第一个是对照组，第二个是处理组** |
| `-S, --species` | 无（必填）| `human` 或 `mouse` |
| `-A, --assay` | RNA | 使用的 assay |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_GSVA_enrichment.r -I 05_annotation/05.seu_annotated.rds \
    -AC celltype -SC orig.ident -GC group -G control,treat -S mouse \
    -P 08. -O 08_GSVA
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>GSVA_score.csv` | 全部细胞类型的 GSVA 分数长表：`celltype`、`pathway`、各样本的分数、`geneset`（该通路的基因列表）|
| `<前缀>GSVA_DE_result.csv` | 逐细胞类型的通路差异结果（limma 的 logFC / P.Value / adj.P.Val / AveExpr / t 等 + `pathway`、`celltype`、`contrast`）|
| `<前缀>GSVA_heatmap_<细胞类型>.png/.pdf` | 该细胞类型的通路热图：行=通路（按数据库与上/下调分组），列=样本（左为对照组、右为处理组），右侧标注数据库来源 |

示例输出：`Images/GSVA_enrichment/`。

### 注意事项

- **组别顺序有含义**：`-G` 的第一个是**对照组**、第二个是**处理组**，对比方向写成 `<处理组> - <对照组>`。热图的列顺序也按此排列（对照组在左）。
- **每组至少 3 个样本**：不足 3 个的细胞类型会被跳过（打印 `insufficient samples per group`）。
- **热图颜色是 GSVA 活性分数，不是差异倍数也不是 p 值**：每列是一个样本、每行是一条通路，颜色只反映该样本在该通路上的分数高低。**热图里的通路不一定显著**——筛选是三层降级（`adj.P.Val < 0.05` → `P.Value < 0.05` → 按 `|logFC|` 排序，日志会打印用了哪一层），且每库每方向只取 `|logFC|` 前 10 条；显著性只能回查 `GSVA_DE_result.csv` 的 adj.P.Val。热图的价值在于直观看组间分离趋势（如对照组一片低、处理组一片高）。
- **每个数据库各取上/下调前 10 条**：所以一张热图最多 3 个库 × 2 个方向 × 10 = 60 行。
- **样本名里的连字符会被换成下划线**：为了与 GSVA 输出矩阵的列名匹配，样本名中的 `-` 统一转成 `_`。
- **细胞类型名含下划线或空格时要注意**：脚本用「列名以 `_细胞类型` 结尾」来定位某个细胞类型的列，若细胞类型名本身以下划线结尾或含特殊字符，可能定位不准。
- **msigdbr 首次运行需要联网**下载基因集，之后会缓存。
- **GO 用的注释库是 GO:BP**（`C5` 集合的 BP 子集），不含 MF 与 CC。
- **输出的是通路活性差异，不是基因差异**：`GSVA_DE_result.csv` 里的 `pathway` 列对应通路的 MSigDB 名称（如 `GOBP_XXX`、`KEGG_XXX`、`REACTOME_XXX`）。


---

## 20 scRNA_genesetScore.r（独立脚本）—— 自定义基因集活性打分

### 功能

**给一份自己指定的基因集打分，看它在各细胞类型（或各细胞）里的活性高低**：读取基因集文件，算出 GSVA / ssGSEA / AddModuleScore / AUCell 四种分数，并画出对应的图。

与另外两个富集脚本的分工：

| 脚本 | 基因集来源 | 回答什么 |
|---|---|---|
| `scRNA_genesetScore.r` | **你自己提供的基因集文件** | 这个基因集在哪些细胞类型 / 细胞里更活跃 |
| `scRNA_GSVA_enrichment.r` | msigdbr（GO:BP / KEGG / REACTOME）| 两组之间哪些**通路活性**有差异 |
| `scRNA_spRNA_enrichment.r` | msigdbr（GO / KEGG）| 差异基因落在哪些**通路上** |

四种打分方法的分工：

| 方法 | 层级 | 特点 |
|---|---|---|
| **GSVA** | 细胞类型（pseudobulk）| 基于基因集在排序中的分布打分 |
| **ssGSEA** | 细胞类型（pseudobulk）| 与 GSVA 同为排序类方法，但走的是单样本积分 |
| **AddModuleScore** | 单细胞 | Seurat 内置，取基因集均值减去对照基因集均值 |
| **AUCell** | 单细胞 | 基于排名曲线下面积，对基因表达量级不敏感 |

前两种把细胞按类型聚合成 pseudobulk 后打分（横向 barplot），后两种逐细胞打分（画在 UMAP 上）。

### 输入

- `-I`：已注释的 Seurat rds（需有 UMAP，用于画单细胞图）
- `-G`：**基因集文本文件，每行一个基因符号**（必填）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| Seurat rds |
| `-G, --geneset` | 无（必填）| 基因集文件，每行一个基因符号 |
| `-AC, --annotation_column` | 无 | 细胞类型列名；脚本用 `Idents(seu) <- <该列>` 后按 Idents 聚合 |
| `-A, --assay` | RNA | 使用的 assay |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_genesetScore.r -I 05_annotation/05.seu_annotated.rds \
    -G geneSet.txt -AC celltype -A RNA -P 09. -O 09_genesetScore
```

基因集文件格式（每行一个符号）：

```
GeneA
GeneB
GeneC
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>GSVA_ssGSEA_score.csv` | 各细胞类型的 GSVA 与 ssGSEA 分数（两列合并）|
| `<前缀>GSVA_score_barplot.png/.pdf` | GSVA 分数横向柱状图：正值红、负值蓝，按分数排序 |
| `<前缀>ssGSEA_score_barplot.png/.pdf` | ssGSEA 分数横向柱状图（统一浅蓝），按分数降序 |
| `<前缀>AddModuleScore_UMAP.png/.pdf` | AddModuleScore 的 UMAP 特征图（灰→红）|
| `<前缀>AUCell_score_UMAP.png/.pdf` | AUCell 分数的 UMAP 特征图（灰→红）|
| `<前缀>seu_genesetScore.rds` | 写入了 `AUCell_score`（以及 AddModuleScore）列的对象 |

示例输出：`Images/genesetScore/`。

### 注意事项

- **细胞类型聚合依赖 `-AC`**：脚本先执行 `Idents(seu) <- annotation_column`，再用 `AggregateExpression`（未指定 `group.by`，即按 Idents 聚合）。因此 `-AC` 必须给出对象里真实存在的列名；不给或写错会导致聚合结果异常。
- **基因集里找不到的基因会被丢弃**：脚本会打印「读入多少个 / 去重后多少个 / 在对象中找到多少个 / 没找到多少个」，并列出缺失的基因名。**基因集与对象物种/基因名版本不一致时缺失会很多**，此时分数不可信。
- **基因集在对象里至少要有一些基因**：全部缺失时 GSVA 会报错。
- **GSVA 与 ssGSEA 的分数不可直接比较**：两者尺度不同，只在同一方法内部横向比较细胞类型。
- **GSVA 分数的正负只表示方向**：正值表示该基因集在排序中偏靠前（相对更活跃），负值表示偏后，不代表「有」或「无」表达。
- **单细胞方法（AddModuleScore / AUCell）需要对象里有 UMAP**，否则出图会报错。
- **AUCell 会覆盖对象里的 `AUCell_score` 列**（若已存在同名列）。
- 运行时会打印 GSVA 的计算进度与基因集匹配信息，属正常提示。


---

## 21 scRNA_DoubleFinder.r（替代方案）—— 用 DoubletFinder 去双细胞

### 功能

**与第 2 章的 `scRNA_dedoublets_step2.r` 做同一件事**（去双细胞 + 可选 DecontX 去环境 RNA），但用**不同的算法**和**不同的处理粒度**。

它是一个**替代方案**，不是主流程的一步。**二选一使用，不要串联**（串联会重复去双细胞）。

### 与主流程 step2 的区别

| | `scRNA_dedoublets_step2.r`（主流程，推荐） | `scRNA_DoubleFinder.r`（本脚本，替代方案） |
|---|---|---|
| 双细胞算法 | **scDblFinder**（Bioconductor）| **DoubletFinder**（CRAN）|
| 处理粒度 | **逐样本**（`samples = -SC`）| **整个合并对象一起** |
| DecontX | 逐样本跑，**逐样本**算 90% 分位阈值 | 整个对象跑，**全局** 90% 分位阈值 |
| 双细胞率 | 可传；不传时 scDblFinder **自动估计** | **必须手填** `-DR`（默认 0.075 是通用经验值，不是你的数据实际的）|
| pK / pN | 不需要 | 需 `paramSweep` 扫 pK + `find.pK` 取最优（这一步较慢）；pN 固定 0.25 |
| 同型双细胞校正 | scDblFinder 内部处理 | 用 `modelHomotypic` 估算后调整期望双细胞数 |
| 输出对象 | `seu_clean.rds`（**只保留 counts 层**，轻量）| `seo_deCont.rds`（保留全部层，含 data/scale/PCA/UMAP/聚类）|
| 汇总记录 | 有 `Dedoublets_summary.txt`（每样本前后细胞数、双细胞统计）| 无，只打印到终端 |
| 缓存跳过 | 有（检测到 `seu_doublets.rds` 就只重画图）| 无 |
| 图尺寸 | 可调 `-W` / `-H` | 写死 |
| 版权年份 | 2026 | 2025 |

### 选择依据

**默认用 step2（scDblFinder）**，理由：
- **逐样本建模**：双细胞率是样本内概念（每个样本的 loading rate 不同）。本脚本把多样本混在一起算期望双细胞数 `nExp = double_rate × 总细胞数`，样本间的 loading 差异被抹平——样本数多、每样本细胞数差异大时更吃亏
- **双细胞率可自动估计**：scDblFinder 在 `dbr = NULL` 时会按数据估计；本脚本的 `-DR` 只是个假设值
- 输出更轻，且带 summary 与缓存

**什么时候用本脚本（DoubleFinder）**：
- 想用**另一种独立算法交叉验证** scDblFinder 的结论——两者都是成熟方法，若某个样本上结论差异很大，值得警惕
- 你更熟悉 DoubletFinder 的参数体系（pK/pN 扫描），需要自己调参

### 输入

- `-I`：QC 之后的 Seurat rds（需含 `orig.ident` 或其它样本信息；本脚本**不按样本拆分**，整体处理）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| 输入 Seurat rds |
| `-D, --decontx` | 关闭 | 是否运行 DecontX 去环境 RNA（`action="store_true"`，**默认关闭**；要跑必须显式加 `-D`）|
| `-PC` | 30 | PCA / 降维 / 聚类用的主成分数 |
| `-R` | 0.1 | 聚类分辨率 |
| `-CA` | Louvain | 聚类算法：`Louvain` 或 `Leiden`（Leiden 需 `-P_env` 指定含 `leidenalg` 的 Python 环境）|
| `-P_env` | 无 | Leiden 用的 Python 环境路径 |
| `-DR` | 0.075 | **期望双细胞率**（双细胞占液滴的比例）|
| `-T` | 10 | 线程数 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 加上 -D 才做 DecontX 去环境 RNA（默认关闭）
./scRNA_DoubleFinder.r -I 01_QC/01.seu_merged.rds -D \
    -PC 30 -R 0.1 -CA Louvain -DR 0.075 -T 20 -P 03. -O 03_deContamination
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>seo_contx.rds` | 跑完 DecontX 的对象（带 `decontX_contamination`、`decontX_class` 列）|
| `<前缀>decontX_umap.png/.pdf` | DecontX 污染程度 FeaturePlot + L-prob/H-prob 分组 DimPlot |
| `<前缀>decontX_violin.png/.pdf` | 按 L-prob/H-prob 分组的 nCount/nFeature/percent.mt 小提琴图 |
| `<前缀>seo_doublets.rds` | 跑完双细胞检测的对象（带 `Double_score`、`Is_Double` 列）|
| `<前缀>doublet_umap.png/.pdf` | 双细胞得分 FeaturePlot + Singlet/Doublet 分组 DimPlot |
| `<前缀>doublet_violin.png/.pdf` | 按 Singlet/Doublet 分组的 nCount/nFeature/percent.mt 小提琴图 |
| `<前缀>seo_deCont.rds` | **最终输出**：过滤掉 Doublet（以及 H-prob）后的对象 |

示例输出：`Images/DoubleFinder/`。

### 注意事项

- **与 step2 二选一，不要串联**：两个脚本都会删双细胞，串联会导致重复过滤、细胞数偏少。
- **`-D` 默认关闭**：与 step2 行为一致——`action="store_true"` 且 `default = FALSE`，要跑 DecontX 必须显式加 `-D`。
- **双细胞率 `-DR` 是假设值**：默认 0.075 是 10x 数据的通用经验值，不是从你的数据估计出来的。样本越多、上样浓度越低，实际双细胞率通常越高。**用 scDblFinder 时可以让它自动估计，是本脚本的主要短板**。
- **不按样本拆分**：所有样本的细胞一起参与建模与过滤，样本间 loading rate 差异被抹平。
- **DecontX 阈值是全局的**：90% 分位在合并后的全部细胞上算，不是每个样本各自算——同一样本内污染程度差异大的数据要注意这点。
- **`find.pK` 会在当前工作目录写一个临时 `pk_plot.png`**（跑完自动删除），不是写进 `-O`；在只读目录里运行会失败。
- **输出对象较重**：`seo_deCont.rds` 保留了 data / scale.data / PCA / UMAP / 聚类等全部层（GB 级），不像 step2 的 `seu_clean.rds` 只留 counts。下一步若要接 step3（归一化 + Harmony），需要留意对象里已有的层。
- **`-D` 会先生成 `seo_contx.rds` 再继续**，即使后续失败也能复用中间结果。
- 运行时会打印每样本过滤前后的细胞数表格，可用于核对。


---

## 22 scRNA_ROC.r（独立脚本）—— 用基因表达做组间判别分析

### 功能

**给一组目标基因（或一个基因集），评估它们能否把两组样本区分开**：逐基因在样本层面算 AUC，画 ROC 曲线，最后给出按 AUC 排序的汇总表。

与库内其它脚本的分工：

| 脚本 | 回答什么 |
|---|---|
| `scRNA_ROC.r` | **判别能力**：某个基因（或基因集）能否把两组分开？AUC 多少？|
| `scRNA_gene_DE.r` | **差异大小**：某基因在两组间表达差多少？差异是否显著？|
| `scRNA_spRNA_cellPro.r` | **构成差异**：某细胞类型的占比在两组间是否不同？|

ROC 是**判别分析**，不回答「差异是否显著」——AUC 高说明两组可分，但不给出 p 值。

### 输入

- `-I`：Seurat rds
- `-SC`：样本身份列（如 `orig.ident`）
- `-G`：**分组信息 CSV**，两列 `sample` 与 `group`，一行一个样本
- `-F`：**基因列表文件**，每行一个基因符号
- `-C`：**要比较的两个组名**（逗号或空格分隔，必须恰好 2 个）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| Seurat rds |
| `-SC, --sample_column` | 无（必填）| 样本身份列名 |
| `-G, --group_info` | 无（必填）| 分组信息 CSV（`sample` + `group` 两列）|
| `-F, --gene_file` | 无（必填）| 基因列表文件，每行一个基因 |
| `-C, --compare` | 无（必填）| 要比较的两个组名 |
| `-GS, --geneset` | FALSE | 改为**基因集打分模式**：用 `AddModuleScore` 把整个基因列表当作一个特征 |
| `-W, --plot_width` | 12 | ROC 图宽（英寸）|
| `-H, --plot_height` | 12 | ROC 图高（英寸）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_ROC.r -I 05_annotation/05.seu_annotated.rds -SC orig.ident \
    -G group.csv -F genelist.txt -C control,treat -W 12 -H 12 -P 09. -O 09_ROC
```

分组信息 CSV 格式：

```
sample,group
S1,control
S2,control
S3,treat
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>AUC_result.csv` | 每个基因（或基因集）的 AUC 与方向，按 AUC 降序；行名即基因名 |
| `<前缀>ROC_coordinates.csv` | 所有 ROC 曲线的坐标（specificity / sensitivity / threshold / FPR / TPR / gene）|
| `<前缀>ROC_<基因名>.png/.pdf` | 每个基因一张 ROC 曲线图，标题含基因名与 AUC |

示例输出：`Images/ROC/`。

`AUC_result.csv` 的 `direction` 列说明哪一组更高：`ctrl_high` 表示第一组更高，`treat_high` 表示第二组更高。

### 注意事项

- **统计单位是样本**：每个样本取该基因的平均表达量作为预测值，因此 AUC 的可信度取决于**每组的样本数**，而不是细胞数。样本数很少时 AUC 波动很大，只能当探索性参考。
- **`-GS` 开关改变的是特征的定义**：不加时每个基因各算一条 ROC；加了 `-GS` 则整个基因列表被 `AddModuleScore` 合成一个分数，只出一条 ROC。两种模式的 `AUC_result.csv` 行含义完全不同（前者行名是基因名，后者行名是 `GeneSet`）。
- **基因必须在对象里**：脚本会打印「找到 / 未找到」的基因清单，未找到的直接剔除。若大部分基因都没找到，说明基因名版本与对象不一致，此时结果不可信。
- **重复基因会被去重**，日志会列出重复的基因名。
- **`-C` 必须是恰好两个组**，否则脚本直接报错停止；组名不在分组文件里也会报错停止。
- **分组文件里的样本必须都能在对象里找到**，否则报错停止（避免静默丢样本）。
- **只分析 `-C` 指定的两组**：`group_info` 里其他组的样本会被排除，不参与计算。
- **AUC 方向由数据自动决定**：`pROC` 会判断哪一组预测值更高并记录在 `direction` 列，不要假设方向一定是「处理组更高」。
- **AUC > 0.7 只是筛选参考**，不是显著性判据；脚本会打印 AUC > 0.7 的基因列表，但真正要下结论请回到 `scRNA_gene_DE.r` 或 `scRNA_pseudobulk_DE.r` 拿差异与 p 值。
- **单个基因的 AUC 对表达量级敏感**：高表达基因即使两组差异不大也可能有较高 AUC，因为 ROC 只看排序不看绝对值。
- **运行时会有 ggplot2 的 `size` → `linewidth` 弃用警告**，属正常提示。


---

## 23 scRNA_cellchat_single.r（独立脚本）—— 单个对象的细胞通讯分析（CellChat）

### 功能

**推断细胞类型之间的通讯信号**：用 CellChat 算法找出「谁在给谁发信号、用什么配受体对、走哪条信号通路」，并给出网络图、通路图与通讯模式。

具体回答三个问题：
1. **哪些细胞类型之间有通讯？**（网络图的连线，粗细 = 通讯强度）
2. **通过什么分子？**（配受体对表 `communication_genes.csv`、信号通路表 `communication_pathways.csv`）
3. **通讯模式有几类？**（`identifyCommunicationPatterns` 的 outgoing pattern，把细胞群按「发出的信号」聚类）

### 输入

- `-I`：已注释的 Seurat rds（需有 `orig.ident` 列——脚本用它构造 `samples`）
- `-AC`：**细胞类型列名**（作为 CellChat 的分组依据）
- `-S`：物种（`human` / `mouse`），决定用哪套配受体数据库与 PPI 网络

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| Seurat rds |
| `-AC, --annotation_column` | 无（必填）| 细胞类型列名 |
| `-S, --species` | 无（必填）| `human` 或 `mouse` |
| `-TP, --type` | triMean | 每个细胞群平均表达的计算方式：`triMean`（默认，出的相互作用更少但更强）/ `truncatedMean`（配合 `-TM`，出的更多）|
| `-TM, --trim` | 0.1 | `truncatedMean` 时两端截去的比例（0–0.25）|
| `-NP, --nPatterns` | 无 | 通讯模式的数量 k，需看 `selectK_outcoming.png` 后确定 |
| `-MC, --min_cell` | 10 | 每个细胞群参与通讯分析的最少细胞数 |
| `-T, --threads` | 24 | 线程数 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 第一步：不带 -NP 跑，生成 selectK 图后看 k 的取值
./scRNA_cellchat_single.r -I 05_annotation/05.seu_annotated.rds \
    -AC celltype -S mouse -TP triMean -TM 0.1 -T 20 -P 10. -O 10_cellChat

# 第二步：带上 -NP 重新跑，生成通讯模式相关图
./scRNA_cellchat_single.r -I 05_annotation/05.seu_annotated.rds \
    -AC celltype -S mouse -TP triMean -TM 0.1 -NP 3 -T 20 -P 10. -O 10_cellChat
```

**必须分两步跑**：脚本在 `-NP` 未给时会先出 `selectK` 图，然后 `stop()` 要求你指定模式数；不指定就跑不下去。

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>communication_genes.csv` | 配受体对层面的通讯表（source / target / ligand / receptor / prob / pval）|
| `<前缀>communication_pathways.csv` | 信号通路层面的通讯表（source / target / pathway_name / prob / pval）|
| `<前缀>cellchat.rds` | CellChat 对象（含网络、中心性、通讯模式），**下次跑同目录会自动复用** |
| `<前缀>circle_overall.png/.pdf` | 总体网络圈图：左 = 通讯数量，右 = 通讯强度 |
| `<前缀>circle_split_count.png/.pdf` | 每个发送方细胞群单独的通讯数量圈图 |
| `<前缀>circle_split_weight.png/.pdf` | 每个发送方细胞群单独的通讯强度圈图 |
| `<前缀>top9_signalPath_chord.png/.pdf` | 概率最高的 9 条信号通路的弦图（3×3 布局）|
| `<前缀><通路名>_LR_distribution.png/.pdf` | 该通路上各配受体对的贡献占比 |
| `<前缀><通路名>_LR_expr.png/.pdf` | 该通路配受体在细胞群中的表达 |
| `<前缀><通路名>_role_network.png/.pdf` | 该通路的信号角色网络（发送者 / 接收者）|
| `<前缀><通路名>_role_scatter.png/.pdf` | 信号角色的 2D 散点（横轴 = 发送角色，纵轴 = 接收角色）|
| `<前缀><LR对名>_circle.png/.pdf` | 贡献最大的 3 个配受体对各自的圈图 |
| `<前缀>allLRs_bubble.png/.pdf` | 全部配受体对的总体气泡图（24×24 英寸）|
| `<前缀>selectK_outcoming.png/.pdf` | 选择模式数 k 的参考图（**第一步生成**）|
| `<前缀>outgoing_patterns_heatmap.png/.pdf` | 通讯模式热图 |
| `<前缀>outgoing_patterns_river.png/.pdf` | 通讯模式河流图 |
| `<前缀>outgoing_patterns_dot.png/.pdf` | 通讯模式点图（各模式的发送/接收特征）|

示例输出：`Images/cellchat/single/`。

### 注意事项

- **必须分两步跑**：不给 `-NP` 时脚本出完 `selectK` 图就 `stop()`，提示你指定模式数。
- **物种必须匹配**：`human` 用 `CellChatDB.human` + `PPI.human`，`mouse` 用 `CellChatDB.mouse` + `PPI.mouse`，传别的值直接报错。
- **`cellchat.rds` 会自动复用**：输出目录里已有该文件时，脚本跳过整个建模流程直接读缓存——**改了 `-TP`/`-TM`/`-MC` 也不会重算**，想重算请先删掉 `cellchat.rds`。
- **`-NP` 决定通讯模式的个数**：k 太小模式粗、太大模式碎；参考 `selectK_outcoming.png` 里 elbow 位置。
- **显著通路不足 9 条时会有 warning**：脚本只画实际找到的通路，格子会空着。
- **`-TP triMean` vs `truncatedMean`**：默认的 `triMean` 只保留三群中至少两群表达的配受体，出的相互作用更少但更强；`truncatedMean` 配合 `-TM` 会保留更多弱信号。
- **`-MC` 是硬门槛**：细胞数少于该值的细胞群不参与通讯分析，通常设 10–20。
- **细胞类型数不宜过多**：CellChat 的网络图与弦图在群数多时可读性会明显下降。
- **输出的 `cellchat.rds` 体积很大**（百 MB 到 GB 级），包含完整 Seurat 对象——因为 `createCellChat` 传了 `meta = seu@meta.data`，但 CellChat 仍会保留部分数据。空间紧张时跑完可以删。
- **单细胞数据做 CellChat 的固有局限**：通讯推断基于「配体和受体是否在各自细胞群表达」的统计推断，**不能证明真实的物理接触**；另外 CellChat 的官方流程更适合 bulk 数据的多个个体重复，单细胞里通常按细胞类型汇总后做，样本层面的重复信息会被丢弃。

---

## 24 scRNA_cellchat_compare.r（独立脚本）—— 两组之间的通讯网络对比

### 功能

**比较两组的细胞通讯网络有什么不同**：把两个 CellChat 对象合并，逐项对比通讯数量、强度、通路组成与差异显著的配受体对。

与单对象脚本的关系：**先跑 `scRNA_cellchat_single.r` 分别得到两组的 `cellchat.rds`，再跑本脚本做对比**。本脚本不重新建模，只做比较。

### 输入

- `-IC`：对照组的 CellChat rds（由 single 脚本生成）
- `-IT`：处理组的 CellChat rds
- `-GN`：两个组名，逗号分隔（如 `control,treatment`），**必须恰好两个**

两个对象的细胞身份必须完全一致，否则脚本报错停止。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-IC, --input_C` | 无（必填）| 对照组 CellChat rds |
| `-IT, --input_T` | 无（必填）| 处理组 CellChat rds |
| `-GN, --group_names` | 无（必填）| 两个组名，逗号分隔 |
| `-TI, --target_idents` | 无 | 只比较这些细胞身份，逗号分隔；不给则全部比较 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_cellchat_compare.r -IC 10_cellChat/10.control_cellchat.rds \
    -IT 10_cellChat/10.treat_cellchat.rds -GN control,treat \
    -P 11. -O 11_cellChat_compare
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>merged_cellchat.rds` | 合并后的 CellChat 对象，**已存在时会自动复用** |
| `<前缀>comparison_barplot.png/.pdf` | 两组的通讯数量与强度对比柱状图 |
| `<前缀>comparison_circleplot.png/.pdf` | 差异网络圈图（左 = 数量，右 = 强度，增粗的边 = 组间差异显著）|
| `<前缀>weight_comparison_heatmap.png/.pdf` | 通讯强度热图 |
| `<前缀>overall_circleplot.png/.pdf` | 2×2 布局：两组各自的数量图与强度图，边宽统一便于对比 |
| `<前缀>role_2d.png/.pdf` | 两组的信号角色 2D 散点并排 |
| `<前缀>pathway_functionalClustering.png/.pdf` | 通路的功能相似性聚类嵌入图 |
| `<前缀>pathway_functionalSimilarity.png/.pdf` | 通路相似性排名热图 |
| `<前缀>weight_DEG_pathways.png/.pdf` | 组间差异显著的通路（堆叠与非堆叠两种）|
| `<前缀>outgoing_pathway_heatmap.png/.pdf` | 两组的 outgoing 通路角色热图并排 |
| `<前缀>LRs_comparison_bubble.png/.pdf` | 配受体对比较气泡图 |
| `<前缀>LRs_comparison_bubble_changed.png/.pdf` | 只显示各组增强的配受体对 |
| `<前缀><通路名>_LRs_expression_comparison.png/.pdf` | 差异最大通路的配受体表达对比小提琴图 |

示例输出：`Images/cellchat/compare/`。

### 注意事项

- **必须先跑 single 脚本**：本脚本的输入是 CellChat 对象，不是 Seurat 对象。
- **两个对象的细胞身份必须完全一致**：脚本会比较 idents 并在不一致时 `stop()`。
- **`-GN` 顺序有意义**：第一个是「对照」、第二个是「处理」，影响 `direction`、图例顺序与「Increased signaling in X」这类标题。
- **`-TI` 只在两个对象都含这些细胞类型时有效**：名字写错会导致该群缺失。
- **`merged_cellchat.rds` 会自动复用**：已存在时直接读缓存。**但缓存里的对象是用当时的 `-TI` 生成的**，改了 `-TI` 想重算请先删掉该文件。
- **`setwd(outpath)` 会改变当前工作目录**：脚本开头会切到输出目录，跑完后当前目录就变了；调用时最好用绝对路径。
- **`rankNet` 的显著性来自 CellChat 的置换检验**，不是常规 p 值；`cutoff.pvalue = 0.05` 是硬编码的。
- **`pathway_functionalSimilarity` / `weight_DEG_pathways` 依赖网络相似性分析**，用的是功能相似性（基于通路共享的配受体），不是基因集富集。
- **`LRs_comparison_bubble` 的 `sources.use = 1`**：只画第一个细胞群作为发送方的比较气泡图；若要换，改脚本里的这个参数。
- **输出同样包含大体积 rds**（`merged_cellchat.rds`），注意磁盘空间。


---

## 25 scRNA_spRNA_Monocle3.r（独立脚本）—— 拟时序与轨迹分析（Monocle3）

### 功能

**给细胞排出一条「发育/转变」的先后顺序，并找出沿这条顺序变化的基因**：用 Monocle3 学习轨迹图、指定根细胞、计算拟时序（pseudotime），再做空间自相关检验与基因模块分析。

具体回答三个问题：
1. **细胞之间有没有内在的转变顺序？**（轨迹图的分支结构）
2. **哪些基因沿这条顺序变化？**（Moran's I 空间自相关检验找出在拟时序上有趋势的基因）
3. **这些基因能分成几个功能模块？**（基因模块分析 + 模块热图）

与其它脚本的分工：

| 脚本 | 回答什么 |
|---|---|
| `scRNA_spRNA_Monocle3.r` | **细胞间的先后顺序**：谁先谁后、沿轨迹什么在变 |
| `scRNA_clustering_step4.r` | **细胞分成几群**（无方向性）|
| `scRNA_subcelltype_clustering.r` | **细胞类型内部分亚型**（无方向性）|

拟时序分析是**有方向**的：根细胞的选择决定了时间轴的起点，选错根会得到完全不同的结论。

### 输入

- `-I`：Seurat rds。**建议输入已跑过 UMAP 的对象**（若没有，脚本会重新做归一化 + PCA + UMAP 并另存一份）
- `-DT`：`SC`（单细胞）或 `SP`（空间转录组）。选 `SP` 时 metadata 必须含 `x`、`y` 两列（空间坐标）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| Seurat rds |
| `-DT, --data_type` | 无（必填）| `SC` 或 `SP` |
| `-C, --cluster_column` | 无 | 指定根细胞的 metadata 列名；**必须与 `-G` 同时给** |
| `-G, --group` | 无 | `-C` 列中代表根细胞群体的取值；**必须与 `-C` 同时给** |
| `-R1, --resolution1` | 1e-5 | Monocle3 聚类分辨率（决定轨迹图上有多少个簇/分支）|
| `-R2, --resolution2` | 0.5 | 基因模块分析的分辨率 |
| `-T, --threads` | 10 | 线程数（Moran's I 与模块分析）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 单细胞
./scRNA_spRNA_Monocle3.r -I 05_annotation/05.seu_annotated.rds -DT SC \
    -C celltype -G <root_celltype> -R1 1e-5 -R2 0.5 -T 20 -P 11. -O 11_monocle3

# 空间转录组（metadata 需含 x、y 列）
./scRNA_spRNA_Monocle3.r -I spatial.seo.rds -DT SP \
    -C celltype -G <root_celltype> -R1 2e-4 -R2 0.5 -T 40 -P 11. -O 11_monocle3
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>seo_nor_umap.rds` | 仅当输入对象缺 UMAP 时才生成（重新归一化 + PCA + UMAP 后的对象）|
| `<前缀>seo_cds.rds` | Monocle3 cell_data_set 对象（含拟时序、轨迹图、基因模块）|
| `<前缀>umap_cluster.png/.pdf` | UMAP 上的轨迹图：按 Monocle3 簇着色，画出轨迹图的分支结构 |
| `<前缀>umap_pseudotime.png/.pdf` | UMAP 上的拟时序着色（viridis 色阶）|
| `<前缀>moran_result.csv` | 每个基因的 Moran's I 统计量（morans_I / p_value / q_value / gene）|
| `<前缀>moran_DEGs_pseudotime.png/.pdf` | 沿拟时序变化最显著的 12 个基因的表达趋势曲线（3 列网格）|
| `<前缀>gene_modules.csv` | 基因模块归属表（基因 → 模块编号）|
| `<前缀>module_heatmap.png/.pdf` | 模块 × 细胞群的热图（按细胞类型聚合，列 z-score，ward.D2 聚类）|
| `<前缀>spatial_pseudotime.png/.pdf` | 空间坐标上的拟时序分布散点图（**仅 `-DT SP` 时生成**）|

示例输出：`Images/Monocle3/`。

### 注意事项

- **根细胞的选择决定一切**：`-C` + `-G` 指定的是拟时序的起点。给错根会让整条时间轴反过来或指向错误方向，**结论会完全改变**。建议先用 `umap_cluster.png` 看懂轨迹结构再选根。**不给 `-C`/`-G` 时脚本会用 Monocle3 的 1 号簇作根**，并给出 warning。
- **`-C` 与 `-G` 必须成对出现**：只给一个会直接报错退出（脚本用 `is.null(cluster_column) != is.null(group)` 判断）。
- **`-R1` 控制轨迹的分支数**：默认 1e-5 对应较粗的划分；数据复杂时可调大（历史项目里用过 2e-4）。调大后簇更多、轨迹图更碎。
- **模块分析需要至少 30 个合格基因**：筛选条件是 `p_value < 0.05 & q_value < 0.1 & morans_I > 0.1`，不足 30 个就跳过模块分析（只打日志，不报错）。这在细胞数少或数据噪声大时很常见。
- **Moran's I 的两个阈值是写死的**（p<0.05、q<0.1），`moran_DEGs_pseudotime` 图固定画 12 个基因——这些在脚本里没有参数可调。
- **`-DT SP` 要求 metadata 有 `x`、`y` 列**，否则脚本在出图前直接 `quit(status = 1)`。
- **输出对象极大**：`seo_cds.rds` 保留了完整表达矩阵，实测单样本可达 **2.8 GB**。磁盘紧张时跑完及时删除。
- **单细胞数据做 Monocle3 的固有局限**：拟时序是「表达谱相似性」推断出的顺序，**不等于真实的时间或谱系关系**；细胞群体的比例变化（组成差异）不会被 Monocle3 建模——那是 scCODA 的领域。
- **空间数据用 `x`、`y` 列做散点**，不做真正的空间邻域计算（与 `spatialDE` / `SPARK` 不同）。
- 运行时若输入对象缺 UMAP 会先重新降维，耗时较长——建议直接喂已跑好 UMAP 的对象。


---

## 26 scRNA_velocityPrepare.r（配套脚本）—— 为 RNA velocity 准备输入对象

### 功能

**把「含 spliced/unspliced 的完整对象」与「实际要分析的对象」对齐，重建成一个同时带三个 assay 的 velocity 输入对象，并转成 h5ad**。

这是 `scRNA_velocity.py` 的**前置步骤**：RNA velocity 需要同时知道每个细胞的成熟 mRNA（spliced）和未剪接 mRNA（unspliced）量，而通常这两层只存在于原始计数对象里，而分析用的是注释后的子集对象。本脚本负责把两者拼起来。

具体做四件事：
1. 用参考对象的**细胞名与基因名**做索引，从完整对象里取出对应的三层 counts（RNA / spliced / unspliced）
2. 建一个新的 Seurat 对象：**metadata 与降维结果全部继承参考对象**，counts 换成完整对象的三层
3. 转成 AnnData：`X` = RNA counts，`layers['spliced']` / `layers['unspliced']` = 对应层，`obsm` 里放 UMAP / PCA / Harmony 坐标
4. 同时保存 rds 与 h5ad 两个格式

### 输入

- `-R`：**参考 Seurat 对象**（分析用的子集，含注释、降维结果）——决定最终保留哪些细胞、哪些基因
- `-F`：**完整 Seurat 对象**（含 `RNA` / `spliced` / `unspliced` 三个 assay 的 counts）

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-R, --ref_seurat` | 无（必填）| 参考对象 rds（分析子集）|
| `-F, --full_seurat` | 无（必填）| 完整对象 rds（含 spliced / unspliced）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
./scRNA_velocityPrepare.r \
    -R 05_annotation/05.seu_annotated.rds \
    -F 01_counts/full_seurat_with_velocity.rds \
    -P 12. -O 12_velocity
```

### 输出

| 文件 | 说明 |
|---|---|
| `<前缀>scRNA_velocity_input.rds` | Seurat 格式的 velocity 输入对象（RNA + spliced + unspliced 三个 assay，继承参考对象的 metadata 与降维）|
| `<前缀>scRNA_velocity_input.h5ad` | AnnData 格式（`X` = RNA counts，`layers['spliced']` / `layers['unspliced']`，`obsm` 含 `X_umap` / `X_pca` / `X_harmony`）|

### 注意事项

- **两个对象的细胞名与基因名必须完全一致**：脚本用 `match()` 对齐并 `stopifnot()` 强制检查，任何一个对不上就直接报错退出（这是好事，避免静默错配）。
- **细胞/基因的取舍完全由参考对象决定**：完整对象里多出来的细胞和基因会被丢掉。
- **`obsm` 里只写参考对象已有的降维**：有 UMAP 就写 `X_umap`，有 PCA 写 `X_pca`，有 Harmony 写 `X_harmony`；没有就不写。**下游 `scRNA_velocity.py` 的邻居图用的是 `X_harmony`**，所以参考对象里最好有 Harmony 结果。
- **`scRNA_velocity.py` 的 subset 模式需要 `X_pca_harmony`**：该脚本在指定 `-TC` 子集时会自己跑 PCA + harmony（`X_pca_harmony`）；不指定 `-TC` 时它直接用 `X_harmony`——此时参考对象**必须**已有 Harmony 坐标，否则会报错。
- **转换 AnnData 时会稠密化**：`as.matrix()` 把三层稀疏矩阵转成 dense，细胞基因数多时内存占用会暴涨。这是这一步的主要瓶颈。
- **需要 `anndataR`**：本机的 `anndataR` 支持 `layers` 参数，新版写法 `adata$layers[["spliced"]] <- ...` 可用。
- 输出的 h5ad 体积约等于三层 counts 的总和（RNA + spliced + unspliced）。

---

## 27 scRNA_velocity.py（独立脚本）—— RNA velocity 分析（scVelo）

### 功能

**推断每个细胞的状态变化方向**：用 scVelo 计算 RNA velocity——每个细胞不只是一个点，而是有「往哪个方向走」的速度矢量。把它叠在 UMAP 上就得到流线图（stream plot），回答「细胞群体的转变方向是什么」。

两个模型可选：
- **stochastic**：更快更稳，只用表达量做局部邻域的马尔可夫链
- **dynamical**（默认）：更精细，先从稳态模型恢复每个基因的转录动力学参数再算速度，能给出 **latent time**（比 Monocle3 的伪时序更贴近真实时间）与 **velocity genes**（哪些基因驱动了转变）

### 输入

- `-I`：**h5ad 文件**（由 `scRNA_velocityPrepare.r` 生成），必须含 `layers['spliced']`、`layers['unspliced']`、`obsm['X_umap']`、`obsm['X_harmony']`
- `-AC`：`adata.obs` 里的细胞注释列名

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input-file` | 无（必填）| velocity 输入 h5ad |
| `-AC, --annotation-column` | 无（必填）| 注释列名 |
| `-M, --mode` | dynamical | 模型：`stochastic` 或 `dynamical` |
| `-TC, --target-cell` | 无 | 只分析这些细胞类型（逗号分隔）；不给则用全部细胞 |
| `-T, --threads` | 4 | 线程数 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 先用 prepare 脚本准备输入，再跑分析
./scRNA_velocity.py -I 12_velocity/12.scRNA_velocity_input.h5ad \
    -AC celltype -M dynamical -TC 'CellTypeA,CellTypeB' -T 40 -P 12. -O 12_velocity/result
```

### 输出

输出目录 `-O` 下：

| 文件 | 说明 |
|---|---|
| `<前缀>dynamical_stream.png/.svg` | dynamical 模型的 velocity 流线图（仅 `-M dynamical`）|
| `<前缀>dynamical_genes.png/.svg` | velocity genes 的表达与速度图（仅 `-M dynamical`）|
| `<前缀>latent_time_umap.png/.svg` | latent time 在 UMAP 上的分布（仅 `-M dynamical`）|
| `<前缀>velocity_pseudotime_umap.png/.svg` | velocity pseudotime 在 UMAP 上的分布（仅 `-M dynamical`）|
| `<前缀>stochastic_stream_pseudotime.png/.svg` | stochastic 模型的流线图 + velocity pseudotime（仅 `-M stochastic`）|
| `<前缀><mode>_velocity_result.h5ad` | 分析后的 AnnData（含 velocity、latent_time、velocity_genes 等结果）|

示例输出：`Images/velocity/`。

### 注意事项

- **必须先跑 `scRNA_velocityPrepare.r`**：本脚本不处理原始计数，只消费准备好的 h5ad。
- **两种模式的图完全不同**：`stochastic` 只出一张 `stochastic_stream_pseudotime.png`；`dynamical` 出三张（stream / genes / latent_time）+ 一张 velocity_pseudotime。**文件名里带模式名，产物不会混淆**。
- **`-TC` 不给 vs 给，处理流程不同**：不给时直接用 `X_harmony` 建邻居图（要求输入 h5ad 里有 Harmony 坐标）；给了时脚本会自己跑 `filter_and_normalize` → PCA → `harmony_integrate` → `neighbors` → `moments` → UMAP，此时用 `X_pca_harmony`。
- **subset 模式下会重算降维**：指定 `-TC` 时 UMAP 会被重算，所以**子集的 UMAP 与原对象的 UMAP 不可直接比较**（坐标轴变了）。
- **dynamical 模式较慢**：`recover_dynamics` 逐基因拟合转录动力学，是整个流程最耗时的一步；细胞多时建议先用 `-TC` 缩小范围。
- **`scv.pl.velocity` 的 `save` 参数**：代码传的是带 `.png` / `.svg` 后缀的完整文件名，而 scVelo 0.3.x 的 `save` 通常按「路径前缀」处理，可能产出 `xxx.png.png` 之类的重名文件——出图后请确认实际文件名。
- **velocity genes 的筛选是两道交集**：`fit_likelihood > 0.1` 的高质量基因 ∩ 各组 rank velocity 的前 3 名，可能为空。
- **输出是 SVG 而非 PDF**：这是少数几个产出 SVG 的脚本；库里 Images 的示例图按惯例只放 PNG。
- **scVelo 对 spliced/unspliced 的质量敏感**：如果原始计数里 intronic 比例不对，velocity 方向会失真。建议先用 `velocyto` / `aCellBender` 等工具重新生成这两层。
- **velocity ≠ 真实的谱系关系**：它推断的是「当前表达状态正朝哪个方向变化」，需要结合生物学解释，不能单独作为谱系结论。


---

## 28 pySCENIC 转录因子调控网络（三脚本一套流程）

本章包含三个脚本，构成一条完整流程，必须按顺序执行：

| 脚本 | 运行环境 | 角色 |
|---|---|---|
| `scRNA_pyscenic_prepare.r` | R | 把 Seurat 对象压成 pySCENIC 能吃的 loom |
| `scRNA_pyscenic.py` | Python | 跑 pySCENIC 的三步核心算法（grn / ctx / aucell）|
| `scRNA_pyscenic_RSS.r` | R | 算 RSS、出图、把调控活性写回 Seurat 对象 |

### 功能

**找出「哪个转录因子在驱动哪一类细胞的转录状态」。**

常规单细胞分析能告诉你「这一类细胞高表达哪些基因」，但回答不了「**是谁在启动这批基因**」。pySCENIC 补的正是这一环：它用**基因组上的 motif 证据**推断转录因子（TF）到靶基因的调控关系，再给每个细胞算每个 TF 的调控活性。

具体回答三个问题：
1. **哪些 TF 调控哪些基因？**（regulon：TF + 有 motif 证据的靶基因）
2. **每个细胞里哪个 TF 在活跃？**（AUC：每细胞 × 每 regulon 的活性分）
3. **哪个 TF 是某类细胞的标志？**（RSS：regulon 特异性分数）

与库里其它脚本的分工：

| 脚本 | 回答什么 |
|---|---|
| `scRNA_pseudobulk_DE.r` | 哪些**基因**在两组间变了 |
| `scRNA_spRNA_enrichment.r` | 这些基因落在哪些**通路**上 |
| `scRNA_genesetScore.r` / `scRNA_GSVA_enrichment.r` | 某个**已知基因集**的活性高低 |
| **pySCENIC（本章）** | 这些基因是**被谁调控**的——从「现象」推到「上游调控者」|

### 原理与核心概念

这一节是理解整条流程的前提。概念按从上到下的依赖顺序排列。

#### 1. 基因是怎么被开/关的

DNA 上每个基因前面有一段**启动子（promoter）**，相当于开关面板。有一类蛋白叫**转录因子（TF, Transcription Factor）**，它结合到启动子上决定这个基因开多大。所以调控关系是：

```
TF  →  结合到靶基因的启动子  →  控制靶基因的表达量
```

#### 2. 什么是 motif

**motif = 一段 6–15 个碱基的短 DNA 序列，是某个 TF 专门识别的「门牌号」。**

- motif 写在 **DNA** 上，不在蛋白上
- 一个 TF 对应一个（或几个）motif，不同 TF 的 motif 不同
- 这是「谁调控谁」的**物理依据**

为什么这点关键：**共表达 ≠ 直接调控**。基因 A 和 B 表达同步，可能只是都被第三个因素控制，A 并不调控 B。**要证明 A 真的调控 B，必须证明 B 的启动子上有 A 的 motif。**

pySCENIC 用到的两个参考文件就是这个意思：

| 文件 | 存的是什么 | 有多大 |
|---|---|---|
| `*.genes_vs_motifs.rankings.feather` | **每个基因的启动子（TSS ±10kb）里，所有 motif 的排名** ← 「每个基因有哪些门牌号」| 约 1 GB |
| `motifs-v9-nr.*.tbl` | **motif → TF 的对应表** ← 「门牌号属于谁」| 约 100 MB |
| `*_tfs.txt` | **该物种的 TF 名单** | 约 10 KB |

#### 3. 什么是 regulon

**regulon = 一个 TF + 它直接调控的那批靶基因。**

构建过程：
1. GRNBoost2 找出「和 TF X 表达高度相关的基因」→ 候选列表（**只有共表达证据**）
2. 逐个查这些候选基因的启动子：**有没有 TF X 的 motif？**
3. 有的留下，没有的扔掉 → 剩下的这批就构成一个 regulon

命名规则：`TF名(+)` 表示激活型（TF 高则靶基因高），`TF名(-)` 表示抑制型。实际跑出来**绝大多数是 `(+)`**，因为抑制型需要额外的负相关判据，默认参数下很少出现。

#### 4. 什么是 AUC（AUCell 打分）

**AUC = 某个 regulon 在【某一个细胞】里有多活跃**，范围 0–1。

算法：
1. 把这个细胞里所有基因**按表达量从高到低排序**
2. 看这个 regulon 的靶基因**排在前面的占多少**
3. 算这条恢复曲线的下面积（AUC）

直观理解：**如果某个 regulon 的靶基因普遍排在这个细胞表达谱的最前面，说明这个 regulon 在这个细胞里很活跃。**

结果是一张矩阵：**细胞 × regulon**。

- **不要和 `scRNA_ROC.r` 的 AUC 混淆**——那是 ROC 曲线下面积、衡量两组判别能力；AUCell 的 AUC 是基因排序恢复曲线的面积、衡量单细胞里的调控活性。两者同名但毫无关系。

#### 5. 什么是 RSS

**RSS（Regulon Specificity Score）= 某个 regulon 有多「专属」于某一类细胞**，范围 0–1。

**RSS 是从 AUC 算出来的**，算法是 **JSD（Jensen-Shannon 散度）**，不是 z-score：

```
RSS = 1 - sqrt(JSD(regulon 的活性分布, 该细胞类型的分布))
```

拆开说：
1. 把某个 regulon 的 AUC **按细胞归一化** → 得到「这个 regulon 的活性在各细胞间的分布」
2. 该细胞类型的分布 = 「哪些细胞属于这一类」（是 = 1，否 = 0）
3. 算两个分布的 **JSD**——越小说明两个分布越像
4. `1 - sqrt(JSD)` → 越像则 RSS 越接近 1

含义：**RSS 接近 1 说明这个 regulon 只在这一类细胞里活跃；接近 0 说明它在所有细胞里活性都差不多。**

**RSS 数值天然偏小**，这是 JSD 的性质决定的，不是数据差。以数百个 regulon × 十余类细胞的常见规模为例：RSS 中位数往往在 0.05 上下，最高也只有 0.3 左右，`RSS > 0.1` 的格子通常只占两三成。所以读 RSS 要**横向比较**（同一列里谁高谁低），不能把 0.3 当成「只有 30% 特异性」。

#### 6. 概念关系一览

| 概念 | 粒度 | 回答 | 谁产生的 |
|---|---|---|---|
| motif | DNA 序列 | 这个 TF 认什么序列 | CisTarget 数据库 |
| regulon | TF + 靶基因集合 | 谁调控哪些基因 | `pyscenic ctx` |
| AUC | **单个细胞** × 单个 regulon | 这个细胞里这个 TF 活跃吗 | `pyscenic aucell` |
| RSS | **一个细胞类型** × 单个 regulon | 这个 TF 是这类细胞的标志吗 | `calcRSS`（脚本 3）|

一句话：**motif 是依据，regulon 是关系，AUC 是原料，RSS 是提炼后的结论。**

### 三个脚本的分工与数据流

```
Seurat 对象（表达矩阵）
    │
    │  [脚本 1] scRNA_pyscenic_prepare.r
    │  过滤低检出基因 → 可选抽细胞 → 写 loom
    ▼
expr_norm.loom
    │
    │  [脚本 2] scRNA_pyscenic.py    ← 内部调 3 次 pyscenic CLI
    │
    ├─ ① pyscenic grn     （共表达推断）      → adj.csv
    ├─ ② pyscenic ctx     （motif 富集）      → regulons.csv
    └─ ③ pyscenic aucell  （逐细胞活性打分）  → SCENIC.loom
    │
    │  [脚本 3] scRNA_pyscenic_RSS.r
    │  算 RSS → 出图 → 活性写回对象
    ▼
rss.csv / 5 张图 / seu_SCENIC.rds
```

**为什么要拆成三个脚本而不是一个**：这条流程跨**两个环境**（R 与 Python），而 shebang 只能指向一个解释器。拆开之后每个脚本的 shebang 指向自己需要的环境，用户不需要手动 `conda activate` 或切环境：

```
scRNA_pyscenic_prepare.r   #!/…/miniforge3/envs/old_base/bin/Rscript
scRNA_pyscenic.py          #!/…/miniforge3/envs/pySCENIC/bin/python
scRNA_pyscenic_RSS.r       #!/…/miniforge3/envs/old_base/bin/Rscript
```

另外三步的**耗时差极大**（grn 十几分钟 → ctx 可能几小时 → aucell 几分钟），分开跑可以在 ctx 失败后只重跑 ctx，不用从头再来。

### 每一步在做什么，为什么要这么做

#### 脚本 1：准备输入

**① 抽表达矩阵**——用 `-L data`（LogNormalize 之后）而不是 counts。GRNBoost2 基于表达量相关性的回归，用归一化后的值可以避免测序深度带来的伪相关。

**② 过滤基因（`-MG`，默认 50）**——只保留至少在 50 个细胞里检出的基因。一个只在十几个细胞里检出的基因，它的「表达量」主要是 dropout 噪声，算出来的关联是假的。

**③ 抽细胞（`-NC`，默认 `0` = 全部保留）**——GRNBoost2 是逐基因跑梯度提升，成本 ≈ 基因数 × 细胞数。全量细胞会到小时级，且细胞数过万后结果基本稳定，抽样边际收益很小。
⚠️ **抽样会漏掉稀有细胞类型**，如果某类细胞只有几十个，抽 2000 后可能只剩个位数，RSS 就不可靠——这种数据应设 `-NC 0` 用全部细胞。

**④ 写 loom**——loom 是 pySCENIC 的输入格式，必须带 `Gene`（行）与 `CellID`（列）两个属性。

#### 脚本 2：三步核心算法

**⑤ `pyscenic grn`——共表达推断**
对每一个 TF，用梯度提升回归（GRNBoost2）拟合「该 TF 的表达 → 各靶基因的表达」，用树模型的特征重要性作为 TF→靶基因的关联强度。**这一步只看共表达，不看 DNA 序列。** 输出 `TF, target, importance` 三列。

**⑥ `pyscenic ctx`——motif 富集，把共表达升级成直接调控**
对每个 TF 的候选靶基因，查它们的启动子上**这个 TF 的 motif 是否显著富集**，只保留富集的靶基因。这一步是 pySCENIC 的核心，也是最慢的一步。
默认开启的 `--mask-dropouts` 对应 CLI 的 `--mask_dropouts`（关闭用 `--no-mask-dropouts`）：把 dropout（表达量为 0）当作**缺失**而非真实的 0。单细胞里的 0 往往是没测到，当成 0 会低估相关性。

**⑦ `pyscenic aucell`——逐细胞活性打分**
用 AUCell 给每个细胞的每个 regulon 打分（算法见上文「什么是 AUC」）。

#### 脚本 3：特异性分析与出图

**⑧ 读 loom**——**必须用 `SCopeLoomR::open_loom`（hdf5r 直读），不能用 Python 的 `loompy`**。`pyscenic aucell` 写出的 loom 不符合 loom v3 规范（多列属性的 dtype 不合法），`loompy` 会报 `Row attribute 'Regulons' dtype ... is not allowed`。

**⑨ 对齐细胞**——loom 里只有抽样后的细胞，Seurat 对象有全部细胞，取交集。

**⑩ 算 RSS**——见上文「什么是 RSS」。

**⑪ 出图**——RSS 总图、逐细胞类型排名图、regulon 活性小提琴图 / Feature 图 / 点图。三张 `regulon_*` 图画的都是 **AUC 活性分**（不是靶基因表达量），具体读法见下文「结果图怎么读」。

**⑫ 把活性写回对象**——把「细胞 × regulon」的 AUC 矩阵作为 metadata 列加进 Seurat 子集并保存 rds。这一步很重要：否则跑完只有图，想换个 TF 再画一次就得重跑整个 pySCENIC（几小时）。

### 结果图怎么读

本章产出 5 张图，它们是**同一份数据的不同视角**，配合起来读：

| 图 | 作用 | 怎么读 |
|---|---|---|
| `RSS_plot` | **找标志物** | 大点 + 深绿（绝对值 + 相对排名）|
| `RSS_oneSet` | **看细节** | 每类细胞 top5 TF 排名 + 有几个核心 TF |
| `regulon_dotplot` | **验证** | 这些 regulon 的活性是否真的集中在该细胞类型（颜色 = 平均 AUC，大小 = 活跃细胞比例）|
| `regulon_violin` | **验证** | 活性分布是否真的集中在一类细胞 |
| `regulon_feature` | **看连续性** | 是一群细胞点亮还是沿轨迹渐变 |

**唯一需要交叉验证的**：**RSS 是相对量**（分布相似度，回答「在这类细胞里是不是排第一」），**AUC 是绝对量**（回答「这类细胞里到底有多活跃」）。RSS 高只说明相对集中，不代表绝对活性高——所以要用 `regulon_dotplot` / `regulon_violin` 看 AUC 的绝对水平来确认。

**另一个必须记住的点**：本章三张 `regulon_*` 图画的全是 **regulon 的 AUCell 活性分（AUC）**，**不是靶基因的表达量**。AUC 是由「靶基因在该细胞表达谱里的排序位置」算出来的聚合量——靶基因表达量都不高、但都挤在排序前列时，AUC 依然可以很高。

#### `RSS_plot`：点的两个编码别搞混

这是最容易读错的一张图——**点的大小和颜色编码的是两个不同的量**：

| 视觉元素 | 编码 | 含义 | 陷阱 |
|---|---|---|---|
| 点的**面积** | **RSS 绝对值** | 这个 regulon 在这类细胞里到底有多专属 | — |
| 点的**颜色** | **Z 值**（= `scale(rss)`，**按细胞类型列标准化**）| 这个 regulon **在这一类细胞内部**、相对其它所有 regulon 有多突出 | 高 Z 但低 RSS = 其实是「矮子里拔将军」|

源码依据（`SCENIC::plotRSS`）：
```r
rssNorm <- scale(rss)          # scale() 默认按【列】= 细胞类型 标准化
rssNorm[rssNorm < 0] <- 0
if (!is.null(zThreshold)) rssSubset[rssSubset < zThreshold] <- 0
dotHeatmap(rss.df, var.size = "RSS", var.col = "Z", ...)
```
实测确认 `scale(rss)` 的**列均值 = 0、列标准差 = 1**（行均值不是 0），所以是**按细胞类型**标准化。

颜色阶梯（脚本设为 `col.low="#cfd9ed"` → `col.mid="#66CC66"` → `col.high="#FFCC33"`）：淡蓝 = 这一列里平平无奇，绿 = 中等，橙黄 = 这一列里最突出。

**找标志 TF 的正确姿势**：找**又大又深绿**的点。五种典型情况（数值为示意）：

```
regulon_A   细胞类型_1   RSS=0.36  Z=4.0    ← 又大又深绿 → 真正的标志候选
regulon_B   细胞类型_2   RSS=0.21  Z=8.6    ← 中等大小 + 最黄
regulon_C   细胞类型_3   RSS=0.29  Z=6.3
regulon_D   细胞类型_4   RSS=0.02  Z=10.7   ← Z 最高但 RSS 极小 → 不要据此下结论
regulon_E   细胞类型_5   RSS=0.04  Z=10.5   ← 同上
```

⚠️ **`-ZW`（zThreshold，默认 3）会把 Z < 3 的格子置 0**，那些格子只显示很小的淡色点——那不代表「没有调控」，只是没进入可视范围。实际数据里**绝大多数格子都低于这个阈值**（常常超过九成），所以图上看不到颜色是常态。

#### `regulon_dotplot`：两个维度分别看什么

行 = 细胞类型，列 = regulon。**两个维度编码两个不同的问题**：

| 视觉元素 | 编码 | 含义 |
|---|---|---|
| 点的**颜色** | **平均 AUC**（原始值，未标准化）| 这个 regulon 在这类细胞里**到底有多活跃**（绝对水平）|
| 点的**大小** | **AUC 超过该 regulon 全局 75 分位的细胞比例** | 这个 regulon 在这类细胞里**有多普遍地活跃**（一致性）|

**为什么不直接用 Seurat `DotPlot` 默认的「AUC > 0 的细胞比例」**：AUCell 的 AUC 在几乎每个细胞里都 > 0，实测该比例中位数 **0.990**、**一半的 regulon 是 1.000**，点大小会全部顶到最大，失去区分度。改用「超过全局 Q75 的比例」后拉开到 **0.00–1.00**。

**「超过全局 Q75 的比例」怎么理解**（容易误解的一点）：

1. 先算这个 regulon 在**全部细胞**上 AUC 的 75 分位 → 得到**一个固定门槛**（每个 regulon 各自一个数，例如 0.042）
2. 再**在每个细胞类型内部**统计「跨过这个门槛的细胞占多少」
3. 所以**全局范围内这个比例恒为 25%**（按 75 分位的定义），但**拆到每个细胞类型就不一样了**

同一个阈值拿到各类细胞里分别算，结果会差很多（下面是示意量级）：
```
细胞类型_1   1.000   ← 全部细胞都跨过门槛
细胞类型_2   0.56
细胞类型_3   0.58
细胞类型_4   0.04    ← 只有 4%
细胞类型_5   0.02
```
差 60 倍，这就是图上的区分度来源。

**读法：`0.25` 是中性线。**
- **明显 > 0.25** = 这个 regulon 在这类细胞里富集
- **≈ 0.25** = 和全体细胞平均水平一样，没有富集
- **接近 0** = 这类细胞里这个 regulon 基本关闭

**找标志物**：颜色偏黄（平均 AUC 高）**且**点大（比例高）。只有颜色黄但点很小，说明少数细胞很高、整体并不普遍。

#### 其余三张

- **`RSS_oneSet`**：每个细胞类型一张散点排序曲线，x 轴 = regulon 按 RSS 排名的名次，y 轴 = RSS 值，top5 带标签框。**曲线左侧的「悬崖」有多陡 = 这一类细胞有几个核心 regulon**；曲线平缓说明特异性分散在很多 regulon 上。
- **`regulon_violin`**：按细胞类型分组的 **AUC** 分布小提琴图，用来确认这个 regulon 的活性是不是**真的只集中在那一类细胞**（一个小提琴高、其余平 = 干净）。
- **`regulon_feature`**：UMAP 上按 **AUC** 着色。**一群细胞突然点亮 = 对应离散细胞类型；沿轨迹渐变 = 连续的分化过程**。

#### 一句话读图顺序

```
1. RSS_plot         → 哪些 TF 是哪些细胞类型的标志（又大又深绿）
2. RSS_oneSet       → 每类细胞有几个核心 TF、排名前几的是谁
3. regulon_dotplot  → 验证：AUC 绝对水平（颜色）与活跃普遍程度（大小），0.25 为中性线
4. regulon_violin   → 验证：活性分布真的集中
5. regulon_feature  → 空间/连续性：是一群细胞点亮还是沿轨迹渐变
```


### 输入

| 脚本 | 需要什么 |
|---|---|
| `scRNA_pyscenic_prepare.r` | 一个 Seurat rds（建议是已注释的、跑过降维的对象）|
| `scRNA_pyscenic.py` | 脚本 1 产出的 loom + **物种对应的三个参考文件** |
| `scRNA_pyscenic_RSS.r` | 脚本 1 用的同一个 Seurat rds + 脚本 2 产出的 `SCENIC.loom` |

**参考文件（`-RD` 指定的目录下，按 `-S` 物种自动选）**：

小鼠（`-S mouse`）：
- `mm_mgi_tfs.txt`
- `mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.mgi-m0.001-o0.0.tbl`

人（`-S human`）：
- `hs_hgnc_tfs.txt`
- `hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather`
- `motifs-v9-nr.hgnc-m0.001-o0.0.tbl`

三个文件可从 CisTarget 数据库与 pySCENIC 仓库下载，加起来约 1.2 GB。

### 参数

**`scRNA_pyscenic_prepare.r`**

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| Seurat rds |
| `-A, --assay` | RNA | 表达矩阵所在 assay |
| `-L, --layer` | data | 用哪一层：`data`（推荐）或 `counts`；`scale.data` 是稠密矩阵，不支持 |
| `-MG, --min_cells` | 50 | 至少多少个细胞检出才保留该基因 |
| `-NC, --n_cells` | **0（全部细胞）** | 抽样细胞数；默认保留全部细胞。细胞数很多、只求速度时可设 2000 等值抽样 |
| `-SD, --seed` | 123 | 抽样随机种子 |
| `--title` | pySCENIC input | 写进 loom 的标题（无短参数，避免与库内 `-T`=线程的约定混淆）|
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

**`scRNA_pyscenic.py`**

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input-loom` | 无（必填）| 脚本 1 产出的 loom |
| `-S, --species` | 无（必填）| `mouse` 或 `human`，决定用哪套参考文件 |
| `-RD, --ref-dir` | ref_data | 三个参考文件所在目录 |
| `-GM, --grn-method` | grnboost2 | `grnboost2`（快）或 `genie3`（慢但另一套算法）|
| `-ST, --steps` | all | 要跑哪些步骤：`grn,ctx,aucell` 任意组合，或 `all` |
| `-T, --threads` | 12 | 传给 pyscenic 的 worker 数；共享服务器上建议不要拉满 |
| `-SE, --seed` | 777 | grn 步骤的随机种子 |
| `--mask-dropouts` / `--no-mask-dropouts` | 开 | 控制是否给 ctx 传 `--mask_dropouts`；默认开，用 `--no-mask-dropouts` 关闭（不开无短参数，因为「默认开」的开关加短参数容易误读）|
| `-F, --force` | 关 | 输出已存在时也强制重跑 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

**`scRNA_pyscenic_RSS.r`**

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input_rds` | 无（必填）| 与脚本 1 相同的 Seurat rds |
| `-L, --input_loom` | 无（必填）| 脚本 2 产出的 `SCENIC.loom` |
| `-AC, --annotation_column` | 无（必填）| 细胞类型列名 |
| `-TF, --target_tfs` | 空 | 单独画哪些 regulon，逗号分隔（如 `TF1(+),TF2(+)`）；空则每类细胞自动取 RSS 前 `-NT` 个 |
| `-NT, --n_top` | 2 | 每类细胞取前几个 regulon 做点图与汇总热图 |
| `-OS, --one_set` | 空 | **只影响逐细胞类型排名图**（2×2 布局，逗号分隔）；空则取 max RSS 最高的 4 类。汇总热图始终覆盖全部细胞类型 |
| `-ZW, --z_threshold` | 3 | 传给 `plotRSS` 的 zThreshold |
| `-TH, --rss_threshold` | 0.1 | 传给 `plotRSS` 的 thr；低于此值的格子留白 |
| `-W, --width` | 10 | **仅** RSS 总图宽度（英寸）；多面板 regulon 图自动定尺寸 |
| `-H, --height` | 12 | **仅** RSS 总图高度（英寸）；多面板 regulon 图自动定尺寸 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# ① 准备 loom
./scRNA_pyscenic_prepare.r -I 05_annotation/05.seu_annotated.rds \
    -A RNA -L data -MG 50 -P 01. -O 12_pyscenic

# ② 跑三步核心算法（参考文件放在 ./ref_data）
./scRNA_pyscenic.py -I 12_pyscenic/01.expr_norm.loom -S mouse \
    -RD ref_data -T 12 -P 01. -O 12_pyscenic

# ②' 只重跑 ctx（比如 ctx 中途失败）
./scRNA_pyscenic.py -I 12_pyscenic/01.expr_norm.loom -S mouse \
    -RD ref_data -ST ctx -T 12 -P 01. -O 12_pyscenic

# ③ 算 RSS 并出图
./scRNA_pyscenic_RSS.r -I 05_annotation/05.seu_annotated.rds \
    -L 12_pyscenic/01.SCENIC.loom -AC celltype \
    -TF 'TF1(+),TF2(+)' -NT 2 -P 02. -O 12_pyscenic/RSS
```

不加 `-TF` 时，脚本会自动挑「每类细胞 RSS 前 `-NT` 个」的 regulon 来画小提琴图/Feature/Dot。

### 输出

**`scRNA_pyscenic_prepare.r`**

| 文件 | 说明 |
|---|---|
| `<前缀>expr_norm.loom` | pySCENIC 输入（`Gene` = 行，`CellID` = 列）|

**`scRNA_pyscenic.py`**

| 文件 | 说明 |
|---|---|
| `<前缀>adj.csv` | 共表达邻接表（`TF, target, importance`）|
| `<前缀>regulons.csv` | regulon 表（每个 TF × motif × 靶基因列表）|
| `<前缀>SCENIC.loom` | AUCell 结果（含 `Regulons` / `RegulonsAUC` / `RegulonsAUCThreshold`）|

**`scRNA_pyscenic_RSS.r`**

| 文件 | 说明 |
|---|---|
| `<前缀>rss.csv` | RSS 矩阵（regulon × 细胞类型）|
| `<前缀>RSS_plot.png/.pdf` | RSS 总图：行 = regulon、列 = 细胞类型，点大小 = RSS、颜色 = Z |
| `<前缀>RSS_oneSet.png/.pdf` | 逐细胞类型的 RSS 排名（2×2 布局）|
| `<前缀>regulon_violin.png/.pdf` | 目标 regulon 的活性小提琴图（按细胞类型分组）；面板自动排成接近正方形的网格 |
| `<前缀>regulon_feature.png/.pdf` | 目标 regulon 活性投到降维图上；同样自动排成接近正方形的网格 |
| `<前缀>regulon_dotplot.png/.pdf` | 目标 regulon 的活性点图：颜色 = 平均 AUC，点大小 = AUC 超过该 regulon 全局 75 分位的细胞比例；宽度与高度随 regulon / 细胞类型数自动调整 |
| `<前缀>seu_SCENIC.rds` | 已把「细胞 × regulon」活性矩阵加进 metadata 的 Seurat 子集 |

示例输出：`Images/pyscenic/`。

### 注意事项

- **三个脚本必须按顺序跑**，中间产物（loom）不能跨物种复用。
- **Seurat v5 的拆分层会自动合并**：若 assay 里是 `counts.1` / `counts.2` 这种拆分层（merge 后常见），脚本 1 会自动 `JoinLayers` 再取矩阵，不需要事先手动合并。
- **图件尺寸自适应**：`regulon_violin` / `regulon_feature` 按面板数自动排成**接近正方形**的网格（面板数 n 时取 `ncol = ceil(sqrt(n))`），不会因为 regulon 多而变成一条长条；`regulon_dotplot` 宽度随 regulon 数加宽。`-W` / `-H` **只作用于 RSS 总图**。
- **5 张结果图各有分工，配合着读**：找标志物看 `RSS_plot`（**又大又深绿**的点）、看排名看 `RSS_oneSet`、用 `regulon_dotplot`（颜色 = 平均 AUC、大小 = 活跃细胞比例）与 `regulon_violin` 交叉验证绝对活性、用 `regulon_feature` 看连续性。详见「结果图怎么读」一节。
- **三张 `regulon_*` 图画的是 AUC，不是靶基因表达量**：`regulon_violin` / `regulon_feature` / `regulon_dotplot` 用的都是写回 metadata 的 **AUCell 活性分**。AUC 由「靶基因在该细胞表达排序中的位置」聚合而来，与靶基因的绝对表达量不是一回事。
- **`regulon_dotplot` 没有沿用 Seurat `DotPlot` 的默认统计量**：默认的「AUC > 0 的细胞比例」几乎恒为 100%（实测中位数 0.990、一半 regulon 为 1.000），点大小会失去信息。改用「AUC 超过该 regulon 全局 75 分位的细胞比例」，该比例在全细胞范围恒为 **0.25**，在单个细胞类型内才拉开差异，**0.25 即中性线**。
- **`RSS_plot` 的大小与颜色编码不同的量**：大小 = RSS 绝对值，颜色 = Z（**按细胞类型标准化**，即在该类细胞内部比较 regulon）。**高 Z 但小点 = 矮子里拔将军，不能当标志物**。
- **环境由 shebang 决定**：`scRNA_pyscenic.py` 必须用 pySCENIC 环境的 python 执行。若用别的 python 跑，脚本会在找 `pyscenic` CLI 时直接报错（CLI 是按 `sys.executable` 同目录定位的，不在 PATH 里）。
- **`pyscenic` CLI 不支持 `python -m pyscenic`**（该包没有 `__main__`），必须调 bin 目录下的可执行文件。
- **ctx 是最慢的一步**，可能几小时。输出已存在时默认跳过；要重跑加 `-F`，或只跑 `-ST ctx`。
- **`SCENIC.loom` 不能用 Python 的 `loompy` 读**——不符合 loom v3 规范（多列属性 dtype 不合法）。只能在 R 里用 `SCopeLoomR::open_loom` 读。这是流程「必须回 R」的根本原因。
- **`-NC` 默认不抽样**（保留全部细胞），避免静默丢数据。它的主要用途是**加速**：细胞数过万时 pySCENIC 会跑很久，可设 `-NC 2000` 等值抽样——但抽样会漏掉稀有细胞类型（某类只有几十个细胞时可能只剩个位数，RSS 不可靠），这类数据不要抽样。
- **regulon 数量取决于数据**：质控差、细胞数少时 `ctx` 可能只输出极少甚至零个 regulon；零个时 `aucell` 无内容可打分，脚本会给出提示。
- **`(+)` / `(-)` 是调控方向**：实际几乎全是 `(+)`（激活型）。抑制型需要额外的负相关判据，默认参数下很少出现。
- **RSS 数值天然偏小**（JSD 性质），要横向比较而非看绝对值；`plotRSS` 的 `-TH` 阈值以下的格子会被留白，看起来「空」是正常的。
- **RSS ≠ 显著性**：它是分布相似度，没有 p 值。判断某个 TF 是不是某类细胞的标志，看它在那一列的 RSS 是否明显高于其它列。
- **`calcRSS` 内部用 zThreshold 做可视化分组**，与统计检验无关。
- **AUCell 的 AUC 与 ROC 的 AUC 无关**，别混用。
- **输出体积**：`adj.csv` 可达上百 MB（实测 300 细胞时 47 MB，2000 细胞时约 68 MB），`seu_SCENIC.rds` 数十 MB，注意磁盘空间。
- **pySCENIC 版本差异**：本流程按 pySCENIC 0.12.x（Python 3.7）验证；新版本 CLI 参数可能变化，尤其是 `ctx` 的 `--mask_dropouts` 默认值在 0.9.16 之后改过（旧版默认开、新版默认关），所以脚本显式传参而不依赖默认值。


---

## 29 scRNA_hdWGCNA.r（独立脚本）—— 高维共表达网络分析（hdWGCNA）

### 功能

**在单个细胞类型内部找出「哪些基因一起变化」，并把它们归成模块，再从每个模块里挑出最核心的 hub 基因。**

回答三个问题：

1. **这个细胞类型里，基因可以分成几组共表达模块？**（模块 = 一群表达高度同步的基因，通常对应同一条通路或同一个调控程序）
2. **每个模块里谁是「hub」？**（模块内连接度最高的基因，往往是这个程序的关键调控者）
3. **这些模块在不同样本/分组之间活性有没有差别？**（用模块特征值 hMEs 在细胞上打分）

**与普通 WGCNA 的区别**：单细胞数据 dropout 严重，直接做 WGCNA 会得到充满噪声的模块。hdWGCNA 先把细胞按「细胞类型 × 样本」聚成 **metacell**（每 25 个邻近细胞合成一个伪细胞），在 metacell 上建网络，从而大幅降噪。所以它**必须先指定一个目标细胞类型**，一次只分析一类。

**与其它脚本的分工**：

| 脚本 | 回答什么 |
|---|---|
| `scRNA_clustering_step4.r` | 细胞之间有哪几群（cluster）|
| `scRNA_gene_DE.r` / `scRNA_pseudobulk_DE.r` | 两组之间哪些基因变了 |
| `scRNA_GSVA_enrichment.r` / `scRNA_genesetScore.r` | 已知通路的活性高低 |
| **`scRNA_hdWGCNA.r`** | **在已知通路之外，数据自己能把基因分成哪几组、谁在组里最核心** |

### 核心概念

整条流程可以串成一句话：

```
共表达（哪些基因同涨同落）
  → 取幂（把弱相关压掉、强相关保留）
  → TOM（两个基因共享多少邻居）
  → 切模块（一组一起变化的基因）
  → kME（模块里谁最核心 = hub）
  → ME（把整个模块压成一个数，给每个细胞打分）
```

#### 术语速查

| 术语 | 一句话含义 |
|---|---|
| **metacell** | 把 k 个邻近细胞平均成一个「伪细胞」，用来抹平 dropout |
| **共表达** | 两个基因在细胞之间同涨同落 |
| **soft power（软阈值）** | 把相关系数取幂 `abs(r)^power`，压低弱相关、放大强相关 |
| **TOM** | 基因之间的拓扑重叠度：两个基因共享的邻居越多，值越高 |
| **模块（module）** | 一群表达高度相关的基因，用颜色命名；`grey` = 未归入任何模块 |
| **kME** | 某个基因与它所在模块整体活性的相关度；越高说明这个基因越能代表该模块 |
| **ME** | module eigengene，模块特征值：把整个模块压成**一个数**，代表该模块在某细胞里的活跃程度 |
| **hME** | 校正后的 ME：去掉与其它模块的重叠，可跨模块比较 |
| **hub gene** | 模块里 kME 最高的基因 |

#### 1. 共表达是什么

两个基因在细胞之间「同涨同落」就叫共表达：

```
细胞:      c1   c2   c3   c4   c5
基因A:     5.2  0.1  4.8  0.3  5.5
基因B:     4.9  0.2  5.1  0.4  5.2    ← 和 A 几乎同步
基因C:     0.1  4.7  0.2  5.0  0.3    ← 和 A 无关
```

A 和 B 高度共表达。**共表达的基因通常因为被同一个转录因子、或同一条通路共同调控**，所以共表达模式能揭示「程序」。

**为什么模块有意义**：一组成员基因表达同步 → 找到模块 = 找到这个细胞类型里**在运作的生物学程序**；找到 hub = 找到这个程序的**关键节点**。

#### 2. 为什么必须先做 metacell

单细胞测序有大量基因测不到（表达量低就记为 0），这叫 dropout。后果是**真实的共表达被这些假的 0 稀释**，相关系数被低估，直接做 WGCNA 会得到充满噪声的模块。

解法：把 k 个「表达谱相近」的细胞**平均成一个伪细胞**（`k = 25`，即 25 个细胞合成 1 个 metacell）。

好处有三条：
- dropout 被平均掉（一个细胞漏了，其它 24 个会补上）
- 细胞数大幅减少，计算量从「几万细胞」降到「几百 metacell」，WGCNA 才跑得动
- 因为是在**同一细胞类型、同一样本内部**聚类，细胞类型之间的差异没有被抹掉

`max_shared = 10` 允许不同 metacell 共用最多 10 个细胞；`min_cells = 100` 表示某样本里该细胞类型不足 100 个就跳过（100 个细胞只能合成 4 个 metacell，太少）。

#### 3. 软阈值：取幂到底改变了什么

普通相关系数 r 在 −1 到 1 之间。WGCNA 不直接用，而是取幂：

```
邻接强度 = abs(r) ^ β        β = soft power
```

效果（以 β=6 为例）：

| 真实相关 r | r^6 | 结果 |
|---|---|---|
| 0.30（弱，多半是噪声）| 0.0007 | **几乎归零** |
| 0.60（中等）| 0.047 | 很低 |
| 0.90（强）| 0.53 | 保留 |

**取幂的作用就是压低弱相关、放大强相关**，把网络变稀疏。

β 怎么选：希望网络成为「无标度网络」（少数基因连接极多、多数基因连接很少，这是真实生物网络的特征）。看 `softPowers.png` 左上角，**第一个让 R² 越过 0.8 的 β 就是合适的值**。默认值 `3` 只是占位，**每份数据都必须重新看图确定**。

- β 太小 → 弱相关都保留 → 网络过密、模块糊在一起
- β 太大 → 只剩极强相关 → 模块碎成小块、大量基因掉进 grey

#### 4. TOM：为什么不直接用相关系数

```
TOM_ij = (两个基因共同邻居的权重之和 + 直接权重) / 归一化分母
```

含义：**两个基因如果共享很多「邻居基因」，它们的连接就更可信。**

为什么需要：两个基因可能碰巧直接相关（比如都受细胞大小影响），但它们如果在网络里没有共同邻居，说明不在同一个程序里。TOM 把「是否共享邻居」这层结构考虑进去，比单独看相关系数稳健得多。脚本把 TOM 存到 `TOM/` 目录，文件大小随基因数平方增长，可能很大。

#### 5. 模块是怎么切出来的

有了 TOM，把基因按「距离 = 1 − TOM」做**层次聚类**（就是 `Dendrogram.png` 那棵树），再用动态剪切树把树枝剪成几段——**每一段就是一个模块**。

- **模块数不是固定的**，由算法根据数据结构决定，同一份数据换 β 就会变
- 模块用颜色命名（turquoise / blue / brown / …）
- **`grey` 是「没归入任何模块」的基因**，通常是噪声或孤立基因，所有分析都排除它
- `gene_modules.csv` 里只列出**已剔除 grey** 的基因，这些才是下游分析的基数

模块大小差异通常很大（大的几百个基因、小的只有几十个），这不是异常，反映的是各程序涉及的基因数量不同。

#### 6. kME 与 hub 基因

```
kME = 某个基因的表达 与 它所在模块的整体活性 的相关系数
```

- kME 接近 1 → 这个基因几乎能代表整个模块
- kME 接近 0 → 这个基因虽在模块内，但只是边缘成员

**hub 基因 = 模块里 kME 最高的那些。** 想拿 hub 清单，从 `gene_modules.csv` 里按 `kME_<模块名>` 列降序取前 10 个（`plotKMEs.png` 每个面板方框里标的就是这 10 个）。

#### 7. ME：怎么把几百个基因变成一个数

**先说清楚：被压缩的不是基因，而是「一个细胞在这个模块上的表达数据」。**

一个模块有 494 个基因 → **一个细胞**在这 494 个基因上就有 **494 个表达值** → ME 把它们换成 **1 个数值**。基因本身没有消失。

```
细胞 c1 在这个模块上：
  GeneA=5.2  GeneB=4.9  GeneC=3.1  ...  第494个基因=2.8
  └───────────── 494 个数字 ─────────────┘
                     ↓  压缩
                  ME = 6.63        ← 1 个数字
```

**严格说，ME 是模块所有基因的第一主成分「得分」**——每个细胞一个投影值，而不是第一主成分本身（那是每个基因一个权重）。eigengene 直译「特征基因」，相当于虚构一个基因，它在每个细胞里的表达值就等于整个模块的综合水平。

**类比**：细胞 = 学生，模块里的基因 = 这个学生的各门课成绩，ME = 「综合能力分」。494 门课没法直接比谁强，算个综合分就比得了。

**用 4 个基因的小例子实算一遍**（真实计算，非示意）：

原始表达矩阵（行=基因，列=细胞），4 个基因同涨同落：

```
          c1    c2    c3    c4    c5
GeneA    5.2   0.1   4.8   0.3   5.5
GeneB    4.9   0.2   5.1   0.4   5.2
GeneC    3.1   0.1   2.9   0.2   3.3
GeneD    2.2   0.1   2.1   0.1   2.4
```

第一步，给每个基因一个「权重」（由主成分分析从数据里学出来，全为正 = 朝同一方向变化）：

```
GeneA  +0.639    GeneB  +0.611    GeneC  +0.379    GeneD  +0.274
```

第二步，每个细胞算出 1 个数：

```
c1   +3.025
c2   -4.816
c3   +2.789
c4   -4.528
c5   +3.531
```

**这 5 个数就是 5 个细胞的 ME。** c2 最负、c1 最正，正对应「c2 的 4 个基因都很低、c1 都很高」。

第三步，验证它确实代表了原来那 4 个基因：

```
ME 与 GeneA 的相关性  r = +0.999
ME 与 GeneB 的相关性  r = +0.999
ME 与 GeneC 的相关性  r = +1.000
ME 与 GeneD 的相关性  r = +0.999
第一主成分解释的方差 = 99.8%
```

**一个数值，和原来 4 个基因的相关系数全在 0.999 以上**——信息基本没丢。「压缩」不是丢信息，而是抓住它们共同的那一个趋势。

**为什么用主成分而不是直接取平均**：模块里若混进一个噪声基因（表达忽高忽低），简单平均会被它拖偏。主成分的权重**从数据里学**——与整体趋势一致的基因权重高，不稳定的基因权重被压到很低甚至为 0，所以更稳健。

**导出后的样子**：`MEs.csv` 每行一个细胞、每列一个模块（含 grey），**每个格子一个 ME**——每一列都是由该模块的全部基因压缩来的一个数。

#### 8. hME：为什么还要再校正一次

模块之间会互相「串味」。比如 `NEW_1` 和 `NEW_2` 的基因有重叠，那 `NEW_1` 的 ME 里就混进了 `NEW_2` 的信号——你看到 `NEW_1` 高，可能其实是被 `NEW_2` 带上去的。

**hME 的解法**：把其它模块的影响回归掉，只留这个模块独有的部分。所以 **hME 才能跨模块直接比较**（「`NEW_1` 和 `NEW_3` 哪个更活跃」这种问题必须用 hME）。

`MEs.csv` 与 `hMEs.csv` 行列结构完全相同（同样的细胞、同样的模块），只是数值不同。

拿同一个细胞对照两张表会发现：**方向一致、hME 的绝对值略小**。这是因为 hME 去掉了从其它模块「借来」的那部分信号，所以更保守。做跨模块比较、或想排除模块间串扰时，用 hME。

#### 9. 三个必须记住的边界

1. **模块 ≠ 通路**：模块是数据驱动的共表达群，可能对应一条通路，也可能只是「一组高表达基因」或技术因素造成的。必须拿基因列表去做富集才知道它是什么。
2. **hub ≠ 功能最重要**：hub 只是网络里连接最多，不代表功能上最关键。要结合 `gene_DE` 的差异倍数和富集结果判断。
3. **ME / hME 只能比相对高低**：绝对值没有意义（不同模块尺度不同、正负号是算法任意定的）。可以说「这个细胞的 hME 比那个细胞高」，不能说「hME = 6.63 代表活性 6.63」。

### 输入

一个 **Seurat rds 对象**，且必须满足：

| 要求 | 原因 |
|---|---|
| 已注释（有细胞类型列，由 `-AC` 指定）| 脚本按细胞类型筛选目标群体 |
| 有样本列（由 `-SC` 指定）| metacell 按「细胞类型 × 样本」分组构建 |
| **有 `harmony` 降维**（`seu@reductions$harmony`）| metacell 构建时 `reduction = "harmony"` 硬编码，缺了直接失败 |
| 已归一化 | 脚本会再跑一次 `NormalizeData(LogNormalize)` |

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| Seurat rds 路径 |
| `-AC, --annotation_column` | 无（必填）| 细胞类型所在的 `meta.data` 列名 |
| `-SC, --sample_column` | 无（必填）| 样本所在的 `meta.data` 列名 |
| `-TC, --target_cell` | 无（必填）| **目标细胞类型名**（`-AC` 列里的一个取值）|
| `-SP, --soft_power` | **3** | WGCNA 软阈值。**先看 `softPowers.png` 再定**，见注意事项 |
| `-A, --assay` | RNA | 用哪个 assay 的表达矩阵 |
| `-T, --threads` | 8 | WGCNA 线程数 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录（不存在会自动创建）|

### 用法示例

```bash
# 第一步：先用默认 soft power 跑一遍，重点看 softPowers.png
./scRNA_hdWGCNA.r -I 05.seu_annotated.rds -AC celltype -SC orig.ident \
    -TC 'CellTypeA' -A RNA -SP 3 -T 8 -P 01. -O 13_hdWGCNA

# 第二步：根据 softPowers.png 重新指定 power（图会覆盖，但 TOM 会重算）
./scRNA_hdWGCNA.r -I 05.seu_annotated.rds -AC celltype -SC orig.ident \
    -TC 'CellTypeA' -A RNA -SP 6 -T 8 -P 01. -O 13_hdWGCNA
```

**如果要分析多个细胞类型，就换 `-TC` 跑多次**，每次给不同的 `-P` 前缀，避免产物互相覆盖。

### 输出

**图（每张都出 PNG + PDF）**

| 文件 | 内容 |
|---|---|
| `<前缀>metacell_<细胞类型列名>umap.png` | metacell 按细胞类型着色的 UMAP（检查 metacell 是否保留了细胞类型的结构）|
| `<前缀>metacell_<样本列名>umap.png` | metacell 按样本着色的 UMAP（检查有没有残留样本批次效应）|
| `<前缀>softPowers.png` | **软阈值诊断图**：4 个面板分别是 scale-free topology fit、mean / median / max connectivity 对 power 的曲线，黑圈标出脚本选用的 power |
| `<前缀>Dendrogram.png` | **基因聚类树 + 模块颜色条**：树的分支是基因间的距离，底部色条是模块归属（灰色 = 未分入模块）|
| `<前缀>plotKMEs.png` | 各模块内基因的 kME 排序（找 hub 用）|
| `<前缀>cell_module.png` | 每个模块的 hME 打分投到 UMAP 上（看模块在哪些细胞里活跃）|
| `<前缀>cell_hubGenes.png` | 每个模块的 **hub 基因打分**（UCell）投到 UMAP 上 |
| `<前缀><目标细胞>_radar.png` | 雷达图：各模块 hME 按样本比较，一眼看出哪个样本哪个模块高 |
| `<前缀>module_correlation.png` | 模块之间的相关性（模块特征值两两相关）|
| `<前缀>module_celltype_dotplot.png` | 细胞类型 × 模块的 dotplot（hME 均值），看模块的细胞类型特异性 |

**表**

| 文件 | 内容 |
|---|---|
| `<前缀>gene_modules.csv` | **基因 → 模块**归属表（已剔除 grey 模块），含 kME、module、color 等列 |
| `<前缀>top10_hubgenes.csv` | 每个模块 kME 前 10 的 hub 基因 |
| `<前缀>MEs.csv` | 细胞 × 模块的**原始**特征值 |
| `<前缀>hMEs.csv` | 细胞 × 模块的**校正后**特征值（跨模块可比）|

**对象与中间文件**

| 文件 | 内容 |
|---|---|
| `<前缀>seu_hdWGCNA.rds` | **最终 Seurat 对象**，所有 hdWGCNA 结果都在里面（模块、kME、hMEs、metacell 等）|
| `<outpath>/TOM/<目标细胞>_TOM/` | TOM 中间文件（`*_TOM.rda`），**体积大**；后续重跑同参数会复用，`overwrite_tom=TRUE` 时重算 |

示例输出：`Images/hdWGCNA/`。

### 结果图怎么读

十张图分两类：**前 4 张是「质检图」，用来判断这次运行的结果可不可信；后 6 张是「结果图」，用来回答生物学问题。**
质检图不对，后面的结果就不用看了——所以必须按顺序读。

| # | 图 | 类别 | 作用 | 怎么读 |
|---|---|---|---|---|
| 1 | `metacell_<AC值>umap.png` | 质检 | metacell 是否保留了细胞类型结构 | 各细胞类型仍能分开 = 正常 |
| 2 | `metacell_<SC值>umap.png` | 质检 | metacell 有无残留样本批次效应 | 各样本充分混合 = 正常 |
| 3 | `softPowers.png` | 质检 | 软阈值 β 选得对不对 | 看 R² 第一个越过 0.8 的点 |
| 4 | `Dendrogram.png` | 质检 | 模块切得合不合理 | 分支清晰、grey 占比不过大 |
| 5 | `plotKMEs.png` | 结果 | 每个模块里 hub 的位置 | 顶部有明显「断崖」= 模块内聚性好 |
| 6 | `module_celltype_dotplot.png` | 结果 | 模块的细胞类型特异性 | 目标细胞类型的那一行整体偏红 = 模块属于它 |
| 7 | `cell_module.png` | 结果 | 模块活性在细胞间的分布 | 连续渐变 = 状态过渡；突然点亮 = 离散亚群 |
| 8 | `cell_hubGenes.png` | 结果 | 只用 hub 基因算的活性 | 与第 7 张形状一致 = 模块信号确实来自 hub |
| 9 | `module_correlation.png` | 结果 | 模块之间的协同/互斥 | 强正相关 = 同一程序的两个片段 |
| 10 | `<目标细胞>_radar.png` | 结果 | 模块活性在样本间的差异 | 某个样本某个模块明显凸出 = 候选生物学差异 |

#### 质检图（先看这 4 张，不对就别往下看）

**① `metacell_<AC值>umap.png` —— metacell 按细胞类型着色**

横竖轴是 metacell 的 UMAP 坐标，每个点是一个 metacell，颜色是细胞类型。

- **正常**：不同细胞类型仍各自成团、彼此分开
- **异常**：所有类型混成一坨 → metacell 把类型差异也平均掉了，后面找出来的模块没有意义（多半是 `k` 太大或 `-TC` 选错了）

**② `metacell_<SC值>umap.png` —— metacell 按样本着色**

同一批 metacell，改按样本上色。

- **正常**：各样本颜色均匀混杂
- **异常**：样本各自成团 → 残留批次效应。这时找出的模块可能反映的是「样本差异」而不是生物学程序（回去检查 `scRNA_Nor_harmony_step3.r` 的整合结果）

**这两张必须一起看**：类型分开 + 样本混匀 = metacell 建对了。这是整个流程的地基。

**③ `softPowers.png` —— 软阈值诊断（4 个面板）**

| 面板 | 内容 |
|---|---|
| 左上 **scale free topology fit** | 横轴 β、纵轴 R²。**这是唯一需要看的** |
| 另三张 | mean / median / max connectivity 随 β 的变化（辅助参考）|

- 图中有一条 **0.8 的灰色虚线**
- 黑圈标出脚本当前用的 `-SP`
- **正确做法**：找**第一个让 R² 越过 0.8** 的 β，用那个值

| 现象 | 含义 | 后果 |
|---|---|---|
| β 太小 | 弱相关全被保留 | 网络过密、模块糊在一起、模块数偏少但每个都很大 |
| β 太大 | 只剩极强相关 | 模块碎成很多小块、大量基因掉进 grey |
| R² 一直上不去 0.8 | 数据本身不满足无标度假设 | 常见于细胞数太少或异质性太强，结论要谨慎 |

**④ `Dendrogram.png` —— 基因聚类树 + 模块颜色条**

上半部分是基因的层次聚类树（纵轴是基因间距离 = 1 − TOM），下半部分是色条，每段颜色对应一个模块。

- **正常**：树枝在高层分成几个清晰的簇，色条与树枝的簇大体对齐
- **grey 占比过大**（灰色段特别长）→ 说明大量基因没被分进任何模块，可能是 β 太大或数据噪声太高
- grey 段往往占色条的很大一部分，这是常见现象；注意**只有非 grey 基因才进入下游分析**

#### 结果图（按 5 → 10 的顺序读）

**⑤ `plotKMEs.png` —— 每个模块内基因按 kME 排序**

每个面板是一个模块，横轴是模块内基因的排名，纵轴是 kME。

- **看曲线顶部的形状**：最高处与后面明显拉开（一个「断崖」），说明这个模块有明确的 hub、内聚性好
- **曲线整体平缓**：说明模块内没有特别突出的基因，模块可能只是「一群差不多相关的基因」拼起来的
- 一个模块**同时有好几个并列最高**（不是单一峰）→ 可能实际上是两个程序被合并了

**⑥ `module_celltype_dotplot.png` —— 细胞类型 × 模块**

行 = 细胞类型，列 = 模块，颜色 = hME 均值（红高蓝低），点大小 = 表达比例。

- **这是判断「模块归谁」的主要依据**
- **正常**：每个模块在某一两类细胞里明显偏红（模块有细胞类型特异性）
- **异常**：所有模块在所有细胞类型上都是浅色 → 模块没有特异性，可能只是反映了「测序深度」或「细胞大小」这类技术因素
- **重点看 `-TC` 那一行**：这一行里哪些列偏红，就是属于这个细胞类型的模块，优先拿去做富集

**⑦ `cell_module.png` —— 每个模块的 hME 投到细胞 UMAP 上**

每个面板一个模块，颜色 = hME。**注意：这里的点是原来的细胞，不是 metacell。**

- **连续渐变**（从一块区域平滑过渡到另一块）→ 这个模块活性随某种状态连续变化，可能对应分化或代谢梯度
- **突然点亮**（某群细胞整块亮、周围暗）→ 这个模块是某个亚群特有的，跟一个离散的细胞状态对应
- **到处都亮** → 模块不特异，回去看第 6 张图确认

**⑧ `cell_hubGenes.png` —— 只用 hub 基因算的活性**

和第 7 张形状一样，但数值来自**每个模块 kME 前 25 个基因**（UCell 打分），而不是整个模块。

- **用途是验证**：如果第 7 张和第 8 张形状基本一致，说明**模块的信号确实集中在 hub 基因上**，hub 的选择可信
- **两者形状差别很大** → 模块的活性不是由少数 hub 驱动的，而是由一大批边缘基因撑着。这种情况下 hub 基因的生物学意义要打折扣

**⑨ `module_correlation.png` —— 模块之间的相关性**

模块特征值两两相关的热图。

- **强正相关**（深红）→ 两个模块很可能属于**同一个生物学程序的不同片段**，被算法切成了两块。可以合并起来一起做富集
- **强负相关**（深蓝）→ 两个程序此消彼长，常见于「分化前状态 vs 分化后状态」
- **接近 0**（浅色）→ 两个模块独立，是两条不相关的过程

**⑩ `<目标细胞>_radar.png` —— 模块活性的样本比较**

雷达图，每个轴是一个模块，每条线是一个样本，半径是该样本里目标细胞的 hME 均值。

- **用途**：看模块活性在样本（或分组）之间有没有差别
- **某个样本在某个轴上明显凸出** → 那个样本里这个模块特别活跃，是候选的生物学差异
- ⚠️ **这只是可视化，没有统计检验**。要下结论必须用 `hMEs.csv` 自己做检验（而且要以样本为单位，不能以细胞为单位，否则是伪重复）

#### 一句话读图顺序

```
【先质检，任何一张不对就停下】
1. metacell_<AC>umap.png       → 细胞类型还分得开吗？
2. metacell_<SC>umap.png       → 样本混匀了吗？
3. softPowers.png              → β 是否是第一个 R²>0.8 的点？
4. Dendrogram.png              → 模块切得干净吗？grey 是否过多？

【再看结果】
5. plotKMEs.png                → 每个模块有明确 hub 吗？
6. module_celltype_dotplot.png → 模块归哪个细胞类型（重点看 -TC 那一行）
7. cell_module.png             → 模块活性是连续渐变还是突然点亮
8. cell_hubGenes.png           → 与第 7 张对照，验证信号来自 hub
9. module_correlation.png      → 哪些模块是同一程序的两半
10. <目标细胞>_radar.png       → 样本间有无差异（仅可视化，需另外检验）
```

#### 常见异常信号对照

| 看到什么 | 大概率原因 | 怎么处理 |
|---|---|---|
| 第 1 张图细胞类型混成一坨 | `k` 太大，或 `-TC` 选错 | 减小 `k`（需改代码，硬编码 25）|
| 第 2 张图样本各自成团 | 上游去批次不充分 | 回 `scRNA_Nor_harmony_step3.r` 检查整合 |
| 第 3 张图 R² 始终到不了 0.8 | 数据细胞数太少 / 异质性过高 | 结论降级为探索性 |
| 第 4 张图 grey 占一大半 | β 太大，或基因过滤太宽松 | 调小 `-SP`；或提高 `fraction`（硬编码 0.05）|
| 模块数只有 2–3 个且都很大 | β 太小，网络过密 | 调大 `-SP` |
| 第 6 张图全是浅色 | 模块反映技术因素 | 检查是否该做 `regress_out`；结论谨慎 |
| 第 7、8 张形状不一致 | hub 不具代表性 | 不要单凭 hub 基因下结论 |


### 注意事项

- **输入对象必须有 `harmony` 降维**：第 156 行 `reduction = "harmony"` 是硬编码的，没有这个 reduction 会在构建 metacell 时直接失败。本库主流程 `scRNA_Nor_harmony_step3.r` 的产物就带 `harmony`，可以直接用。
- **`-TC` 只接受一个细胞类型**，且必须是 `-AC` 列里真实存在的取值。脚本不会自动校验，写错会在 `SetDatExpr` 处报「群组不存在」。
- **soft power 必须先看图**：默认 `3` 只是这份测试数据的值。正确做法是跑到 `softPowers.png` 后，找**scale-free topology fit 第一个越过 0.8（图中虚线）**的 power，再回填 `-SP`。power 太小 → 网络过于稠密、模块糊成一团；太大 → 模块全碎、只剩 grey。
- **目标细胞类型在每个样本里至少要 100 个细胞**（`min_cells = 100` 硬编码）。细胞数不足的样本会被跳过；如果所有样本都不足，metacell 建不起来会报错。稀有细胞类型不适合做这个分析。
- **metacell 的两张 UMAP 分别按 `-AC` 和 `-SC` 着色**，文件名也用对应的列名：`metacell_<AC值>umap.png` 按细胞类型着色（看 metacell 是否保留类型结构），`metacell_<SC值>umap.png` 按样本着色（看有没有残留批次效应）。两张图对照着看：类型分得开、样本混得匀，说明 metacell 建得合适。
- **雷达图只画 `-TC` 那一个细胞类型**：脚本按 `-AC` 列筛出该类型的细胞再画，所以 `-AC` 与 `-TC` 必须指向同一套注释。
- **几个关键数值写死在代码里，没有对应参数**：`fraction = 0.05`（基因需在 ≥5% 的 metacell 中表达）、`k = 25`（每个 metacell 由 25 个细胞合成）、`max_shared = 10`、`min_cells = 100`、`n_genes = 25`（给细胞打模块分时，每个模块取 kME 前 25 个基因）。想改只能改代码。
- **⚠️ 脚本里有两个「25」，含义完全不同，别混**：
  - 第 157 行 `k = 25` —— **metacell 粒度**：每 25 个邻近细胞平均成 1 个伪细胞（影响降噪程度与网络构建）
  - 第 259 行 `n_genes = 25` —— **hub 打分粒度**：给每个细胞算模块分时，只用该模块 kME 最高的 25 个基因（决定 `cell_hubGenes.png`）
  
  两者改动的后果互不相干：改 `k` 会改变 metacell、进而改变整个网络与模块划分；改 `n_genes` **只影响 hub 打分这一路结果**（`cell_hubGenes.png` 和对象里存的 score），不改变模块划分，也不改变 `MEs` / `hMEs`（后者由 `ModuleEigengenes` 用整个模块算出）。
- **`plotKMEs.png` 方框里标的是 kME 前 10 个基因**（不是 25），这是绘图函数的内部设定，与上面的 `n_genes` 无关。
- **会重新做一次归一化**：第 142 行 `NormalizeData(...)`，所以输入对象即使是归一化过的也会再归一化一遍（重复 log 化会改变尺度）。若输入已是 log 归一化数据，建议先确认这一步是否符合预期。
- **`seu_hdWGCNA.rds` 会很大**：metacell 对象与全部模块信息都塞在里面，体积随细胞数增长，可能达到 GB 级；`TOM/` 目录同样很大。磁盘要留够。
- **`overwrite_tom = TRUE`**：每次重跑都会重算 TOM，换 `-SP` 时这是必要的；但如果只是想重画图，会白等一次 TOM 计算。
- **第 294 行的日志信息写错了**：实际保存的是 `cell_hubGenes.png`，但打印的提示是 `cell_module.png`。不影响产物，只是日志会误导。
- **radar 图依赖 `ggradar`，hub 基因打分依赖 `UCell`**：两个包缺失时脚本会尝试 `install.packages`，离线环境会失败。
- **模块颜色会重命名**：脚本用 `ResetModuleNames(new_name = "<目标细胞>_NEW_")` 把模块改成 `<细胞类型>_NEW_<颜色>`，所以 `gene_modules.csv` 里的模块名和 `Dendrogram.png` 里的颜色不是同一套写法，对照时注意。
- **hub 基因只是「模块内连接度最高」**，不等于功能最重要；要结合 `gene_DE` / 富集结果看。


---

## 30 scRNA_scTenifoldKnk.r（独立脚本）—— 单细胞虚拟基因敲除（scTenifoldKnk）

### 功能

**在计算机里把某个基因「敲掉」，预测哪些基因会跟着受影响——不需要做真实的敲除实验。**

回答三个问题：

1. **敲掉目标基因，哪些基因的调控关系会明显改变？**（`diffRegulation` 表）
2. **这些基因富集在哪些通路上？**（`KO_effect_plot` 的网络节点配色与图例）
3. **受影响最强的基因是哪些？**（FC 排序 top20）

**与其它脚本的分工**：

| 脚本 | 回答什么 | 需要什么 |
|---|---|---|
| `scRNA_pseudobulk_DE.r` | 两组真实样本之间哪些基因表达变了 | 有分组的多组样本 |
| `scRNA_hdWGCNA.r` | 某细胞类型里基因分成哪些共表达模块 | 单一细胞类型、足够细胞数 |
| **`scRNA_scTenifoldKnk.r`** | **如果拿掉这个基因，网络会怎么变** | **一份表达矩阵 + 一个目标基因** |

**最大的特点：不需要敲除组数据。** 只有野生型（WT）单细胞数据也能做，靠的是从数据里重建调控网络，再在网络上做「虚拟敲除」，比较敲除前后的网络差异。适合**没有条件做基因敲除、但想知道某个基因影响力**的场景。

### 核心概念

#### 虚拟敲除（in silico knockout）是什么

不是改数据，而是**改网络**：

```
① 用 WT 单细胞数据重建基因调控网络      → WT 网络
② 把目标基因在网络里的所有连接置零        → KO 网络（这就是「敲除」）
③ 比较两个网络的差异                      → 哪些基因受影响
```

所以它测的是「**这个基因在网络拓扑上有多重要**」，不是「表达量会变多少」。

#### 为什么要建 10 个子网络（`nc_nNet`）

单细胞网络噪声大，一次建网不可靠。做法是**随机抽若干细胞建多个子网络**，再做张量分解把共同结构提出来、压掉噪声。

- `nc_nNet = 10`：建 10 个子网络
- `nc_nCells = 500`：每个子网络随机抽 500 个细胞

**这 10 个子网络还有第二个用途**：既当噪声平滑，又当**重复**——用它们之间的离散程度估计 Z 值和 p 值（相当于做了 10 次重复的置换检验）。所以子网络数不能太少。

#### 输出表里四列的含义

| 列 | 含义 |
|---|---|
| `gene` | 被检验的基因 |
| `distance` | WT 网络与 KO 网络之间这个基因的**流形距离**——距离越大，说明敲除后这个基因在网络里的位置变化越大 |
| `Z` | 该距离相对 10 个子网络间波动的**标准化偏离程度**。Z 越大越显著 |
| `FC` | 敲除前后距离的**变化倍数**，脚本按它排序取 top20 |
| `p.value` / `p.adj` | 由子网络间的分布算出的显著性，`p.adj` 是校正后（脚本按 `p.adj` 判显著）|

**注意这张表的逻辑和常规 DE 表完全不同**：DE 表比的是「两组表达量」，这里比的是「**网络拓扑**」。所以**一个基因完全可以表达量不变、但被判定为显著受影响**。

### 输入

一个 **Seurat rds 对象**，脚本会自己取两样东西：

| 取什么 | 从哪取 | 用途 |
|---|---|---|
| **counts 矩阵** | `GetAssayData(seu, layer = "counts")` | **必须是原始 counts**，不能是归一化数据 |
| **HVG 列表** | `VariableFeatures(seu)` | 只取高变基因 + 目标基因参与建网 |

**前提**：对象的 `VariableFeatures` 必须已经算过（跑过 `FindVariableFeatures`，即主流程 step4 之后的产物）。

脚本会自动把目标基因加入 HVG 列表（`unique(c(target_gene, HVGs))`），并**先校验目标基因是否在 counts 矩阵里**，不在则直接报错退出。

### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| Seurat rds 路径 |
| `-TG, --target_gene` | 无（必填）| **要敲除的基因名**，必须在 counts 矩阵里 |
| `-T, --threads` | 8 | 并行线程数 |
| `-p_val` | 0.05 | **`p.adj` 阈值**，判定「显著受影响基因」 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |

### 用法示例

```bash
# 单基因敲除
./scRNA_scTenifoldKnk.r -I 05.seu_annotated.rds -TG GeneX \
    -T 8 -p_val 0.05 -P 01. -O 14_scTenifoldKnk

# 多个基因要敲除，就换 -TG 跑多次（每次都换前缀，避免覆盖）
./scRNA_scTenifoldKnk.r -I 05.seu_annotated.rds -TG GeneY \
    -T 8 -p_val 0.05 -P 02. -O 14_scTenifoldKnk
```

**耗时**：建 10 个子网络 + 张量分解，通常是分钟到十几分钟级别，取决于细胞数和基因数（HVG 数量直接决定网络维度）。

### 输出

**表**

| 文件 | 内容 |
|---|---|
| `<前缀>diffRegulation_raw.csv` | **全部被检验基因**的结果表（`gene, distance, Z, FC, p.value, p.adj`）。**包含目标基因本身**（其 p 值为 0）|
| `<前缀>diffRegulation_sig.csv` | 只保留 `p.adj < -p_val` 的基因，**已剔除目标基因** |
| `<前缀>scTenifoldKnk_result.rds` | **完整分析结果对象**（`scTenifoldKnk()` 的返回值），可用来换 `q` 值重画网络图而不必重跑 |

**图（每张都出 PNG + PDF）**

| 文件 | 内容 |
|---|---|
| `<前缀>KO_effect_plot.png` | **敲除效应网络图**：中心是目标基因，向外连到受影响的基因；边的颜色区分正/负影响；节点颜色代表基因所属的富集通路，图例给出通路名与 FDR |
| `<前缀>diff_genes_barplot.png` | 按 `FC` 排序的 **top20 受影响基因**横向条形图 |
| `<前缀>Z_logp_scatter.png` | **Z 值 vs −log10(p.adj) 散点图**，标注高 Z 且显著的基因 |
| `<前缀>pie.png` | 显著 vs 不显著基因的**占比饼图** |

示例输出：`Images/scTenifoldKnk/`。

### 结果图怎么读

四张图分工不同，建议按 **① 饼图 → ② 效应网络 → ③ 条形图 → ④ 散点图** 的顺序看。

| # | 图 | 作用 | 怎么读 |
|---|---|---|---|
| 1 | `pie.png` | 先看「这次敲除到底有没有效果」 | 红色（显著）扇区**很窄 = 影响面小**，说明这个基因在网络里比较边缘；只有极少甚至 0 个显著基因时脚本会直接报错退出 |
| 2 | `KO_effect_plot.png` | **最重要的一张**，看「影响了谁、走的什么通路」 | 中心大圈 = 目标基因；放射到外圈 = 受影响基因。**节点颜色 = 该基因所属的富集通路**，一个节点有多种颜色说明它同时属于多条通路；图例里每条通路带 FDR |
| 3 | `diff_genes_barplot.png` | 看「影响最强的 20 个基因」 | 横轴 FC，越长受影响越强。**注意这是网络距离的变化倍数，不是表达量的 fold change** |
| 4 | `Z_logp_scatter.png` | 看整体分布与离群点 | 横轴 Z、纵轴 −log10(p.adj)；**蓝色竖虚线标 Z=2**，红色横虚线标 p 阈值。右下角弯上去的那条「弧」是绝大多数基因（Z 小、不显著）；真正值得注意的只有**越过两条线的点**，图上会标出基因名 |

**读 `KO_effect_plot` 的三个要点**：

1. **看中心基因连出去多少条边**——边多说明这个基因在网络里是枢纽，敲掉它影响面大
2. **看边的颜色**——红蓝分别代表敲除后上调/下调关联（脚本用 `plotKO(q = 0.99)` 只画最强的那些边，所以图上不是全部连接）
3. **看节点上的通路配色**——如果一个基因同时挂在多条通路颜色上（饼状节点），说明它被多条受影响通路共享，往往是关键的「交叉点」

### 注意事项

- **必须用原始 counts**：脚本取的是 `layer = "counts"`。如果对象里没有 counts 层（比如只保留了归一化数据），会直接失败。主流程 `scRNA_QC_step1.r` / `scRNA_Nor_harmony_step3.r` 的产物都带 counts。
- **必须先算过 HVG**：脚本用 `VariableFeatures(seu)` 选基因，没跑过 `FindVariableFeatures` 会得到空列表，只剩目标基因一个基因建网，结果无意义。
- **一次只能敲一个基因**：`-TG` 只接受单个基因名。多个基因要跑多次，且**每次换 `-P` 前缀**，否则产物会互相覆盖。
- **目标基因自身不出现在 `diffRegulation_sig.csv` 里**：原始表保留它（p 值恒为 0，因为是它自己），但下游分析前脚本会把它剔除。
- **没有显著基因时脚本会终止**：若 `p.adj < -p_val` 的基因数为 0，脚本 `stop()` 退出，连图都不会生成。遇到这种情况多半是这个基因在网络里太边缘，或 `-p_val` 设得太严。
- **`-p_val` 是「校正后 p 值」的阈值**，不是原始 p 值。
- **`Z_logp_scatter` 的标注条件跟随 `-p_val`**：图上标注的是 `abs(Z) > 2` 且 `p.adj < -p_val` 的基因，与红色横虚线（同一条阈值）完全一致，不会出现「标注点落在红线下方」的矛盾。
- **几个关键参数写死在代码里**：`nc_nNet = 10`（子网络数）、`nc_nCells = 500`（每网络细胞数）、`qc_mtThreshold = 0.1`（线粒体比例上限）、`qc_minLSize = 300`（每个细胞最少检出基因数）、`plotKO(q = 0.99)`（只画最强 99% 分位的边）。想改只能改代码。
- **细胞数不能太少**：每个子网络要抽 `nc_nCells = 500` 个细胞，`qc` 还会先过滤掉线粒体高、检出基因少的细胞。**细胞数不足 500（过滤后）时抽样会耗尽**，脚本可能报错或结果不稳。细胞类型特异的分析建议先把目标细胞类型单独 subset 出来再跑。
- **`plotKO` 会联网做富集**：脚本启动时设了 `setEnrichrSite("https://maayanlab.cloud/Enrichr")`，网络不通时这一步会失败或返回空富集。离线环境要提前考虑。
- **完整结果另存为 rds**：`<前缀>scTenifoldKnk_result.rds` 保存了 `scTenifoldKnk()` 的完整返回对象。想换 `plotKO()` 的 `q` 值重画网络图时，`readRDS` 后直接重画即可，不用再跑一遍十几分钟的分析。
- **`as.data.frame()` 会把矩阵稠密化**：第 119 行把稀疏 counts 转成 data.frame，基因数（HVG）越多内存占用越大。全库 2000 个 HVG 量级一般可接受。
- **FC 不是表达量倍数**：它是网络距离的变化倍数，数值可以很大（成百上千），**不要和 DE 的 log2FC 混着解释**。


---

## Scanpy 路线

本库有**两条并列的分析路线**：

| | Seurat 路线（第 1–28 章）| **Scanpy 路线（本章）** |
|---|---|---|
| 语言 | R | Python |
| 中间产物 | `.rds`（Seurat 对象）| `.h5ad`（AnnData 对象）|
| 脚本位置 | 库根目录 | **`scanpy/` 子目录** |
| 主流程 | `scRNA_QC_step1.r` → `dedoublets_step2.r` → `Nor_harmony_step3.r` → `clustering_step4.r` → `annotate_step5.r` | `scRNA_pipeline_scanpy.py`（一步做完 QC + 整合 + 聚类 + marker）→ `scRNA_annotate_scanpy.py` |

**两条路线产物格式不同，不能混用**：Seurat 路线的脚本读 rds，Scanpy 路线的脚本读 h5ad。需要跨路线时用 `scRNA_h5ad2rds.r` / `scRNA_h5ad_seurat_convert.r` 转换。

**本章脚本位于 `scanpy/` 子目录**，调用时路径要带上目录名：

```bash
./scanpy/scRNA_pipeline_scanpy.py -I raw.h5ad ...
```

### 路线总览

| 顺序 | 脚本 | 作用 | 输入 → 输出 |
|---|---|---|---|
| ① | `scRNA_pipeline_scanpy.py` | 质控、去双细胞、归一化、整合、聚类、找 marker（CPU）| h5ad → `preprocessed.h5ad` |
| ② | `scRNA_annotate_scanpy.py` | 画 marker dotplot、把 cluster 映射成细胞类型 | `preprocessed.h5ad` → `annotated.h5ad` |
| ③ | `scRNA_DE_scanpy.py` | 差异表达分析（Wilcoxon / t-test / logreg）+ 出图 | h5ad → `DE_analyzed.h5ad` |
| ④ | `scRNA_pipeline_scanpy_GPU.py` | ① 的 GPU 加速版（rapids-singlecell）| 同 ① |

①②③ 是**顺序依赖**的：① 的产物喂给 ②，② 的产物喂给 ③。④ 是 ① 的替代实现，**只跑其中之一**。

### 与 Seurat 主线的步骤对照

| Seurat 路线 | Scanpy 路线 | 说明 |
|---|---|---|
| `scRNA_QC_step1.r` | pipeline 的 Step 1–2 | 质控过滤 |
| `scRNA_dedoublets_step2.r` | pipeline 的 Step 3 | 去双细胞（scDblFinder ↔ **Scrublet**）|
| `scRNA_Nor_harmony_step3.r` | pipeline 的 Step 4–8 | 归一化 + HVG + PCA + Harmony |
| `scRNA_clustering_step4.r` | pipeline 的 Step 9–12 | 邻接图 + **Leiden** + UMAP + marker |
| `scRNA_annotate_step5.r` | `scRNA_annotate_scanpy.py` | 细胞类型注释 |

**区别**：Seurat 路线拆成 5 个 stepN 脚本、每步产物落盘可单独重跑；Scanpy 路线把前四步**合并在一个 pipeline 脚本**里，只有一个中间产物 `preprocessed.h5ad`。

### 环境与依赖

Scanpy 路线统一用 `scanpy` 环境。**四个脚本的 shebang 保持各自原样**：`scRNA_pipeline_scanpy.py` 与 `scRNA_pipeline_scanpy_GPU.py` 用 `#!/usr/bin/env python`（依赖 `PATH` 上有对应解释器），`scRNA_annotate_scanpy.py` 与 `scRNA_DE_scanpy.py` 用绝对路径 `/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python`。用 env 写法的两个脚本若直接 `./` 执行报解释器不对，显式调用即可：

```bash
/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python scanpy/scRNA_pipeline_scanpy.py -I raw.h5ad ...
```

核心依赖：`scanpy`、`anndata`、`pandas`、`numpy`、`matplotlib`、`harmonypy`（Harmony）、`igraph`（Leiden 用 `flavor='igraph'`，因此**不需要** `leidenalg`）。

额外注意：
- **Scrublet 步骤**（Step 3）依赖 `scrublet` 包。该包缺失时 pipeline 会在这一步失败——如果不需要去双细胞，只能自行注释掉这一步或先装包。
- **GPU 版**额外需要 `cupy` + `rapids-singlecell` + **CUDA 显卡**；三者缺一不可，纯 CPU 环境请用 CPU 版。

---

### 一、scRNA_pipeline_scanpy.py —— 质控、整合与聚类（CPU）

#### 功能

**一条命令走完单细胞数据的前四步预处理**：质控 → 去双细胞 → 归一化 → 高变基因 → 标准化 → PCA → Harmony 去批次 → 邻接图 → Leiden 聚类 → UMAP → 找 marker，最后存成 `preprocessed.h5ad` 供下游注释使用。

对应 Seurat 路线的 `step1`–`step4`。

#### 输入

一个 **h5ad 文件**（AnnData）。`adata.X` 应为原始 counts（整数）；`adata.obs` 若含批次列可交给 `-BK` 用于去批次与分步去双细胞。

#### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| 输入 h5ad 路径 |
| `-O, --outpath` | 当前目录 | 输出目录 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-minG, --min_genes` | **None** | 每个细胞最少基因数 |
| `-maxG, --max_genes` | **None** | 每个细胞最多基因数 |
| `-minU, --min_counts` | **None** | 每个细胞最少 UMI 数 |
| `-maxU, --max_counts` | **None** | 每个细胞最多 UMI 数 |
| `-minC, --min_cells` | **None** | 每个**基因**至少在多少个细胞中表达才保留（这是一条基因层面的阈值）|
| `-maxMT, --max_mt` | **None** | 线粒体基因占比上限（%）|
| `-maxHB, --max_hb` | **None** | 血红蛋白基因占比上限（%）|
| `-BK, --batch_key` | None | 批次列名（`adata.obs` 的列）；给了就去批次 + 分步去双细胞，并在 UMAP 上单独出一张图 |
| `-nPC, --n_pcs` | 30 | 聚类用的 PC 数（注意 PCA 实际算了 50 个主成分，这里只取前 nPC 个）|
| `-nHVG, --n_top_genes` | 2000 | 高变基因数量 |
| `-res, --resolution` | 1.0 | 最终 Leiden 聚类分辨率 |
| `-dpi, --dpi` | 300 | 出图 dpi |

#### 用法示例

**必须分两次跑**（因为 7 个 QC 阈值默认全为 None）：

```bash
# 第一次：只给输入输出，脚本会画完 QC 图后退出
./scanpy/scRNA_pipeline_scanpy.py -I raw.h5ad -BK sample -O 01_scanpy -P 01.

# 看完 before_QC_violin.png / before_QC_scatter.png 定阈值，第二次带上全部 7 个参数
./scanpy/scRNA_pipeline_scanpy.py -I raw.h5ad -BK sample -O 01_scanpy -P 01. \
    -minG 200 -maxG 2500 -minU 500 -maxU 10000 -minC 3 -maxMT 5 -maxHB 5
```

#### 输出

| 文件 | 说明 |
|---|---|
| `<前缀>before_QC_violin.png/.pdf` | 过滤前 4 个 QC 指标的小提琴图（有 `-BK` 时按批次分组）|
| `<前缀>before_QC_scatter.png/.pdf` | 过滤前 3 张散点图（基因数 vs UMI、mt vs UMI、hb vs UMI）|
| `<前缀>after_QC_violin.png/.pdf` | 过滤后同上 |
| `<前缀>after_QC_scatter.png/.pdf` | 过滤后同上 |
| `<前缀>HVGs.png/.pdf` | 高变基因 | 
| `<前缀>pca_loadings.png/.pdf` | PCA 载荷 |
| `<前缀>pca_variance_ratio.png/.pdf` | 各主成分方差贡献 |
| `<前缀>pca.png/.pdf` | PCA 散点 |
| `<前缀>umap_clusters.png/.pdf` | **6 个分辨率并排**的 UMAP（0.2 / 0.5 / 0.8 / 1.0 / 1.2 / 1.5）|
| `<前缀>umap_leiden.png/.pdf` | 最终分辨率的 Leiden 聚类 UMAP（标签画在点上）|
| `<前缀>umap_<batch_key>.png/.pdf` | 按批次着色的 UMAP（仅 `-BK` 时）|
| `<前缀>marker_genes.png/.pdf` | 各 cluster top20 marker 图 |
| `<前缀>all_markers.csv` | 全部 cluster 的 marker 表（cluster / 基因 / logFC / p 值等）|
| `<前缀>preprocessed.h5ad` | **最终对象**（含降维、聚类、marker 结果）|

示例输出：`Images/scanpy/pipeline/`。GPU 版产物结构相同。

#### 注意事项

- **两段式 QC 是设计，不是 bug**：7 个 QC 阈值任一为 None 时，脚本画完 `before_QC_*` 图后**以退出码 0 结束**，并打印缺少哪些参数 + 一条示例命令。补齐后重跑才会继续。所以**第一次跑一定会停在这里**。
- **`-minC` 与其他 6 个参数的方向不同**：`-minC` 是「基因至少在几个细胞里表达」（基因过滤），其余 6 个是「细胞要满足什么条件」（细胞过滤）。
- **HVG 用的是 `counts` 层**（`layer="counts"`），不是 log 化后的数据；给了 `-BK` 会按批次选高变基因。
- **Leiden 用 `flavor='igraph'`**，所以不需要 `leidenalg` 包。
- **会跑 6 个分辨率**（0.2–1.5）用于对比，`-res` 只决定最终保存哪一个。数据大时这 6 次聚类是主要耗时之一。
- **输出 h5ad 里三个矩阵分别是**：`adata.X` = **z-score**（`sc.pp.scale(max_value=10)` 就地覆盖，取值约 −3 ~ 10）；`adata.layers['counts']` = 原始 counts；`adata.layers['log1p']` = log 归一化表达量；`adata.raw.X` = 全基因的 log1p。**想在 `.X` 上画表达量或算均值时，请显式改用 `layers['log1p']` 或 `raw`**，否则拿到的是 z-score。
- 该 h5ad 会较大：`X` 是稠密的 float64 z-score 矩阵（未压缩时约为 `细胞数 × 基因数 × 8` 字节），加上 counts / log1p / raw 三份表达数据。细胞数很多时需要相应内存；若机器内存吃紧，可自行把 `sc.pp.scale` 改为只对 HVG 做。
- **`-BK` 写错列名在读入数据后立即退出**：脚本在 Step 1 之后校验该列是否存在，不存在则打印可用列名并 `sys.exit(1)`，不会白跑 QC。
- **PCA 算了 50 个主成分，但聚类只用前 `-nPC` 个**（默认 30）——`-nPC` 超过 50 无意义。
- **Scrublet 依赖 `scrublet` 包**，缺失会在 Step 3 失败。

---

### 二、scRNA_annotate_scanpy.py —— 细胞类型注释

#### 功能

**把聚类的编号变成有生物学含义的细胞类型名**。分两个阶段：
1. **注释前**（`-MB`）：按给定的「细胞类型 → marker 基因」表画 dotplot，用来判断每个 cluster 是什么
2. **注释后**（`-MA`）：按给定的「cluster → 细胞类型」表把标签写进 `adata.obs['celltype']`，再出 UMAP 与 dotplot

对应 Seurat 路线的 `scRNA_annotate_step5.r`。

#### 输入

- `-I`：一个已聚类、已算 UMAP 的 h5ad（通常是 pipeline 的 `preprocessed.h5ad`）
- `-MB`：**注释前的 marker 表**，2 列 TSV，无表头：`CellType` + `marker基因（逗号分隔）`
- `-MA`：**注释后的映射表**，3 列 TSV，无表头：`Cluster` + `marker基因（逗号分隔）` + `CellType`

#### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| 已聚类的 h5ad |
| `-MB, --markers_before` | None | 注释前 marker 表（2 列 TSV）|
| `-MA, --markers_after` | None | 注释后映射表（3 列 TSV）|
| `-CK, --cluster_key` | leiden | cluster 标签所在的 `adata.obs` 列 |
| `-O, --outpath` | 当前目录 | 输出目录 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-S, --split` | 关 | 按细胞类型拆开，每种类型单独出一张 dotplot（marker 很多时用）|
| `-CO, --celltype_order` | None | 指定细胞类型在图上的顺序（逗号或空格分隔）|
| `--dendrogram` / `--no-dendrogram` | 显示树状图 | 控制注释前 dotplot 是否画树状图。默认显示，用 `--no-dendrogram` 关闭（该开关无短参数）|
| `-VR, --var_group_rotation` | 45 | 细胞类型标签旋转角度 |
| `-DW, --dotplot_width` | None（自动）| dotplot 宽度（英寸）|
| `-DH, --dotplot_height` | None（自动）| dotplot 高度（英寸）|
| `-dpi, --dpi` | 300 | 出图 dpi |

#### 用法示例

```bash
# 阶段一：画 marker dotplot 判断 cluster 身份
./scanpy/scRNA_annotate_scanpy.py -I 01_scanpy/01.preprocessed.h5ad \
    -MB markers.tsv -P 02. -O 02_annotation

# 阶段二：写回注释并出图
./scanpy/scRNA_annotate_scanpy.py -I 01_scanpy/01.preprocessed.h5ad \
    -MA annotation.tsv -CO 'CellTypeA,CellTypeB,CellTypeC' -P 02. -O 02_annotation
```

#### 输出

**`-MB` 模式**

| 文件 | 说明 |
|---|---|
| `<前缀><cluster_key>_umap.png/.pdf` | cluster 的 UMAP |
| `<前缀><cluster_key>_dotplot.png/.pdf` | 各 cluster 的 marker dotplot（带树状图）|
| `<前缀>markers_dotplot/` | 仅 `-S`：该子目录下每种细胞类型一张 dotplot |

**`-MA` 模式**

| 文件 | 说明 |
|---|---|
| `<前缀>umap_celltype.png/.pdf` | 按 celltype 着色的 UMAP（标签画在点上）|
| `<前缀>dotplot_celltype.png/.pdf` | 各 celltype 的 marker dotplot |
| `<前缀>annotated.h5ad` | **最终对象**（`adata.obs['celltype']` 已写入）|

示例输出：`Images/scanpy/annotate/`。

#### 注意事项

- **`-MB` 与 `-MA` 是两阶段，不是二选一**：先用 `-MB` 看图，再用 `-MA` 写回。两个都给时按顺序都会执行。
- **两个表都是 TSV 且无表头**，列顺序固定；`-MB` 是 2 列、`-MA` 是 3 列，给错列数会直接报错。
- **`-MA` 里没有出现的 cluster 会被填成 `"Unknown"`**，脚本会打印注释覆盖率（`Annotated cells: x / y`）——覆盖率明显不足说明映射表漏了 cluster。
- **关闭树状图用 `--no-dendrogram`**：默认会画树状图（`sc.tl.dendrogram` 预计算），marker 很多时可以让 dotplot 更紧凑。
- **`-MA` 输出的 dotplot 不带树状图**（代码里写死 `dendrogram=False`），只有 `-MB` 阶段的带。
- **`-CO` 只影响图的顺序**，不改变注释结果。
- **输出文件名里的 cluster 名会做字符替换**（去掉特殊字符），带点号或斜杠的 cluster 名要注意。

---

### 三、scRNA_DE_scanpy.py —— 差异表达分析

#### 功能

**找出每个 cluster 的 marker 基因**：对指定 cluster 列做组间差异检验（Wilcoxon / t-test / logreg），出 marker 概览图与逐 cluster 的 dotplot，并把结果存表。

与本库 R 侧脚本的区别：`scRNA_pseudobulk_DE.r` 是**样本级** pseudobulk + DESeq2（可外推的组间结论）；本脚本是**细胞级** `rank_genes_groups`（marker 发现用途，有伪重复问题）。两者用途不同，不要互相替代。

#### 输入

一个已聚类的 h5ad（`adata.obs` 里有 cluster 列，且已归一化 + log 化）。

#### 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-I, --input` | 无（必填）| 已聚类的 h5ad |
| `-CK, --cluster_key` | leiden | cluster 标签所在的列 |
| `-M, --method` | wilcoxon | 检验方法，只接受 `wilcoxon` / `t-test` / `logreg`；写错由 argparse 直接拒绝并列出可选值 |
| `-n_genes` | 20 | marker 概览图里每个 cluster 显示的基因数 |
| `-top_n` | 10 | 逐 cluster dotplot 里取的 top marker 数 |
| `-P, --prefix` | 空 | 输出文件名前缀 |
| `-O, --outpath` | 当前目录 | 输出目录 |
| `-dpi, --dpi` | 300 | 出图 dpi |
| `-F, --force` | 关 | 强制重算 DE（默认复用已有结果）|

#### 用法示例

```bash
./scanpy/scRNA_DE_scanpy.py -I 02_annotation/02.annotated.h5ad \
    -CK celltype -M wilcoxon -n_genes 20 -top_n 10 -P 03. -O 03_DE
```

#### 输出

| 文件 | 说明 |
|---|---|
| `<前缀>marker_genes.png/.pdf` | 各 cluster top `-n_genes` 个 marker 的图 |
| `<前缀>clusters_DE/` | 子目录 |
| `<前缀>clusters_DE/<前缀>all_markers.csv` | 全部 cluster 的 marker 表 |
| `<前缀>clusters_DE/<前缀>cluster_<id>_dotplot.png/.pdf` | 每个 cluster 的 top `-top_n` marker dotplot |
| `<前缀>DE_analyzed.h5ad` | 含 `adata.uns['rank_genes_groups']` 的对象 |

示例输出：`Images/scanpy/DE/`。

#### 注意事项

- **默认会复用已有结果**：若 `adata.uns['rank_genes_groups']` 已存在且 `groupby` 与当前 `-CK` 一致，跳过计算直接出图。改 `-CK` 会自动重算；想强制重算加 `-F`。
- **`-M` 有取值校验**：只接受 `wilcoxon` / `t-test` / `logreg`，写错会在参数解析阶段直接报 `invalid choice` 并列出可选值。
- **结果是细胞级检验**：p 值会被细胞数放大，**不能当作组间差异的统计证据**——那是 `scRNA_pseudobulk_DE.r` 的职责。
- **输出子目录名含前缀**：`<前缀>clusters_DE/`，注意别与其它输出混淆。
- **cluster 名里的特殊字符会被替换**后再用于文件名。

---

### 四、scRNA_pipeline_scanpy_GPU.py —— GPU 加速版

#### 功能

**与 `scRNA_pipeline_scanpy.py` 完全相同的 13 步流程**，只是把 `scanpy`（`sc.*`）换成 `rapids-singlecell`（`rsc.*`），在 CUDA 显卡上运行以缩短耗时。

**两者只跑一个**，不要串联。输出格式一致，下游脚本通用。

#### 输入

同 CPU 版：一个 h5ad（`adata.X` 为原始 counts）。

#### 参数

与 CPU 版**完全相同**（`-I -O -P -minG -maxG -minU -maxU -minC -maxMT -maxHB -BK -nPC -nHVG -res -dpi`），**额外多一个**：

| 参数 | 默认值 | 说明 |
|---|---|---|
| `-raw, --save_raw` | 关 | 输出 h5ad 时保留**全部基因**；默认只保留高变基因（可显著减小文件）|

#### 用法示例

```bash
./scanpy/scRNA_pipeline_scanpy_GPU.py -I raw.h5ad -BK sample -O 01_scanpy_gpu -P 01. \
    -minG 200 -maxG 2500 -minU 500 -maxU 10000 -minC 3 -maxMT 5 -maxHB 5
```

#### 输出

与 CPU 版一致（`before_QC_*` / `after_QC_*` / `HVGs` / `pca*` / `umap_clusters` / `umap_leiden` / `umap_<batch>` / `marker_genes` / `all_markers.csv` / `preprocessed.h5ad`）。

#### 注意事项

- **必须有 CUDA 显卡**：脚本启动时先 `import cupy` 并检查 GPU，缺 `cupy` 或没有可用 GPU 会**立即报错退出**，不会退回 CPU 模式。
- **依赖 `cupy` + `rapids-singlecell`**，还需与 CUDA 版本匹配。
- **`-raw` 决定 h5ad 里存多少基因**：默认只存高变基因（文件小、但下游只能用这批基因）；需要全部基因时加 `-raw`。
- **数值结果与 CPU 版可能有细微差异**：GPU 与 CPU 的浮点运算顺序不同，聚类边界附近的细胞可能落到不同 cluster，属于正常现象。
- **无 GPU 环境下无法验证**：本脚本在本机未做过运行验证，仅通过语法检查。


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
- velocity 输入准备（`scRNA_velocityPrepare.r`）：`anndataR`（AnnData 转换 + `write_h5ad`）、`Seurat`
- 虚拟敲除（`scRNA_scTenifoldKnk.r`）：`scTenifoldKnk`（GitHub `cailab-tamu/scTenifoldKnk`）、`Seurat`、`tidyverse`、`patchwork`、`igraph`、`ggrepel`、`enrichR`（`plotKO` 会联网访问 Enrichr）
- 共表达网络（`scRNA_hdWGCNA.r`）：`hdWGCNA`（GitHub `smorabit/hdWGCNA`）、`WGCNA`、`UCell`、`ggradar`、`Seurat`、`harmonypy`（metacell 与降维依赖 harmony）、`tidyverse`、`patchwork`
- pySCENIC 输入准备与 RSS（`scRNA_pyscenic_prepare.r` / `scRNA_pyscenic_RSS.r`）：`SCopeLoomR`（写/读 loom）、`SCENIC`（`calcRSS` / `plotRSS` / `getAUC`）、`gridExtra`、`Seurat`、`tidyverse`
- 亚型聚类（`scRNA_subcelltype_clustering.r`）：`clustree`、`harmony`、`reticulate`（Leiden 时需 `leidenalg`）、`patchwork`
- 差异表达火山图与 pseudobulk DE（`scRNA_pseudobulk_DE.r`）：`ggrepel`（标签避让）、`RColorBrewer`、`scales`；细胞类型 > 66 时另需 `viridis`
- pseudobulk 差异表达（`scRNA_pseudobulk_DE.r`）：`DESeq2`（Bioconductor）、`Seurat`、`Matrix`、`tidyverse`、`ggrepel`、`RColorBrewer`、`scales`
- 细胞比例分析（`scRNA_spRNA_cellPro.r`）：`ggalluvial`（桑基图）、`ggpubr`（显著性标注）、`patchwork`（亚型图拼接）、`tidyverse`、`Seurat`、`RColorBrewer`、`scales`
- 自动注释（`scRNA_SingleR.r`）：`SingleR`、`SingleCellExperiment`、`scater`（参考对象的 logNormCounts）、`celldex`（内置 human/mouse 参考）、`Seurat`、`tidyverse`
- 通路富集（`scRNA_spRNA_enrichment.r` + 配套 `enrichment_analysis.r`）：主脚本只准备输入文件，实际富集由配套脚本执行（需要 `clusterProfiler`、`enrichplot`、`GseaVis`、`org.Hs.eg.db`/`org.Mm.eg.db`）
- 通路活性打分（`scRNA_GSVA_enrichment.r`）：`GSVA`、`msigdbr`（GO:BP / KEGG MEDICUS / REACTOME 基因集）、`limma`（组间比较）、`ComplexHeatmap`（热图）、`Seurat`、`tidyverse`
- 自定义基因集打分（`scRNA_genesetScore.r`）：`GSVA`（GSVA + ssGSEA）、`AUCell`、`Seurat`（AddModuleScore）、`tidyverse`
- 去双细胞（`scRNA_dedoublets_step2.r` / `scRNA_DoubleFinder.r`）：`scDblFinder` / `DoubletFinder`、`celda`（DecontX）、`BiocParallel`、`Seurat`、`tidyverse`、`patchwork`、`reticulate`（Leiden 时）
- 判别分析（`scRNA_ROC.r`）：`pROC`（AUC / ROC）、`tidyverse`、`Seurat`
- 细胞通讯（`scRNA_cellchat_single.r` / `scRNA_cellchat_compare.r`）：`CellChat`、`NMF`（通讯模式聚类）、`ggalluvial`、`patchwork`、`ComplexHeatmap`（outgoing 通路热图）、`Seurat`、`tidyverse`
- 拟时序与轨迹（`scRNA_spRNA_Monocle3.r`）：`monocle3`、`SeuratWrappers`（`as.cell_data_set`）、`pheatmap`（模块热图）、`igraph`、`leidenbase`、`Rcpp`、`Seurat`、`tidyverse`

**Python 侧**（`scRNA_spRNA_h5ad_preprocess.py`、`scRNA_composition_sccoda.py`、`scRNA_velocity.py`、`scRNA_pyscenic.py`）：

- `scCODA`、`tensorflow` / `tensorflow-probability`
- `anndata`（注意与 `zarr` / `numpy` 的版本匹配：`anndata < 0.12` 拒 zarr ≥3，`anndata ≥ 0.13` 要求 numpy ≥2.1）
- h5ad 预处理：`scanpy`、`anndata`、`scipy`、`numpy`、`pandas`
- 组成分析：`numpy`、`pandas`、`matplotlib`、`seaborn`
- RNA velocity（`scRNA_velocity.py`）：`scanpy`、`scvelo`、`anndata`、`numpy`、`pandas`、`h5py`、`matplotlib`（建议用 `scanpy` 环境运行）
- pySCENIC（`scRNA_pyscenic.py`）：`pyscenic`（0.12.x，含 `arboreto` / `ctxcore` / `primp`）、`loompy`、`dask`、`pandas`、`numpy`（须用独立的 pySCENIC 环境，Python 3.7）
- 桑基图需 `Rscript` 与 R 的 `ggalluvial`（缺失时该图跳过，其余产物不受影响）

**解释器**：脚本 shebang 用 `/usr/bin/env` 解析 `Rscript` / `python3`，运行前请确保所用环境在 `PATH` 上（例如 `export PATH=/path/to/env/bin:$PATH`）；若不在 PATH 上，直接 `./<script>` 会报解释器未找到，此时用 `Rscript <script>.r …` / `/path/to/env/bin/python <script>.py …` 显式调用即可。

R 脚本启动时会检查缺失的依赖并提示安装命令；Python 侧缺包会在导入处报错。
