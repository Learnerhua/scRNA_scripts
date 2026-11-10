# scRNA-seq Analysis Scripts

A general-purpose, **continuously growing** toolkit for single-cell RNA-seq (scRNA-seq) analysis. It combines a step-numbered **Seurat v5 main pipeline** — per-sample QC and merging, doublet removal, normalization and batch integration, clustering and differential expression, cell-type annotation — with **standalone scripts** for data preparation and format conversion (h5ad count cleaning, h5ad → Seurat rds, h5ad ⇄ rds) and for downstream analyses (target-gene expression visualization, per-sample × per-cell-type statistics, compositional differential abundance, …), plus the companion R helper that renders the composition sankey. Every script takes its inputs as command-line paths (10x directories, `.h5`, `rds`, `h5ad`, tables, …) and writes its results to files (data objects, tables, figures, …), so any stage can be run on its own or chained into a full pipeline. **New scripts are added over time.**

[中文说明 (Chinese)](README_zh.md) · Full manual (Chinese): [`USAGE.md`](USAGE.md)

## Scripts

| Script | Type | Main function | Input → Output | Example outputs |
|---|---|---|---|---|
| `scRNA_spRNA_h5ad_preprocess.py` | standalone (data prep) | Makes sure the matrix is raw integer counts (uses `X`, falls back to `raw.X`, aborts if neither looks like counts), renames genes from a chosen `var` column and sums duplicate gene names; writes both `X` and `layers["counts"]` | `.h5ad` (directory / single file) → `.h5ad` with `X` + `layers["counts"]` | — |
| `scRNA_h5ad2rds.r` | standalone (conversion) | Batch `.h5ad` → Seurat `.rds`; two backends (`anndataR`, pure R, default; `sceasy` via reticulate). The input must carry `layers["counts"]` | `.h5ad` (directory / single file) → `.rds` | — |
| `scRNA_h5ad_seurat_convert.r` | standalone (conversion) | Single-file bidirectional conversion, direction chosen from the extension: `.h5ad` → `.rds` and `.rds` → `.h5ad` (factors → character first; assay / layer / compression selectable) | `.h5ad` ⇄ `.rds` | — |
| `scRNA_mtx2rds.r` | standalone (data prep) | Builds a Seurat object from **Broad Single Cell Portal** exports: dense matrix (`-E`) or `.mtx` + cell/gene files (`-M/-C/-G`), plus a required metadata table (`-m`, cell-name column defaults to `NAME`, a leading `TYPE` row is dropped). Genes get `_`→`-` and `make.unique`; non-integer input keeps only the `data` layer | matrix + metadata → `seurat.rds` | — |
| `scRNA_subcelltype_clustering.r` | standalone (refinement) | Sub-clusters every cell type of an annotated object separately: per-cell-type re-embedding (normalize → HVG → PCA → Harmony) plus a multi-resolution clustree, then re-clusters at the resolution you recorded per cell type, runs sub-celltype DE (markers, top-10, dotplots) and writes the `sub_celltype` labels back to the original object | annotated rds + resolution table → sub-celltype rds + per-subtype DE | `Images/subcelltype_clustering/` |
| `scRNA_pseudobulk_DE.r` | standalone (statistics) | **Finds the genes that differ between two groups within each cell type.** Raw counts are summed per (cell type × sample) and `-G1` vs `-G2` is tested at sample level with DESeq2 — one independent test per cell type, so the result reads as "in this cell type, these genes differ". Units with too few cells and cell types with too few samples per group are skipped. Outputs a long DE table plus a unit provenance table, and a volcano plot | annotated rds (counts layer) → DE csv + volcano plot | `Images/pseudobulk_DE/` |
| `scRNA_QC_step1.r` | pipeline `step1` | Per-sample QC with sample-specific thresholds (quantiles of nFeature_RNA / nCount_RNA + hard cap on percent.mt, optional percent.HB for red-blood-cell removal), then merge QC-passing samples; gene-level `min.cells` filtering is applied **after** merging | 10x directory / `.h5` list (CSV) → merged Seurat object (counts only) + QC violin/scatter plots | `Images/01_QC/` |
| `scRNA_dedoublets_step2.r` | pipeline `step2` | Per-sample doublet detection with **scDblFinder**; optional **DecontX** ambient-RNA decontamination (also per sample) | merged rds → singlet-only rds (`seu_clean.rds`) | `Images/02_deContamination/` |
| `scRNA_Nor_harmony_step3.r` | pipeline `step3` | NormalizeData (LogNormalize) → FindVariableFeatures (vst) → ScaleData → RunPCA → **Harmony** batch correction | clean rds → normalized + integrated rds (`seu_harmony.rds`) | `Images/03_Nor_Harmony/` |
| `scRNA_clustering_step4.r` | pipeline `step4` | Two stages: multi-resolution clustree exploration (omit `-Res`) → clustering + UMAP + `FindAllMarkers` DE (pass `-Res`); Louvain or Leiden | integrated rds → clustered rds + `clusters_DE/` marker tables and plots | `Images/04_clustering/`, `Images/04_clustering_louvain/` |
| `scRNA_annotate_step5.r` | pipeline `step5` | Cell-type annotation, two exclusive modes: `-MB` = marker-list exploration dotplots (writes no rds); `-MA` = cluster → cell-type mapping (writes `seu_annotated.rds`, plus celltype UMAP / dotplot) | clustered rds + marker TSV → annotated rds + plots | `Images/05_annotation/` |
| `scRNA_gene_expr.r` | standalone | Average expression of target genes in every cell type (CSV) plus 5 plot types: Ridge / Violin / Feature / Dot / Heatmap | annotated rds → expression CSV + plots | `Images/gene_expr/` |
| `scRNA_gene_DE.r` | standalone | For one target gene: sample × cell-type mean-expression table, per-celltype z-scored heatmap, one barplot per cell type, and sample-level boxplots with Wilcoxon test + significance bars (each point is a sample, not a cell) | annotated rds → CSV + heatmap / bar / box plots | `Images/gene_DE/` |
| `scRNA_composition_sccoda.py` | standalone | Compositional differential-abundance analysis with **scCODA** (Dirichlet-multinomial Bayesian model); `--deterministic` gives bit-reproducible MCMC | `.h5ad` / sample × cell-type count table / per-cell obs table → effect tables + sankey / boxplot / dispersion / forest plots | `Images/composition_sccoda/` |
| `scRNA_spRNA_cellPro.r` | standalone (statistics) | **How much of each cell type is in each sample?** Counts and proportions per sample and per group, a 100 % stacked bar chart, a group→cell-type sankey plot, and a per-cell-type comparison of the two groups (Wilcoxon, one point per sample). Optionally repeats the comparison for sub-cell types | annotated rds → proportion tables + bar / sankey / comparison plots | `Images/cellPro/` |
| `scRNA_SingleR.r` | standalone (annotation) | **Assigns cell-type names automatically from a reference dataset.** Matches every cluster of the query object against a reference (celldex built-in human/mouse via `-S`, or your own Seurat/RData reference via `-RDS`/`-RD`) with SingleR, writes the cluster→cell-type mapping back into the object, and draws the score heatmap plus a UMAP coloured by the inferred labels | clustered rds (with `seurat_clusters`) + reference → mapping csv + annotated rds + figures | `Images/SingleR/` |
| `scRNA_spRNA_enrichment.r` | standalone (pathways) | **Turns differential genes into pathway-level conclusions.** Per cell type it runs GO/KEGG enrichment (using the genes detected in that cell type as background) and GSEA (all genes ranked by log2FC), reporting pathway p/q values, NES, dotplots and barplots. Prepares the inputs and calls the external `enrichment_analysis` command | annotated rds + DE result csv → enrichment tables + figures | `Images/enrichment/` |
| `scRNA_GSVA_enrichment.r` | standalone (pathways) | **Compares pathway activity between two groups without needing a DE result.** Aggregates counts per (sample × cell type), normalises to log2(CPM+1), scores every pathway with GSVA (GO:BP, KEGG MEDICUS and REACTOME gene sets), then compares the scores between the two groups with limma and draws a per-cell-type pathway heatmap (control on the left) | annotated rds → pathway score table + pathway DE table + heatmaps | `Images/GSVA_enrichment/` |
| `scRNA_genesetScore.r` | standalone (pathways) | **Scores a gene set of your own choosing.** Reads a text file of gene symbols and scores it four ways — GSVA and ssGSEA at the cell-type level, AddModuleScore and AUCell at single-cell level — then plots the cell-type barplots and the single-cell UMAP features. Useful for asking which cell types express a gene set of interest | annotated rds + gene set file → score table + 4 kinds of figures | `Images/genesetScore/` |
| `enrichment_analysis.r` | helper (**required** by the script above) | Runs the actual GO/KEGG enrichment and GSEA with clusterProfiler: takes a target-gene list (`-G`) plus an optional background set (`-E`), or a ranked all-gene table for GSEA (`-GSEA` + `-GM`). Writes the enrichment tables, dotplots and barplots. Also usable on its own for bulk RNA | gene list / ranked table → enrichment tables + figures | — |
| `scRNA_ROC.r` | standalone (statistics) | **How well can a set of target genes tell two groups apart?** Scores each gene (or, with `-GS`, one combined gene-set score) per sample, computes the AUC of each gene as a classifier, draws one ROC curve per gene, and writes the AUCs ranked in a summary table. This is a discrimination analysis, not a significance test — pair it with `scRNA_gene_DE.r` or `scRNA_pseudobulk_DE.r` when you need p values | annotated rds + group csv + gene list → AUC table + ROC curves | `Images/ROC/` |
| `scRNA_cellchat_single.r` | standalone (cell communication) | **Which cell types talk to each other, and through which ligands and receptors?** Runs the CellChat workflow on one annotated object: builds the interaction network, lists the significant ligand–receptor pairs and signalling pathways, draws the overall and per-sender circle plots, the top-9 pathway chord diagrams, and clusters the cell groups into outgoing communication patterns. Needs two runs — first without `-NP` to get the `selectK` plot, then again with the number of patterns you read off it | annotated rds → communication tables + network / pathway / pattern figures | `Images/cellchat/single/` |
| `scRNA_cellchat_compare.r` | standalone (cell communication) | **How does the communication network differ between two groups?** Merges two CellChat objects (produced by the script above) and compares interaction counts, interaction strength, pathway composition and the ligand–receptor pairs that differ significantly between the groups | two CellChat rds → merged object + comparison figures | `Images/cellchat/compare/` |
| `scRNA_DoubleFinder.r` | standalone (alternative) | **Alternative doublet-removal implementation.** Same job as the step-2 script but using DoubletFinder instead of scDblFinder, and modelling the whole merged object at once rather than sample by sample. The expected doublet rate must be given (`-DR`) instead of being estimated, and the output keeps all layers. Use it to cross-check the step-2 result with an independent algorithm — never run both on the same data | merged rds → clean rds | `Images/DoubleFinder/` |
| `sankey_composition.r` | helper (**required** by the script above) | Renders the composition sankey (alluvial) plot with `ggalluvial`. Called automatically by `scRNA_composition_sccoda.py`; **must live in the same directory** as it | composition CSV → sankey PNG/PDF | ↑ |

## Pipeline overview

```text
optional upstream preparation — when starting from .h5ad instead of 10x outputs

  raw .h5ad ──▶ scRNA_spRNA_h5ad_preprocess.py ──▶ scRNA_h5ad2rds.r ──▶ Seurat rds
                 (integer counts + layers["counts"])                    │
                                                                        ▼ (feeds step1)

main pipeline — scripts carry the stepN prefix and are chained through rds files

        │  step1 scRNA_QC_step1.r           per-sample QC + merge
        ▼
        │  step2 scRNA_dedoublets_step2.r   doublets (scDblFinder) [+ DecontX]
        ▼
        │  step3 scRNA_Nor_harmony_step3.r  normalize + HVG + PCA + Harmony
        ▼
        │  step4 scRNA_clustering_step4.r   clustree → cluster + UMAP + DE
        ▼
        │  step5 scRNA_annotate_step5.r     cluster → cell type
  seu_annotated.rds ──┬──▶ scRNA_gene_expr.r   target genes across cell types
                      └──▶ scRNA_gene_DE.r     per-sample × per-celltype expression

standalone — independent of the pipeline numbering
  scRNA_composition_sccoda.py   .h5ad / count table → compositional DA (+ sankey_composition.r)
```

## Highlights

- **Per-sample adaptive QC** — thresholds come from each sample's own quantiles, so samples with very different sequencing depth stay comparable.
- **Normalization happens once, at step3** — genes are filtered before normalization, avoiding inconsistent denominators.
- **Re-plot-only caching** — if a key output already exists, most scripts skip the heavy computation and only redraw the figures (change plot size in seconds; delete the rds to force a rerun).
- **Statistics at sample level, not cell level** — boxplot comparisons in `scRNA_gene_DE.r` use one point per sample, avoiding pseudoreplication. The test is two-sided by default (`--alt` switches to one-sided) and the stats table carries both `p_value` and `p_adj`, where `p_adj` is BH-corrected across the cell types of that gene. Note the hard floor of the two-sided rank-sum test: with 3 samples per group the smallest attainable p is 0.1, so nothing can reach significance.
- **Reproducible scCODA runs** — `--deterministic` makes repeated runs bit-identical.
- **Faithful plots** — the composition sankey is rendered by `ggalluvial` through R instead of being re-implemented in Python.

## Requirements

**R side** (the pipeline steps and the standalone R scripts)

- core: `Seurat` (v5), `tidyverse` (dplyr / ggplot2 / tidyr …), `patchwork`, `argparse`
- doublets: `scDblFinder`, `BiocParallel`; optional DecontX: `celda`
- integration / exploration: `harmony`, `clustree`
- optional Leiden clustering: `reticulate` (the scripts call `py_require("leidenalg")`)
- `scRNA_gene_DE.r` also uses `reshape2`
- sankey (`sankey_composition.r`): `ggalluvial`, `RColorBrewer`, `scales`, `ggplot2`, `dplyr`, `tidyr`
- format conversion (`scRNA_h5ad2rds.r`, `scRNA_h5ad_seurat_convert.r`): `anndataR` (default backend); optional `sceasy` (via reticulate, whose Python environment needs `numpy` + `anndata`)
- matrix import (`scRNA_mtx2rds.r`): `Matrix` (`readMM`), `data.table` (`fread`)
- sub-celltype clustering (`scRNA_subcelltype_clustering.r`): `clustree`, `harmony`, `reticulate` (`leidenalg` for Leiden), `patchwork`
- pseudobulk DE (`scRNA_pseudobulk_DE.r`): `DESeq2` (Bioconductor), `Seurat`, `Matrix`, `tidyverse`, `ggrepel`, `RColorBrewer`, `scales`; `viridis` when there are more than 66 cell types
- cell proportion analysis (`scRNA_spRNA_cellPro.r`): `ggalluvial` (sankey), `ggpubr` (significance labels), `patchwork` (sub-celltype panels), `tidyverse`, `Seurat`, `RColorBrewer`, `scales`
- discrimination analysis (`scRNA_ROC.r`): `pROC` (AUC / ROC), `tidyverse`, `Seurat`
- cell-cell communication (`scRNA_cellchat_single.r` / `scRNA_cellchat_compare.r`): `CellChat`, `NMF` (communication patterns), `ggalluvial`, `patchwork`, `ComplexHeatmap` (outgoing pathway heatmap), `Seurat`, `tidyverse`
- automatic annotation (`scRNA_SingleR.r`): `SingleR`, `SingleCellExperiment`, `scater`, `celldex` (built-in human/mouse reference), `Seurat`, `tidyverse`
- pathway enrichment (`scRNA_spRNA_enrichment.r` + companion `enrichment_analysis.r`): the main script only prepares the inputs; the enrichment itself is done by the companion script, which needs `clusterProfiler`, `enrichplot`, `GseaVis` and `org.Hs.eg.db` / `org.Mm.eg.db`
- pathway activity scoring (`scRNA_GSVA_enrichment.r`): `GSVA`, `msigdbr` (GO:BP / KEGG MEDICUS / REACTOME gene sets), `limma` (group comparison), `ComplexHeatmap` (heatmaps), `Seurat`, `tidyverse`
- custom gene set scoring (`scRNA_genesetScore.r`): `GSVA` (GSVA + ssGSEA), `AUCell`, `Seurat` (AddModuleScore), `tidyverse`

**Python side** (`scRNA_spRNA_h5ad_preprocess.py`, `scRNA_composition_sccoda.py`)

- `scCODA`, `tensorflow` / `tensorflow-probability`
- `anndata` (mind the zarr / numpy compatibility: `anndata < 0.12` rejects zarr ≥ 3, `anndata ≥ 0.13` requires numpy ≥ 2.1)
- h5ad preprocessing: `scanpy`, `scipy`, `numpy`, `pandas`
- composition analysis: `pandas`, `matplotlib`, `seaborn`
- the sankey plot additionally needs `Rscript` with `ggalluvial`; if missing, that plot is skipped and everything else is unaffected

**Interpreters** — shebangs are `#!/usr/bin/env Rscript` and `#!/usr/bin/env python3`, so the environment that holds these packages must be on `PATH` (e.g. `export PATH=/path/to/env/bin:$PATH`) or you should call the interpreter explicitly: `Rscript script.r …`, `/path/to/env/bin/python script.py …`.

## Quick start

```bash
# optional: start from .h5ad
./scRNA_spRNA_h5ad_preprocess.py -I h5ad_dir -V gene_symbols -P 00. -O 00_preprocess
./scRNA_h5ad2rds.r               -I 00_preprocess -P 00. -O 00_rds

# optional: Broad Single Cell Portal export (dense matrix + metadata)
./scRNA_mtx2rds.r                -E exprMatrix.tsv.gz -m metadata.tsv -CC NAME -P 00. -O 00_rds

# optional: sub-cluster each cell type of an annotated object (stage 1 -> clustree)
./scRNA_subcelltype_clustering.r -I 05_annotation/05.seu_annotated.rds -AC celltype -P 06. -O 06_subtype

# optional: which genes differ between two groups within each cell type (pseudobulk, DESeq2)
./scRNA_pseudobulk_DE.r -I 05_annotation/05.seu_annotated.rds -CC celltype -SC orig.ident -GC group -G1 treat -G2 control -P 07. -O 07_pseudobulk

# optional: cell-type composition per sample + comparison between two groups
./scRNA_spRNA_cellPro.r -I 05_annotation/05.seu_annotated.rds -SC orig.ident -CC celltype -GC group -GO control,treat -EC Unknown -P 08. -O 08_cellPro

# optional: how well do target genes separate two groups (AUC / ROC curves)
./scRNA_ROC.r -I 05_annotation/05.seu_annotated.rds -SC orig.ident -G group.csv -F genelist.txt -C control,treat -P 09. -O 09_ROC

# optional: cell-cell communication (run twice: first without -NP, then with it)
./scRNA_cellchat_single.r -I 05_annotation/05.seu_annotated.rds -AC celltype -S mouse -TP triMean -TM 0.1 -NP 3 -T 20 -P 10. -O 10_cellChat

# optional: automatic cell-type labels from a reference dataset (SingleR)
./scRNA_SingleR.r -I 04_clustering/04.seu_res0.6_clustered.rds -S mouse -P 05. -O 05_SingleR

# optional: GO/KEGG enrichment + GSEA on the DE results (per cell type)
./scRNA_spRNA_enrichment.r -I 05_annotation/05.seu_annotated.rds -DER 06_pseudobulk/06.KO_vs_WT_DE_raw.csv -CC celltype -GC group -G control,treat -A RNA -S mouse -IC Unknown -P 07. -O 07_enrichment

# optional: pathway activity comparison between two groups (GSVA, no DE result needed)
./scRNA_GSVA_enrichment.r -I 05_annotation/05.seu_annotated.rds -AC celltype -SC orig.ident -GC group -G control,treat -S mouse -P 08. -O 08_GSVA

# optional: score your own gene set (GSVA / ssGSEA / AddModuleScore / AUCell)
./scRNA_genesetScore.r -I 05_annotation/05.seu_annotated.rds -G geneSet.txt -AC celltype -A RNA -P 09. -O 09_genesetScore

# main pipeline
./scRNA_QC_step1.r          -I input_list.csv -DT 10X -P 01. -O 01_QC
./scRNA_dedoublets_step2.r  -I 01_QC/01.seu_afterQC.rds  -P 02. -O 02_deContamination
./scRNA_Nor_harmony_step3.r -I 02_deContamination/02.seu_clean.rds -P 03. -O 03_Nor_Harmony

# step4 has two stages: exploration first, then the final resolution
./scRNA_clustering_step4.r  -I 03_Nor_Harmony/03.seu_harmony.rds -P 04. -O 04_clustering
./scRNA_clustering_step4.r  -I 03_Nor_Harmony/03.seu_harmony.rds -Res 0.7 -P 04. -O 04_clustering

# step5: explore with a marker list, then write the final annotation
./scRNA_annotate_step5.r    -I 04_clustering/04.seu_res0.7_clustered.rds -MB markers.tsv    -P 05. -O 05_annotation
./scRNA_annotate_step5.r    -I 04_clustering/04.seu_res0.7_clustered.rds -MA annotation.tsv -P 05. -O 05_annotation

# standalone scripts
./scRNA_gene_expr.r         -I 05_annotation/05.seu_annotated.rds -G GENE1,GENE2,GENE3 -P 06. -O 06_gene_expr
./scRNA_gene_DE.r           -I 05_annotation/05.seu_annotated.rds -G GENE1 -CT celltypeA,celltypeB -P 07. -O 07_gene_DE
./scRNA_composition_sccoda.py -I input.h5ad -S sample -CT celltype -GC group --deterministic -P 08. -O 08_composition

# single-file format round trip
./scRNA_h5ad_seurat_convert.r -I 05_annotation/05.seu_annotated.rds -O 09_export -P export.
```

Every script is self-documenting: run `<script> -h` for the full option list.

### Input list for step1

```csv
path,sample,group
/path/to/sample1/filtered_feature_bc_matrix,sample1,control
/path/to/sample2/filtered_feature_bc_matrix,sample2,treat
```

`-DT 10X` expects a 10x directory per sample; `-DT h5` expects one `.h5` file per sample. `sample` must be unique and `group` is stored in the object metadata.

## Conventions

- Scripts whose filename carries a `stepN` prefix form the main pipeline and are meant to be chained in order; the remaining scripts are standalone and can be applied to any suitable input.
- Stages communicate through files — intermediate objects can be reused, and each stage's output prefix (e.g. `01.`, `02.`) matches its step number.
- Existing key outputs are treated as a cache (re-plot only); delete the corresponding `rds` to force a full re-run.
- For Harmony, the batch variable must be the sample/donor column — never the biological group column (step3 rejects `-BC group` so the group effect is not regressed out).
- `step5` does not recompute DE; it reuses the `clusters_DE/` tables produced by `step4`.
- When starting from a **Broad Single Cell Portal** export, `scRNA_mtx2rds.r` assembles the Seurat object directly (cell-name column defaults to `NAME`; a leading `TYPE` row is skipped).
- `scRNA_subcelltype_clustering.r` refines an annotated object: stage 1 writes one clustree per cell type, stage 2 applies the resolution you recorded per cell type (CSV: `celltype`, `resolution`) and runs sub-celltype DE.
- `scRNA_pseudobulk_DE.r` answers "which genes differ between two groups within each cell type". It aggregates raw counts per (cell type × sample) and runs DESeq2, so `-SC` must be the **biological replicate**. Cell-level tests treat every cell as a replicate and give inflated p-values; the pseudobulk table is the one whose conclusions hold for new samples.
- Cell-cell communication is a two-step affair: `scRNA_cellchat_single.r` models one object (run it **first without `-NP`** to get the `selectK` plot, read the number of patterns off it, then run again with `-NP`) and produces a `cellchat.rds`; `scRNA_cellchat_compare.r` takes two such objects and compares them. Both objects must have exactly the same cell identities.
- `scRNA_ROC.r` is a discrimination analysis: the unit is the sample, and AUC says how well a gene separates the two groups — it is not a significance test. Use `scRNA_gene_DE.r` (one gene, sample-level Wilcoxon) or `scRNA_pseudobulk_DE.r` (per cell type, DESeq2) when you need p values.
- `scRNA_spRNA_cellPro.r` answers "how much of each cell type is in each sample, and do the proportions differ between two groups". The unit is the sample (`-SC`), not the cell. It is a fast description plus a simple Wilcoxon test; for a formal composition test with credible intervals use `scRNA_composition_sccoda.py` instead.
- `scRNA_SingleR.r` gives cell-type names automatically by matching each cluster against a reference dataset. The input must already carry `seurat_clusters`, so it expects the output of step 4. Treat its labels as a first guess to be cross-checked against markers — it complements, rather than replaces, the marker-based annotation of step 5.
- Doublet removal comes in two flavours — `scRNA_dedoublets_step2.r` (scDblFinder, per-sample, estimates the doublet rate, outputs counts only) and `scRNA_DoubleFinder.r` (DoubletFinder, whole object, doublet rate given by hand). Pick one, never both on the same data; the second is for cross-checking with an independent algorithm.
- `scRNA_spRNA_enrichment.r` turns a DE table into pathway statements: GO/KEGG enrichment on the significant genes, GSEA on all genes ranked by `-GM`. Its DE input must carry `gene`, `cluster`, `p_val_adj` and `avg_log2FC`, so the output of `scRNA_pseudobulk_DE.r` feeds straight in. The enrichment itself runs in a separate background process — the script finishing does not mean the results are ready.
- `scRNA_GSVA_enrichment.r` answers the same question without a DE table: it scores pathway activity per cell type and compares the two groups. The order of `-G` matters — the first group is the control (it ends up on the left of the heatmap), the second is the treatment. Pathways shown in a heatmap are not necessarily significant: the script falls back from `adj.P.Val < 0.05` to `P.Value < 0.05` to plain `|logFC|` ranking, so check `adj.P.Val` in the DE table before quoting a pathway.
- `scRNA_genesetScore.r` scores a gene set you provide, not a pathway library. It reports four scores for the same set (GSVA, ssGSEA, AddModuleScore, AUCell) at two levels (cell type, single cell), and prints how many of your genes were found in the object — check that count first, since a gene-set / object mismatch makes the scores meaningless.
- When starting from `.h5ad`, run `scRNA_spRNA_h5ad_preprocess.py` first: `scRNA_h5ad2rds.r` needs a `layers["counts"]` layer, and a Seurat object whose only layer is named `X` cannot be used by the pipeline (step1 reads `layer = "counts"`).

## Visualization examples

The `Images/` directory contains example outputs for the plotting scripts, one subdirectory per script that produces figures — see the "Example outputs" column above. All figures are produced by the scripts themselves and are meant to show the output format and layout.

**PNG only.** Each script writes the same figure twice (`.png` and `.pdf`), but the example set keeps **only the PNG** and the matching PDFs have been removed, to keep the repository light and browser-friendly. So although the per-script output tables below list `.png/.pdf`, you will find just the `.png` in `Images/`; re-run the script if you need the PDF.

## Documentation

- [`USAGE.md`](USAGE.md) — complete manual (Chinese): per-script parameters, output tables, caching behaviour and known pitfalls.
- Each script's `-h` output — authoritative and always in sync with the code.

## Built with

This library was developed with the assistance of the **Hermes agent** ([Nous Research](https://hermes-agent.nousresearch.com/docs)) together with large language models including **DeepSeek** and **MiniMax**. Every script was reviewed, actually run and validated by the author before being included here.

## Roadmap

The library is under active development — **more scripts will be added over time** (additional analysis stages and downstream utilities). Suggestions and requests are welcome.

---

Copyright © 2026 OYJH. All Rights Reserved.
