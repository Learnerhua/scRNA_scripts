# scRNA-seq Analysis Scripts

A general-purpose, **continuously growing** toolkit for single-cell RNA-seq (scRNA-seq) analysis. It combines a step-numbered **Seurat v5 main pipeline** — per-sample QC and merging, doublet removal, normalization and batch integration, clustering and differential expression, cell-type annotation — with **standalone downstream scripts** (target-gene expression visualization, per-sample × per-cell-type statistics, compositional differential abundance, …), plus the companion R helper that renders the composition sankey. Every script takes its inputs as command-line paths (10x directories, `.h5`, `rds`, `h5ad`, tables, …) and writes its results to files (data objects, tables, figures, …), so any stage can be run on its own or chained into a full pipeline. **New scripts are added over time.**

[中文说明 (Chinese)](README_zh.md) · Full manual (Chinese): [`USAGE.md`](USAGE.md)

## Scripts

| Script | Type | Main function | Input → Output | Example outputs |
|---|---|---|---|---|
| `scRNA_QC_step1.r` | pipeline `step1` | Per-sample QC with sample-specific thresholds (quantiles of nFeature_RNA / nCount_RNA + hard cap on percent.mt, optional percent.HB for red-blood-cell removal), then merge QC-passing samples; gene-level `min.cells` filtering is applied **after** merging | 10x directory / `.h5` list (CSV) → merged Seurat object (counts only) + QC violin/scatter plots | `Images/01_QC/` |
| `scRNA_dedoublets_step2.r` | pipeline `step2` | Per-sample doublet detection with **scDblFinder**; optional **DecontX** ambient-RNA decontamination (also per sample) | merged rds → singlet-only rds (`seu_clean.rds`) | `Images/02_deContamination/` |
| `scRNA_Nor_harmony_step3.r` | pipeline `step3` | NormalizeData (LogNormalize) → FindVariableFeatures (vst) → ScaleData → RunPCA → **Harmony** batch correction | clean rds → normalized + integrated rds (`seu_harmony.rds`) | `Images/03_Nor_Harmony/` |
| `scRNA_clustering_step4.r` | pipeline `step4` | Two stages: multi-resolution clustree exploration (omit `-Res`) → clustering + UMAP + `FindAllMarkers` DE (pass `-Res`); Louvain or Leiden | integrated rds → clustered rds + `clusters_DE/` marker tables and plots | `Images/04_clustering/`, `Images/04_clustering_louvain/` |
| `scRNA_annotate_step5.r` | pipeline `step5` | Cell-type annotation, two exclusive modes: `-MB` = marker-list exploration dotplots (writes no rds); `-MA` = cluster → cell-type mapping (writes `seu_annotated.rds`, plus celltype UMAP / dotplot) | clustered rds + marker TSV → annotated rds + plots | `Images/05_annotation/` |
| `scRNA_gene_expr.r` | standalone | Average expression of target genes in every cell type (CSV) plus 5 plot types: Ridge / Violin / Feature / Dot / Heatmap | annotated rds → expression CSV + plots | `Images/gene_expr/` |
| `scRNA_gene_DE.r` | standalone | For one target gene: sample × cell-type mean-expression table, per-celltype z-scored heatmap, one barplot per cell type, and sample-level boxplots with Wilcoxon test + significance bars (each point is a sample, not a cell) | annotated rds → CSV + heatmap / bar / box plots | `Images/gene_DE/` |
| `scRNA_composition_sccoda.py` | standalone | Compositional differential-abundance analysis with **scCODA** (Dirichlet-multinomial Bayesian model); `--deterministic` gives bit-reproducible MCMC | `.h5ad` / sample × cell-type count table / per-cell obs table → effect tables + sankey / boxplot / dispersion / forest plots | `Images/composition_sccoda/` |
| `sankey_composition.r` | helper (**required** by the script above) | Renders the composition sankey (alluvial) plot with `ggalluvial`. Called automatically by `scRNA_composition_sccoda.py`; **must live in the same directory** as it | composition CSV → sankey PNG/PDF | ↑ |

## Pipeline overview

```text
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
- **Statistics at sample level, not cell level** — boxplot comparisons in `scRNA_gene_DE.r` use one point per sample, avoiding pseudoreplication.
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

**Python side** (`scRNA_composition_sccoda.py`)

- `scCODA`, `tensorflow` / `tensorflow-probability`
- `anndata` (mind the zarr / numpy compatibility: `anndata < 0.12` rejects zarr ≥ 3, `anndata ≥ 0.13` requires numpy ≥ 2.1), `numpy`, `pandas`, `matplotlib`, `seaborn`
- the sankey plot additionally needs `Rscript` with `ggalluvial`; if missing, that plot is skipped and everything else is unaffected

**Interpreters** — shebangs are `#!/usr/bin/env Rscript` and `#!/usr/bin/env python3`, so the environment that holds these packages must be on `PATH` (e.g. `export PATH=/path/to/env/bin:$PATH`) or you should call the interpreter explicitly: `Rscript script.r …`, `/path/to/env/bin/python script.py …`.

## Quick start

```bash
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

## Visualization examples

The `Images/` directory contains example outputs for every script (PNG + PDF pairs), one subdirectory per script — see the "Example outputs" column above. All figures are produced by the scripts themselves and are meant to show the output format and layout.

## Documentation

- [`USAGE.md`](USAGE.md) — complete manual (Chinese): per-script parameters, output tables, caching behaviour and known pitfalls.
- Each script's `-h` output — authoritative and always in sync with the code.

## Built with

This library was developed with the assistance of the **Hermes agent** ([Nous Research](https://hermes-agent.nousresearch.com/docs)) together with large language models including **DeepSeek** and **MiniMax**. Every script was reviewed, actually run and validated by the author before being included here.

## Roadmap

The library is under active development — **more scripts will be added over time** (additional analysis stages and downstream utilities). Suggestions and requests are welcome.

---

Copyright © 2026 OYJH. All Rights Reserved.
