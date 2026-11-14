#!/usr/bin/env python
import sys, argparse, warnings, os, logging
from datetime import datetime
import time
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)


def check_gpu():
    """Check GPU availability and print device info."""
    try:
        import cupy as cp
    except ImportError:
        print("Error: cupy is not installed. This pipeline requires a CUDA GPU environment.", file=sys.stderr)
        sys.exit(1)

    if not cp.is_available():
        print("Error: No GPU available. This pipeline requires a CUDA-capable GPU.", file=sys.stderr)
        sys.exit(1)

    props = cp.cuda.runtime.getDeviceProperties(0)
    gpu_name = props['name'].decode() if isinstance(props['name'], bytes) else props['name']
    gpu_mem = round(props['totalGlobalMem'] / 1024**3, 1)
    print("GPU available  : True")
    print("GPU name       :", gpu_name)
    print("GPU memory     : {} GB".format(gpu_mem))
    return cp


def main():
    parser = argparse.ArgumentParser(
        description="GPU-accelerated single-cell RNA-seq analysis pipeline (rapids-singlecell). "
                    "QC -> Doublet detection -> Normalization -> Feature selection -> "
                    "PCA -> Batch correction -> Clustering -> Marker gene identification"
    )

    # Input/Output arguments
    parser.add_argument("-I", "--input", required=True, metavar="H5AD_FILE",
                        help="Path to input h5ad file")
    parser.add_argument("-O", "--outpath", type=str, default=".", metavar="OUTPUT_DIR",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="", metavar="PREFIX",
                        help="Prefix for output files (default: none)")

    # QC arguments (no default values - user must specify after viewing QC plots)
    parser.add_argument("-minG", "--min_genes", type=int, default=None, metavar="INT",
                        help="Minimum number of genes per cell")
    parser.add_argument("-maxG", "--max_genes", type=int, default=None, metavar="INT",
                        help="Maximum number of genes per cell")
    parser.add_argument("-minU", "--min_counts", type=int, default=None, metavar="INT",
                        help="Minimum UMI counts per cell")
    parser.add_argument("-maxU", "--max_counts", type=int, default=None, metavar="INT",
                        help="Maximum UMI counts per cell")
    parser.add_argument("-minC", "--min_cells", type=int, default=None, metavar="INT",
                        help="Minimum cells per gene")
    parser.add_argument("-maxMT", "--max_mt", type=float, default=None, metavar="FLOAT",
                        help="Maximum mitochondrial percentage")
    parser.add_argument("-maxHB", "--max_hb", type=float, default=None, metavar="FLOAT",
                        help="Maximum hemoglobin percentage")

    # Batch correction arguments
    parser.add_argument("-BK", "--batch_key", type=str, default=None, metavar="STR",
                        help="Batch key in adata.obs for batch correction (default: None)")

    # Clustering arguments
    parser.add_argument("-nPC", "--n_pcs", type=int, default=30, metavar="INT",
                        help="Number of PCs for clustering (default: 30)")
    parser.add_argument("-nHVG", "--n_top_genes", type=int, default=2000, metavar="INT",
                        help="Number of highly variable genes (default: 2000)")
    parser.add_argument("-res", "--resolution", type=float, default=1.0, metavar="FLOAT",
                        help="Leiden clustering resolution (default: 1.0)")

    # Output arguments
    parser.add_argument("-raw", "--save_raw", action="store_true",
                        help="Save all genes in adata.X (default: only HVGs)")
    parser.add_argument("-dpi", "--dpi", type=int, default=300, metavar="INT",
                        help="DPI for saved figures (default: 300)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file   = args.input
    outpath      = args.outpath
    prefix       = args.prefix
    min_genes    = args.min_genes
    max_genes    = args.max_genes
    min_counts   = args.min_counts
    max_counts   = args.max_counts
    min_cells    = args.min_cells
    max_mt       = args.max_mt
    max_hb       = args.max_hb
    batch_key    = args.batch_key
    n_pcs        = args.n_pcs
    n_top_genes  = args.n_top_genes
    resolution   = args.resolution
    save_raw     = args.save_raw
    dpi          = args.dpi

    # ------------------------------------------------------------------ #
    # Setup environment                                                    #
    # ------------------------------------------------------------------ #
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams["font.sans-serif"] = "Arial"
    sc.set_figure_params(dpi=80, dpi_save=dpi, figsize=(8, 8), facecolor="white")
    sc.settings.verbosity = 1
    os.makedirs(outpath, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n============= scRNA-seq Analysis Pipeline (GPU) ==============\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix)
    print("Batch key       :", batch_key)
    print("Clustering      : n_pcs={}, n_top_genes={}, resolution={}".format(
          n_pcs, n_top_genes, resolution))
    print("Save all genes  :", save_raw)

    # ------------------------------------------------------------------ #
    # Check GPU environment                                                #
    # ------------------------------------------------------------------ #
    cp = check_gpu()

    # Initialize RMM allocator for better GPU memory management
    try:
        import rmm
        from rmm.allocators.cupy import rmm_cupy_allocator
        rmm.reinitialize(managed_memory=False, pool_allocator=True, devices=0)
        cp.cuda.set_allocator(rmm_cupy_allocator)
        print("RMM allocator   : enabled (pool_allocator=True)")
    except ImportError:
        print("RMM allocator   : not available, using default")

    # Import rapids_singlecell
    try:
        import rapids_singlecell as rsc
    except ImportError:
        print("Error: rapids_singlecell is not installed.", file=sys.stderr)
        sys.exit(1)

    print("\nCopyright (c) 2026 OYJH. All Rights Reserved.")
    print("\n==============================================================\n")

    pipeline_start = time.time()

    # ------------------------------------------------------------------ #
    # Step 1: Load data                                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 1: Loading data\n")
    adata = ad.read_h5ad(input_file)
    print("Data shape: {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # Fail fast: an unknown batch_key would otherwise only surface at the
    # doublet-detection step, after the QC filtering work has been done.
    if batch_key is not None and batch_key not in adata.obs.columns:
        print("Error: batch_key '{}' not found in adata.obs.columns".format(batch_key), file=sys.stderr)
        print("Available columns:", list(adata.obs.columns), file=sys.stderr)
        sys.exit(1)

    # ------------------------------------------------------------------ #
    # Step 2: Quality control                                              #
    # ------------------------------------------------------------------ #
    print(">>> Step 2: Quality control\n")

    adata.var['mt'] = adata.var_names.str.upper().str.startswith('MT-')
    adata.var["hb"] = adata.var_names.str.upper().str.match(r"^HB(?!P)")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

    # Convert to float32 if needed (GPU only supports float32/float64)
    if hasattr(adata.X, 'dtype'):
        if adata.X.dtype not in [np.float32, np.float64, np.bool_]:
            print("Converting data type: {} -> float32".format(adata.X.dtype))
            if hasattr(adata.X, 'toarray'):
                adata.X = adata.X.astype(np.float32)
            else:
                adata.X = adata.X.astype(np.float32)

    # Move data to GPU
    rsc.get.anndata_to_GPU(adata)
    print("Data moved to GPU")

    rsc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'hb', 'ribo'])

    # Move back to CPU for plotting
    rsc.get.anndata_to_CPU(adata)

    # Plot QC metrics before filtering
    features = ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_hb']

    if batch_key is not None and batch_key in adata.obs.columns:
        # Plot with batch grouping
        print("Plotting QC metrics grouped by: {}\n".format(batch_key))

        # Violin plots grouped by batch (2x2 grid, no stripplot)
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
        for ax, feature in zip(axes.flatten(), features):
            sc.pl.violin(adata, keys=feature, groupby=batch_key, stripplot=False,
                         rotation=45, ax=ax, show=False)
        plt.tight_layout()
        plt.savefig("{}/{}before_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        # Scatter plots colored by batch (shared legend at bottom)
        n_batches = len(adata.obs[batch_key].astype('category').cat.categories)
        # Use enough colors for all batches
        if n_batches <= 20:
            colors = plt.cm.tab20(range(n_batches))
        elif n_batches <= 40:
            colors = np.vstack([plt.cm.tab20(range(20)), plt.cm.tab20b(range(n_batches - 20))])
        else:
            colors = plt.cm.hsv(np.linspace(0, 1, n_batches, endpoint=False))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 7))
        batches = adata.obs[batch_key].astype('category')
        batch_categories = batches.cat.categories

        # Plot each scatter
        for i, (x, y, title) in enumerate([
            ('total_counts', 'n_genes_by_counts', 'n_genes vs total_counts'),
            ('total_counts', 'pct_counts_mt', 'pct_mt vs total_counts'),
            ('total_counts', 'pct_counts_hb', 'pct_hb vs total_counts')
        ]):
            for j, batch in enumerate(batch_categories):
                mask = batches == batch
                ax[i].scatter(adata.obs.loc[mask, x], adata.obs.loc[mask, y],
                             c=[colors[j]], s=5, alpha=0.6, label=batch if i == 0 else "")
            ax[i].set_xlabel(x)
            ax[i].set_ylabel(y)
            ax[i].set_title(title)

        # Add shared legend at bottom
        ncol = min(n_batches, 8)
        fig.legend(loc='lower center', ncol=ncol,
                   bbox_to_anchor=(0.5, -0.05), frameon=False, fontsize=8)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.12 + 0.02 * (n_batches // 8))
        plt.savefig("{}/{}before_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()
    else:
        # Plot without grouping
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
        for i, feature in enumerate(features):
            sc.pl.violin(adata, keys=feature, stripplot=False, ax=axes.flatten()[i], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}before_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax[0], show=False)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax[1], show=False)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb', ax=ax[2], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}before_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}before_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

    # Check if all QC parameters are provided
    qc_params = {
        'min_genes': min_genes,
        'max_genes': max_genes,
        'min_counts': min_counts,
        'max_counts': max_counts,
        'min_cells': min_cells,
        'max_mt': max_mt,
        'max_hb': max_hb
    }
    missing_params = [k for k, v in qc_params.items() if v is None]

    if missing_params:
        print("QC plots saved. Please specify the following parameters:")
        for p in missing_params:
            print("  --{} / -{}".format(p, {
                'min_genes': 'minG', 'max_genes': 'maxG',
                'min_counts': 'minU', 'max_counts': 'maxU',
                'min_cells': 'minC', 'max_mt': 'maxMT', 'max_hb': 'maxHB'
            }[p]))
        print("\nExample: python {} -I {} -O {} -minG 200 -maxG 2500 -minU 500 -maxU 10000 -minC 3 -maxMT 5 -maxHB 5\n".format(
              sys.argv[0], input_file, outpath))
        sys.exit(0)

    print("QC thresholds: genes [{}, {}], counts [{}, {}], mt<={}, hb={}\n".format(
          min_genes, max_genes, min_counts, max_counts, max_mt, max_hb))

    # Apply filtering
    print("Before filtering: {} cells x {} genes".format(adata.n_obs, adata.n_vars))
    adata = adata[
        (adata.obs['n_genes_by_counts'] >= min_genes) &
        (adata.obs['n_genes_by_counts'] <= max_genes) &
        (adata.obs['total_counts'] >= min_counts) &
        (adata.obs['total_counts'] <= max_counts) &
        (adata.obs['pct_counts_mt'] <= max_mt) &
        (adata.obs['pct_counts_hb'] <= max_hb), :
    ].copy()

    # Move back to GPU and filter genes
    rsc.get.anndata_to_GPU(adata)
    rsc.pp.filter_genes(adata, min_cells=min_cells)
    rsc.get.anndata_to_CPU(adata)
    print("After filtering : {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # Plot QC metrics after filtering
    if batch_key is not None and batch_key in adata.obs.columns:
        # Violin plots grouped by batch (2x2 grid, no stripplot)
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
        for ax, feature in zip(axes.flatten(), features):
            sc.pl.violin(adata, keys=feature, groupby=batch_key, stripplot=False,
                         rotation=45, ax=ax, show=False)
        plt.tight_layout()
        plt.savefig("{}/{}after_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        # Scatter plots colored by batch (shared legend at bottom)
        n_batches = len(adata.obs[batch_key].astype('category').cat.categories)
        # Use enough colors for all batches
        if n_batches <= 20:
            colors = plt.cm.tab20(range(n_batches))
        elif n_batches <= 40:
            colors = np.vstack([plt.cm.tab20(range(20)), plt.cm.tab20b(range(n_batches - 20))])
        else:
            colors = plt.cm.hsv(np.linspace(0, 1, n_batches, endpoint=False))

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 7))
        batches = adata.obs[batch_key].astype('category')
        batch_categories = batches.cat.categories

        # Plot each scatter
        for i, (x, y, title) in enumerate([
            ('total_counts', 'n_genes_by_counts', 'n_genes vs total_counts'),
            ('total_counts', 'pct_counts_mt', 'pct_mt vs total_counts'),
            ('total_counts', 'pct_counts_hb', 'pct_hb vs total_counts')
        ]):
            for j, batch in enumerate(batch_categories):
                mask = batches == batch
                ax[i].scatter(adata.obs.loc[mask, x], adata.obs.loc[mask, y],
                             c=[colors[j]], s=5, alpha=0.6, label=batch if i == 0 else "")
            ax[i].set_xlabel(x)
            ax[i].set_ylabel(y)
            ax[i].set_title(title)

        # Add shared legend at bottom
        ncol = min(n_batches, 8)
        fig.legend(loc='lower center', ncol=ncol,
                   bbox_to_anchor=(0.5, -0.05), frameon=False, fontsize=8)
        plt.tight_layout()
        plt.subplots_adjust(bottom=0.12 + 0.02 * (n_batches // 8))
        plt.savefig("{}/{}after_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()
    else:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 10))
        for i, feature in enumerate(features):
            sc.pl.violin(adata, keys=feature, stripplot=False, ax=axes.flatten()[i], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}after_QC_violin.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_violin.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 5))
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', ax=ax[0], show=False)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', ax=ax[1], show=False)
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_hb', ax=ax[2], show=False)
        plt.tight_layout()
        plt.savefig("{}/{}after_QC_scatter.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}after_QC_scatter.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
        plt.close()

    # ------------------------------------------------------------------ #
    # Step 3: Doublet detection                                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 3: Doublet detection (Scrublet)\n")

    rsc.get.anndata_to_GPU(adata)
    if batch_key is None:
        rsc.pp.scrublet(adata)
    elif batch_key in adata.obs.columns:
        rsc.pp.scrublet(adata, batch_key=batch_key)
    else:
        print("Error: batch_key '{}' not found in adata.obs.columns".format(batch_key), file=sys.stderr)
        sys.exit(1)

    n_doublets = adata.obs["predicted_doublet"].sum()
    print("Detected {} doublets ({:.2f}%)".format(n_doublets, 100*n_doublets/adata.n_obs))
    adata = adata[adata.obs["predicted_doublet"] == False].copy()
    rsc.get.anndata_to_GPU(adata)
    print("After removal  : {} cells\n".format(adata.n_obs))

    # ------------------------------------------------------------------ #
    # Step 4: Normalization                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 4: Normalization\n")

    # Move to CPU before converting to dense array to avoid GPU OOM
    rsc.get.anndata_to_CPU(adata)
    adata.layers['counts'] = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X.copy()
    rsc.get.anndata_to_GPU(adata)

    rsc.pp.normalize_total(adata, target_sum=1e4)
    rsc.pp.log1p(adata)

    # Move to CPU before converting to dense array to avoid GPU OOM
    rsc.get.anndata_to_CPU(adata)
    adata.layers['log1p'] = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X.copy()
    rsc.get.anndata_to_GPU(adata)

    print("Normalized (target_sum=1e4) and log1p transformed\n")

    # ------------------------------------------------------------------ #
    # Step 5: Feature selection                                            #
    # ------------------------------------------------------------------ #
    print(">>> Step 5: Feature selection (HVGs)\n")

    if batch_key is None:
        rsc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='seurat_v3')
    elif batch_key in adata.obs.columns:
        rsc.pp.highly_variable_genes(adata, batch_key=batch_key, n_top_genes=n_top_genes, flavor='seurat_v3')
    else:
        print("Error: batch_key '{}' not found in adata.obs.columns".format(batch_key), file=sys.stderr)
        sys.exit(1)

    n_hvgs = adata.var['highly_variable'].sum()
    print("Found {} HVGs out of {} genes\n".format(n_hvgs, adata.n_vars))

    # Plot HVGs before filtering (move to CPU for plotting)
    rsc.get.anndata_to_CPU(adata)
    sc.pl.highly_variable_genes(adata, show=False)
    plt.savefig("{}/{}HVGs.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}HVGs.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()
    rsc.get.anndata_to_GPU(adata)

    # Save full data to raw
    adata.raw = adata

    # Filter to HVGs only (unless save_raw is True)
    if not save_raw:
        rsc.pp.filter_highly_variable(adata)
        print("Filtered to {} HVGs in adata.X\n".format(adata.n_vars))
    else:
        print("Keeping all genes in adata.X (save_raw=True)\n")

    # ------------------------------------------------------------------ #
    # Step 6: Scaling                                                      #
    # ------------------------------------------------------------------ #
    print(">>> Step 6: Scaling\n")

    # rsc.pp.regress_out(adata, keys=["total_counts", "pct_counts_mt"])
    rsc.pp.scale(adata, max_value=10)
    print("Data scaled (max_value=10)\n")

    # ------------------------------------------------------------------ #
    # Step 7: PCA                                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 7: PCA\n")

    # PCA always uses HVGs (use_highly_variable=True when save_raw=True)
    if save_raw:
        rsc.tl.pca(adata, n_comps=50, use_highly_variable=True)
        print("Computed 50 PCs using {} HVGs (GPU)\n".format(n_hvgs))
    else:
        rsc.tl.pca(adata, n_comps=50)
        print("Computed 50 PCs (GPU)\n")

    # Move to CPU for plotting
    rsc.get.anndata_to_CPU(adata)

    sc.pl.pca_loadings(adata, show=False)
    plt.savefig("{}/{}pca_loadings.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}pca_loadings.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=False, show=False)
    plt.savefig("{}/{}pca_variance_ratio.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}pca_variance_ratio.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    if batch_key is not None:
        sc.pl.pca(adata, color=[batch_key, "pct_counts_mt"], dimensions=[(0, 1), (0, 1)],
                  ncols=2, size=2, wspace=0.4, show=False)
    else:
        sc.pl.pca(adata, color="pct_counts_mt", dimensions=[(0, 1)], size=2, show=False)
    plt.savefig("{}/{}pca.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}pca.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # ------------------------------------------------------------------ #
    # Step 8: Batch correction (Harmony)                                    #
    # ------------------------------------------------------------------ #
    print(">>> Step 8: Batch correction\n")

    if batch_key is not None:
        rsc.get.anndata_to_GPU(adata)
        rsc.pp.harmony_integrate(adata, key=batch_key, dtype=cp.float32)
        print("Harmony batch correction completed (batch_key={})\n".format(batch_key))
    else:
        print("No batch key specified, skipping batch correction\n")

    # ------------------------------------------------------------------ #
    # Step 9: Neighbor graph                                               #
    # ------------------------------------------------------------------ #
    print(">>> Step 9: Building neighbor graph\n")

    rsc.get.anndata_to_GPU(adata)
    if 'X_pca_harmony' in adata.obsm:
        rsc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca_harmony')
        print("Using X_pca_harmony, n_pcs={}\n".format(n_pcs))
    else:
        rsc.pp.neighbors(adata, n_pcs=n_pcs, use_rep='X_pca')
        print("Using X_pca, n_pcs={}\n".format(n_pcs))

    # ------------------------------------------------------------------ #
    # Step 10: Clustering                                                  #
    # ------------------------------------------------------------------ #
    print(">>> Step 10: Clustering (Leiden)\n")

    for res in [0.2, 0.5, 0.8, 1.0, 1.2, 1.5]:
        rsc.tl.leiden(adata, resolution=res, key_added="leiden_{}".format(res))

    rsc.tl.umap(adata)

    # Move to CPU for plotting
    rsc.get.anndata_to_CPU(adata)

    sc.pl.umap(adata, color=["leiden_0.2", "leiden_0.5", "leiden_0.8",
                             "leiden_1.0", "leiden_1.2", "leiden_1.5"],
               ncols=3, wspace=0.4, hspace=0.3, show=False)
    plt.savefig("{}/{}umap_clusters.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}umap_clusters.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    rsc.get.anndata_to_GPU(adata)
    rsc.tl.leiden(adata, resolution=resolution)
    n_clusters = len(adata.obs["leiden"].cat.categories)
    print("Final resolution={}, {} clusters detected\n".format(resolution, n_clusters))

    # ------------------------------------------------------------------ #
    # Step 11: UMAP visualization                                          #
    # ------------------------------------------------------------------ #
    print(">>> Step 11: UMAP visualization\n")

    rsc.tl.umap(adata)
    rsc.get.anndata_to_CPU(adata)

    # Plot leiden clustering with labels on UMAP
    sc.pl.umap(adata, color=["leiden"],
               title="Leiden (resolution: {})".format(resolution), legend_loc='on data', show=False)
    plt.savefig("{}/{}umap_leiden.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}umap_leiden.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    # Plot batch_key separately if provided
    if batch_key is not None:
        sc.pl.umap(adata, color=[batch_key], title=batch_key, show=False)
        plt.savefig("{}/{}umap_{}.png".format(outpath, prefix, batch_key), dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}umap_{}.pdf".format(outpath, prefix, batch_key), dpi=dpi, bbox_inches='tight')
        plt.close()

    # ------------------------------------------------------------------ #
    # Step 12: Marker genes                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 12: Marker gene identification\n")
    print("Method: Wilcoxon rank-sum test (CPU, searching all genes in adata.raw)\n")

    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", use_raw=True)

    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, ncols=4, show=False)
    plt.savefig("{}/{}marker_genes.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}marker_genes.pdf".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.close()

    clusters = adata.obs["leiden"].cat.categories
    all_markers = pd.concat(
        [sc.get.rank_genes_groups_df(adata, group=cluster).assign(cluster=cluster)
         for cluster in clusters],
        ignore_index=True
    )
    all_markers.to_csv("{}/{}all_markers.csv".format(outpath, prefix), index=False)
    print("Marker genes saved to {}all_markers.csv\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Step 13: Save results                                                #
    # ------------------------------------------------------------------ #
    print(">>> Step 13: Saving results\n")

    output_file = "{}/{}preprocessed.h5ad".format(outpath, prefix)
    adata.write_h5ad(output_file, compression="gzip")
    print("Data saved to: {}".format(output_file))
    print("Final shape  : {} cells x {} genes".format(adata.n_obs, adata.n_vars))

    elapsed = time.time() - pipeline_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("==============================================================")
    print("Pipeline completed successfully!")
    print("==============================================================\n")


if __name__ == "__main__":
    main()
