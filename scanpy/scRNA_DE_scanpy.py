#!/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python
"""
Single-cell RNA-seq differential expression analysis script based on scanpy.

Example usage:
  scRNA_DE_scanpy.py -I preprocessed.h5ad -CK leiden -O DE_output -P lung_
"""
import sys, argparse, warnings, os, logging, time
from datetime import datetime
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)


def main():
    parser = argparse.ArgumentParser(
        description="""
        A script for differential expression analysis of single-cell data using scanpy.
        Example: 
          scRNA_DE_scanpy.py -I preprocessed.h5ad -CK leiden -O DE_output -P lung_
        """
    )

    # Required arguments
    parser.add_argument("-I", "--input", required=True, type=str,
                        help="Path to an h5ad file for AnnData object (already clustered). (Required)")

    # Optional arguments
    parser.add_argument("-CK", "--cluster_key", type=str, default="leiden",
                        help="Key for cluster labels in adata.obs (default: leiden)")
    parser.add_argument("-M", "--method", type=str, default="wilcoxon",
                        choices=["wilcoxon", "t-test", "logreg"],
                        help="Method for DE analysis: 'wilcoxon', 't-test', 'logreg' (default: wilcoxon)")
    parser.add_argument("-n_genes", "--n_genes", type=int, default=20,
                        help="Number of top genes to show in marker plot (default: 20)")
    parser.add_argument("-top_n", "--top_n", type=int, default=10,
                        help="Number of top markers per cluster for dotplot (default: 10)")
    parser.add_argument("-P", "--prefix", type=str, default="",
                        help="Prefix for output files (default: none)")
    parser.add_argument("-O", "--outpath", type=str, default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("-dpi", "--dpi", type=int, default=300,
                        help="DPI for saved figures (default: 300)")
    parser.add_argument("-F", "--force", action="store_true", default=False,
                        help="Force re-run DE analysis even if existing results found (default: False)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file = args.input
    cluster_key = args.cluster_key
    method = args.method
    n_genes = args.n_genes
    top_n = args.top_n
    prefix = args.prefix
    outpath = args.outpath
    dpi = args.dpi
    force = args.force

    # ------------------------------------------------------------------ #
    # Setup environment                                                    #
    # ------------------------------------------------------------------ #
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams["font.sans-serif"] = "Arial"
    sc.set_figure_params(dpi=80, dpi_save=dpi, figsize=(8, 8), facecolor="white")
    sc.settings.verbosity = 1

    # Create output directories
    os.makedirs(outpath, exist_ok=True)
    de_dir = os.path.join(outpath, "{}clusters_DE".format(prefix))
    os.makedirs(de_dir, exist_ok=True)

    # ------------------------------------------------------------------ #
    # Print run information                                                #
    # ------------------------------------------------------------------ #
    print("\n======================== scRNA DE Analysis ========================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("DE output dir   :", de_dir)
    print("Output prefix   :", prefix)
    print("Cluster key     :", cluster_key)
    print("DE method       :", method)
    print("\nCopyright (c) 2026 OYJH. All Rights Reserved.")
    print("\n====================================================================\n")

    de_start = time.time()

    # ------------------------------------------------------------------ #
    # Load data                                                            #
    # ------------------------------------------------------------------ #
    print(">>> Loading data\n")
    adata = sc.read_h5ad(input_file)
    print("Data shape: {} cells x {} genes\n".format(adata.n_obs, adata.n_vars))

    # Check if cluster_key exists
    if cluster_key not in adata.obs.columns:
        print("Error: cluster_key '{}' not found in adata.obs.columns".format(cluster_key),
              file=sys.stderr)
        print("Available columns:", list(adata.obs.columns), file=sys.stderr)
        sys.exit(1)

    # Get cluster information
    clusters = adata.obs[cluster_key].astype('category').cat.categories
    n_clusters = len(clusters)
    print("Found {} clusters: {}\n".format(n_clusters, list(clusters)))

    # ------------------------------------------------------------------ #
    # Differential expression analysis                                      #
    # ------------------------------------------------------------------ #
    print(">>> Differential expression analysis\n")

    # Check if DE analysis already exists
    has_existing_de = 'rank_genes_groups' in adata.uns
    need_run_de = True  # Default: need to run DE analysis

    if has_existing_de and not force:
        # Check if existing results match the cluster_key
        existing_params = adata.uns['rank_genes_groups'].get('params', {})
        existing_groupby = existing_params.get('groupby', '')
        
        if existing_groupby == cluster_key:
            print("Found existing DE analysis results for groupby='{}'.".format(cluster_key))
            print("Using existing results (skip re-analysis).\n")
            need_run_de = False
        else:
            print("Existing DE results are for groupby='{}', not '{}'.".format(
                  existing_groupby, cluster_key))
            print("Auto re-running analysis with new cluster_key.\n")
            need_run_de = True

    if force and has_existing_de:
        print("Force re-run enabled. Re-running analysis.\n")
        need_run_de = True

    if need_run_de:
        print("Running DE analysis with method: {} (CPU)\n".format(method))
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method)

    # ------------------------------------------------------------------ #
    # Plot 1: Marker genes overview                                         #
    # ------------------------------------------------------------------ #
    print(">>> Plot 1: Marker genes overview\n")

    sc.pl.rank_genes_groups(adata, n_genes=n_genes, sharey=False, ncols=4, show=False)
    plt.savefig("{}/{}marker_genes.png".format(outpath, prefix), dpi=dpi, bbox_inches='tight')
    plt.savefig("{}/{}marker_genes.pdf".format(outpath, prefix), bbox_inches='tight')
    plt.close()
    print("Plot 1: {}marker_genes.png(pdf) generated successfully!\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Save all markers to CSV                                               #
    # ------------------------------------------------------------------ #
    print(">>> Saving all markers to CSV\n")

    # Get all markers as DataFrame (scanpy native format)
    all_markers = pd.concat(
        [sc.get.rank_genes_groups_df(adata, group=cluster).assign(cluster=cluster)
         for cluster in clusters],
        ignore_index=True
    )

    # Save all markers
    all_markers_file = "{}/{}all_markers.csv".format(de_dir, prefix)
    all_markers.to_csv(all_markers_file, index=False)
    print("File 1: {}all_markers.csv saved successfully!\n".format(prefix))

    # ------------------------------------------------------------------ #
    # Plot 2: Dotplot for each cluster                                      #
    # ------------------------------------------------------------------ #
    print(">>> Plot 2: Generating dotplots for each cluster\n")

    for cluster in clusters:
        print("Generating dotplot for cluster: {} ...".format(cluster))
        
        # Get top markers for this cluster
        cluster_markers_df = sc.get.rank_genes_groups_df(adata, group=cluster)
        cluster_markers = cluster_markers_df['names'].head(top_n).tolist()
        
        if len(cluster_markers) == 0:
            print("  No markers found for cluster {}, skipping.\n".format(cluster))
            continue
        
        print("  Number of markers for cluster {}: {}".format(cluster, len(cluster_markers)))
        
        # Plot dotplot
        sc.pl.dotplot(
            adata,
            var_names=cluster_markers,
            groupby=cluster_key,
            standard_scale="var",
            cmap="Reds",
            show=False,
            dendrogram=True,
            title="Cluster {}".format(cluster)
        )
        
        cluster_clean = str(cluster).replace(" ", "_")
        plt.savefig("{}/{}cluster_{}_dotplot.png".format(de_dir, prefix, cluster_clean),
                   dpi=dpi, bbox_inches='tight')
        plt.savefig("{}/{}cluster_{}_dotplot.pdf".format(de_dir, prefix, cluster_clean),
                   bbox_inches='tight')
        plt.close()
        print("  Dotplot for cluster {} saved.\n".format(cluster))

    print("DE analysis completed!\n")

    # ------------------------------------------------------------------ #
    # Save results                                                          #
    # ------------------------------------------------------------------ #
    print(">>> Saving results\n")

    if need_run_de:
        output_file = "{}/{}DE_analyzed.h5ad".format(outpath, prefix)
        adata.write_h5ad(output_file, compression="gzip")
        print("Data saved to: {}".format(output_file))
    else:
        print("No new analysis performed, h5ad file not modified. Skip saving.\n")

    elapsed = time.time() - de_start
    print("Total time   : {:.1f} s\n".format(elapsed))

    print("============================ All Done ============================\n")


if __name__ == "__main__":
    main()