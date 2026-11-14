#!/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python
"""
Single-cell RNA-seq annotation script based on scanpy.
Functionality equivalent to scRNA_annotate_step6.r

Example usage:
  # Before annotation - plot markers by cluster
  scRNA_annotate_scanpy.py -I preprocessed.h5ad -MB Markers.txt -DW 12 -DH 8 -P 06. -O 06_annotation

  # After annotation - annotate clusters and plot
  scRNA_annotate_scanpy.py -I preprocessed.h5ad -MA 06_annotation/annotation.txt -P 06. -O 06_annotation
"""
import sys, argparse, warnings, os, logging
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
        A script for manual annotation of single-cell data using scanpy.
        Example: 
          scRNA_annotate_scanpy.py -I preprocessed.h5ad -MB Markers.txt -DW 12 -DH 8 -P 06. -O 06_annotation
          scRNA_annotate_scanpy.py -I preprocessed.h5ad -MA annotation.txt -P 06. -O 06_annotation
        """
    )

    # Required arguments
    parser.add_argument("-I", "--input", required=True, type=str,
                        help="Path to an h5ad file for AnnData object used for scRNA data annotation. (Required)")

    # Marker file arguments (mutually exclusive)
    parser.add_argument("-MB", "--markers_before", type=str, default=None,
                        help="Path to a tsv file used for annotation before cluster labeling. "
                             "Contains 2 columns: CellType (col1) and marker genes separated by commas (col2).")
    parser.add_argument("-MA", "--markers_after", type=str, default=None,
                        help="Path to a tsv file used for annotation after cluster labeling. "
                             "Contains 3 columns: Cluster, Markers (comma-separated), CellType.")

    # Core arguments
    parser.add_argument("-CK", "--cluster_key", type=str, default="leiden",
                        help="Key for cluster labels in adata.obs (default: leiden)")
    parser.add_argument("-O", "--outpath", type=str, default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("-P", "--prefix", type=str, default="",
                        help="Prefix for output files (default: none)")

    # Plotting arguments
    parser.add_argument("-S", "--split", action="store_true", default=False,
                        help="Whether to split the markers by cell type and generate separate dotplots. "
                             "Use when marker count per cell type is large. (default: False)")
    parser.add_argument("-CO", "--celltype_order", type=str, default=None,
                        help="Order of cell types in the dot plot (comma- or space-separated, "
                             "e.g., 'B_cell,T_cell,NK').")
    parser.add_argument("--dendrogram", dest="dendrogram", action="store_true", default=True,
                        help="Show the dendrogram in the dotplot before annotation (default).")
    parser.add_argument("--no-dendrogram", dest="dendrogram", action="store_false",
                        help="Hide the dendrogram in the dotplot before annotation.")
    parser.add_argument("-VR", "--var_group_rotation", type=int, default=45,
                        help="Rotation angle for cell type labels in dotplot (default: 45)")
    parser.add_argument("-DW", "--dotplot_width", type=int, default=None,
                        help="Width of the dotplot (default: auto)")
    parser.add_argument("-DH", "--dotplot_height", type=int, default=None,
                        help="Height of the dotplot (default: auto)")
    parser.add_argument("-dpi", "--dpi", type=int, default=300,
                        help="DPI for saved figures (default: 300)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    input_file = args.input
    markers_before = args.markers_before
    markers_after = args.markers_after
    split = args.split
    dotplot_width = args.dotplot_width
    dotplot_height = args.dotplot_height
    var_group_rotation = args.var_group_rotation
    show_dendrogram = args.dendrogram
    prefix = args.prefix
    outpath = args.outpath
    dpi = args.dpi
    cluster_key = args.cluster_key

    if args.celltype_order:
        celltype_order = [x.strip() for x in args.celltype_order.replace(',', ' ').split()]
    else:
        celltype_order = None

    # Validate arguments
    if split and markers_before is None and markers_after is None:
        print("Error: If you want to split markers by cell type, you must provide -MB parameter.", 
              file=sys.stderr)
        sys.exit(1)

    if (markers_before is None) == (markers_after is None):
        print("Error: Either -MB or -MA must be specified, but not both.", file=sys.stderr)
        sys.exit(1)

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
    print("\n======================== scRNA Annotation ========================\n")
    print("Current time    :", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Working directory:", os.getcwd())
    print("Input file      :", input_file)
    print("Output directory:", outpath)
    print("Output prefix   :", prefix)
    print("Cluster key     :", cluster_key)
    print("\nCopyright (c) 2026 OYJH. All Rights Reserved.")
    print("\n===================================================================\n")

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

    genes_in_data = adata.var_names.tolist()

    # ------------------------------------------------------------------ #
    # Before annotation - plot markers by cluster                          #
    # ------------------------------------------------------------------ #
    if markers_before is not None and not split:
        print(">>> Plotting markers before annotation (combined dotplot)\n")

        # Read marker file (2 columns: CellType, Markers)
        markers_df = pd.read_csv(markers_before, sep="\t", header=None,
                                 names=["CellType", "Markers"], dtype=str)

        # Check genes
        all_markers = []
        for m in markers_df["Markers"]:
            all_markers.extend([g.strip() for g in m.split(",")])
        all_markers = list(set(all_markers))

        existing_genes = [g for g in all_markers if g in genes_in_data]
        missing_genes = [g for g in all_markers if g not in genes_in_data]

        print("Total input genes:", len(all_markers))
        print("Genes found in AnnData:", len(existing_genes))
        print("Genes missing:", len(missing_genes))
        if missing_genes:
            print("Missing genes:", missing_genes[:20], "..." if len(missing_genes) > 20 else "")

        # Deduplicate markers per cell type
        dedup_df = markers_df.drop_duplicates(subset=["CellType", "Markers"])

        # Sort by cell type order if provided
        if celltype_order:
            dedup_df["CellType"] = pd.Categorical(dedup_df["CellType"], categories=celltype_order)
            dedup_df = dedup_df.sort_values("CellType").dropna(subset=["CellType"])
        else:
            dedup_df = dedup_df.sort_values("CellType")

        # Build marker list with deduplication across groups
        markers_list = {}
        seen_genes = set()
        for _, row in dedup_df.iterrows():
            ct = row["CellType"]
            genes = [g.strip() for g in row["Markers"].split(",")]
            new_genes = [g for g in genes if g not in seen_genes and g in genes_in_data]
            if new_genes:
                markers_list[ct] = new_genes
                seen_genes.update(new_genes)

        # Check empty groups
        empty_groups = [ct for ct, genes in markers_list.items() if len(genes) == 0]
        if empty_groups:
            print("Warning: The following cell types have no marker genes after deduplication:",
                  empty_groups)

        if markers_list:
            # Plot UMAP with cluster labels on data
            print("Plotting UMAP with cluster labels...")
            cluster_key_clean = cluster_key.replace(" ", "_")
            sc.pl.umap(adata, color=[cluster_key], legend_loc='on data',
                       title="Clusters ({})".format(cluster_key), show=False)
            plt.savefig("{}/{}{}_umap.png".format(outpath, prefix, cluster_key_clean),
                       dpi=dpi, bbox_inches='tight')
            plt.savefig("{}/{}{}_umap.pdf".format(outpath, prefix, cluster_key_clean),
                       bbox_inches='tight')
            plt.close()
            print("Plot: {}{}_umap.png(pdf) generated successfully!\n".format(prefix, cluster_key_clean))

            # Plot dotplot grouped by cluster
            # Pre-compute dendrogram to avoid warning (only if show_dendrogram is True)
            if show_dendrogram:
                sc.tl.dendrogram(adata, groupby=cluster_key)
            # Build kwargs for dotplot
            dotplot_kwargs = {
                "var_names": markers_list,
                "groupby": cluster_key,
                "standard_scale": "var",
                "cmap": "Reds",
                "show": False,
                "dendrogram": show_dendrogram,
                "var_group_rotation": var_group_rotation
            }
            if dotplot_width is not None and dotplot_height is not None:
                dotplot_kwargs["figsize"] = (dotplot_width, dotplot_height)
            sc.pl.dotplot(adata, **dotplot_kwargs)
            plt.savefig("{}/{}{}_dotplot.png".format(outpath, prefix, cluster_key_clean),
                       dpi=dpi, bbox_inches='tight')
            plt.savefig("{}/{}{}_dotplot.pdf".format(outpath, prefix, cluster_key_clean),
                       bbox_inches='tight')
            plt.close()
            print("Plot: {}{}_dotplot.png(pdf) generated successfully!\n".format(prefix, cluster_key_clean))

    # ------------------------------------------------------------------ #
    # Before annotation - split markers by cell type                       #
    # ------------------------------------------------------------------ #
    if markers_before is not None and split:
        print(">>> Plotting markers before annotation (split by cell type)\n")

        dotplot_dir = os.path.join(outpath, "{}markers_dotplot".format(prefix))
        os.makedirs(dotplot_dir, exist_ok=True)

        # Read marker file
        markers_df = pd.read_csv(markers_before, sep="\t", header=None,
                                 names=["CellType", "Markers"], dtype=str)

        markers_dict = {}
        for _, row in markers_df.iterrows():
            ct = row["CellType"]
            genes = [g.strip() for g in row["Markers"].split(",")]
            markers_dict[ct] = list(set(genes))

        for ct, markers in markers_dict.items():
            print("Processing cell type:", ct)

            existing_genes = [g for g in markers if g in genes_in_data]
            missing_genes = [g for g in markers if g not in genes_in_data]

            print("  Total input genes:", len(markers))
            print("  Genes found:", len(existing_genes))
            print("  Genes missing:", len(missing_genes))
            if missing_genes:
                print("  Missing genes:", missing_genes[:10], "..." if len(missing_genes) > 10 else "")

            if existing_genes:
                # Build kwargs for dotplot
                dotplot_kwargs = {
                    "var_names": existing_genes,
                    "groupby": cluster_key,
                    "standard_scale": "var",
                    "cmap": "Reds",
                    "show": False,
                    "dendrogram": show_dendrogram,
                    "var_group_rotation": var_group_rotation,
                    "title": "Markers for {}".format(ct)
                }
                if dotplot_width is not None and dotplot_height is not None:
                    dotplot_kwargs["figsize"] = (dotplot_width, dotplot_height)
                sc.pl.dotplot(adata, **dotplot_kwargs)
                ct_clean = ct.replace(" ", "_")
                cluster_key_clean = cluster_key.replace(" ", "_")
                plt.savefig("{}/{}{}_{}_dotplot.png".format(dotplot_dir, prefix, ct_clean, cluster_key_clean),
                           dpi=dpi, bbox_inches='tight')
                plt.savefig("{}/{}{}_{}_dotplot.pdf".format(dotplot_dir, prefix, ct_clean, cluster_key_clean),
                           bbox_inches='tight')
                plt.close()

        print("Plots saved in directory: {}markers_dotplot\n".format(prefix))

    # ------------------------------------------------------------------ #
    # After annotation - annotate clusters and plot                        #
    # ------------------------------------------------------------------ #
    if markers_after is not None:
        print(">>> Annotating clusters and plotting\n")

        # Read annotation file (3 columns: Cluster, Markers, CellType)
        annotation_df = pd.read_csv(markers_after, sep="\t", dtype=str)
        annotation_df.columns = ["Cluster", "Markers", "CellType"]

        # Create cluster to celltype mapping
        cluster_to_celltype = dict(zip(annotation_df["Cluster"], annotation_df["CellType"]))

        # Annotate cells
        adata.obs["celltype"] = adata.obs[cluster_key].astype(str).map(cluster_to_celltype)
        adata.obs["celltype"] = adata.obs["celltype"].fillna("Unknown")

        # Check annotation coverage
        annotated = adata.obs["celltype"] != "Unknown"
        print("Annotated cells: {} / {} ({:.1f}%)".format(
              annotated.sum(), adata.n_obs, 100 * annotated.sum() / adata.n_obs))

        # Plot UMAP colored by celltype
        if 'X_umap' in adata.obsm:
            sc.pl.umap(adata, color="celltype", title="UMAP by CellType",
                      legend_loc="on data", show=False)
            plt.savefig("{}/{}umap_celltype.png".format(outpath, prefix),
                       dpi=dpi, bbox_inches='tight')
            plt.savefig("{}/{}umap_celltype.pdf".format(outpath, prefix),
                       bbox_inches='tight')
            plt.close()
            print("Plot: {}umap_celltype.png(pdf) generated successfully!".format(prefix))
        else:
            print("Warning: X_umap not found in adata.obsm, skipping UMAP plot")

        # Prepare markers for dotplot grouped by celltype
        dedup_df = annotation_df.drop_duplicates(subset=["CellType", "Markers"])

        # Sort by cell type order if provided
        if celltype_order:
            dedup_df["CellType"] = pd.Categorical(dedup_df["CellType"], categories=celltype_order)
            dedup_df = dedup_df.sort_values("CellType").dropna(subset=["CellType"])
        else:
            dedup_df = dedup_df.sort_values("CellType")

        # Expand markers and deduplicate
        marker_long = []
        for _, row in dedup_df.iterrows():
            ct = row["CellType"]
            markers = [g.strip() for g in row["Markers"].split(",")]
            for g in markers:
                marker_long.append({"CellType": ct, "Markers": g})

        marker_long_df = pd.DataFrame(marker_long)
        marker_long_df = marker_long_df.drop_duplicates(subset=["Markers"], keep="first")

        # Filter genes present in data
        missing_genes = [g for g in marker_long_df["Markers"].unique() if g not in genes_in_data]
        if missing_genes:
            print("[DotPlot] {} markers NOT found and will be removed:".format(len(missing_genes)))
            print("  ", missing_genes[:20], "..." if len(missing_genes) > 20 else "")

        marker_long_df = marker_long_df[marker_long_df["Markers"].isin(genes_in_data)]

        # Build markers list by celltype
        markers_list = {}
        for ct in dedup_df["CellType"].unique():
            genes = marker_long_df[marker_long_df["CellType"] == ct]["Markers"].tolist()
            if genes:
                markers_list[ct] = genes

        # Set celltype order for plotting
        celltype_levels = list(dedup_df["CellType"].unique())
        adata.obs["celltype"] = pd.Categorical(adata.obs["celltype"], 
                                               categories=celltype_levels[::-1])

        # Plot dotplot grouped by celltype
        if markers_list:
            # Build kwargs for dotplot (no dendrogram for annotated celltype)
            dotplot_kwargs = {
                "var_names": markers_list,
                "groupby": "celltype",
                "standard_scale": "var",
                "cmap": "Reds",
                "show": False,
                "dendrogram": False,
                "var_group_rotation": var_group_rotation
            }
            if dotplot_width is not None and dotplot_height is not None:
                dotplot_kwargs["figsize"] = (dotplot_width, dotplot_height)
            sc.pl.dotplot(adata, **dotplot_kwargs)
            plt.savefig("{}/{}dotplot_celltype.png".format(outpath, prefix),
                       dpi=dpi, bbox_inches='tight')
            plt.savefig("{}/{}dotplot_celltype.pdf".format(outpath, prefix),
                       bbox_inches='tight')
            plt.close()
            print("Plot: {}dotplot_celltype.png(pdf) generated successfully!".format(prefix))

        # Save annotated data
        output_file = "{}/{}annotated.h5ad".format(outpath, prefix)
        adata.write_h5ad(output_file, compression="gzip")
        print("File: {}annotated.h5ad saved successfully!\n".format(prefix))

    print("============================ All Done ============================\n")


if __name__ == "__main__":
    main()