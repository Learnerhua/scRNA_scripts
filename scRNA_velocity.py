#!/oldhome/ouyjh/miniforge3/envs/scanpy/bin/python
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import h5py,os,time,argparse,platform,sys
from datetime import datetime

scv.settings.verbosity = 1
scv.settings.set_figure_params(style = 'scvelo', facecolor='white', dpi=300)
plt.rcParams['savefig.transparent'] = False

def data_load(input_file, annotation_column, target_cell=None):
    '''
    Prepare data for scVelo analysis.
    '''
    print(f"\n>>> Reading {input_file}\n")
    adata = sc.read_h5ad(input_file)
    print("\nBasic information of the loaded adata:\n")
    print(adata)
    print(f"\nSelected annotation column: '{annotation_column}'\n")
    if target_cell is not None:
        print(f"\nTarget cell types:\n")
        for ct in target_cell:
            print(f"- {ct}")
        adata_subset = adata[adata.obs[annotation_column].isin(target_cell)].copy()
        print("\nBasic information of the subset adata:\n")
        print(adata_subset)
    else:
        print("\nNo target cell types specified, using all cells.\n")
        adata_subset = adata.copy()
    print(adata_subset.obs[annotation_column].value_counts())
    return adata_subset

def preprocess_data(adata, target_cell=None):
    if target_cell is not None:
        print(f"\n>>> Preprocessing subset data with target cell types: {target_cell}\n")
        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=20,
            n_top_genes=2000
        )
        sc.tl.pca(adata, n_comps=30)
        sc.external.pp.harmony_integrate(adata, key='orig.ident')
        sc.pp.neighbors(
            adata,
            n_neighbors=30,
            n_pcs=30,
            use_rep='X_pca_harmony'
        )
        scv.pp.moments(adata)
        sc.tl.umap(adata)
    else:
        print(f"\n>>> Preprocessing full data\n")
        scv.pp.filter_and_normalize(
            adata,
            min_shared_counts=20,
            n_top_genes=2000
        )
        sc.pp.neighbors(
            adata,
            n_neighbors=30,
            n_pcs=30,
            use_rep='X_harmony'
        )
        scv.pp.moments(adata)
    return adata

def run_velocity_stochastic(adata, threads):
    """
    Run RNA velocity analysis using the stochastic model.
    """
    scv.tl.velocity(
        adata,
        mode="stochastic"
    )
    scv.tl.velocity_graph(adata, threads)
    return adata

def run_velocity_dynamical(adata, threads):
    """
    Run RNA velocity analysis using the dynamical model.
    """
    # 1. Recover gene-specific dynamics (mandatory for dynamical model)
    scv.tl.recover_dynamics(
        adata,
        n_jobs=threads
    )
    # 2. Compute RNA velocity using the dynamical model
    scv.tl.velocity(
        adata,
        mode="dynamical"
    )
    # 3. Construct velocity graph for downstream visualization
    scv.tl.velocity_graph(
        adata,
        n_jobs=threads
    )
    # 4. Compute latent time (core output of the dynamical model)
    scv.tl.latent_time(adata)
    return adata

def visualize_stochastic(adata, annotation_column, outpath, prefix, plot_n):
    """
    Visualize RNA velocity results using stochastic model.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=annotation_column,
        legend_loc='right margin',
        title="RNA Velocity Stream Plot",
        dpi = 300,
        ax = ax1,
        show=False
    )
    scv.tl.velocity_pseudotime(adata)
    scv.pl.scatter(
        adata,
        color="velocity_pseudotime",
        cmap="gnuplot",
        size=30,
        title="RNA Velocity Pseudotime Plot",
        ax = ax2,
        dpi=300,
        show=False
    )
    plt.tight_layout()
    plt.subplots_adjust(right=0.95)
    plt.savefig(f"{outpath}/{prefix}stochastic_stream_pseudotime.png", dpi=300)
    plt.savefig(f"{outpath}/{prefix}stochastic_stream_pseudotime.svg", bbox_inches="tight")
    plt.close()
    plot_n += 1
    print(f"Plot {plot_n}: {prefix}stochastic_stream_pseudotime.png (svg) has been generated successfully !\n")
    # scv.pl.velocity_embedding(
    #     adata,
    #     basis="umap",
    #     color=annotation_column,
    #     legend_loc='right margin',
    #     title="RNA Velocity Cell Plot",
    #     save=f"{outpath}/{prefix}stochastic_cell.png",
    #     dpi = 300,
    #     figsize = (8, 6),
    #     arrow_size=2,
    #     arrow_length=3,
    #     show=False
    # )
    # scv.pl.velocity_graph(
    #     adata,
    #     basis="umap",
    #     color=annotation_column,
    #     threshold=0.2,
    #     alpha=0.5,
    #     save=f"{outpath}/{prefix}stochastic_graph.png",
    #     dpi=300,
    #     figsize=(8, 6),
    #     show=False
    # )

def visualize_dynamical(adata, annotation_column, outpath, prefix, plot_n):
    """
    Visualize RNA velocity results using dynamical model.
    """
    # dynamical stream plot
    fig, ax = plt.subplots()
    scv.pl.velocity_embedding_stream(
        adata,
        title="RNA Velocity Stream Plot (Dynamical)",
        basis="umap",
        color=annotation_column,
        ax=ax,
        legend_loc="right margin",
        show=False
    )
    plt.tight_layout()
    plt.savefig(f"{outpath}/{prefix}dynamical_stream.png", dpi=300)
    plt.savefig(f"{outpath}/{prefix}dynamical_stream.svg", bbox_inches="tight")
    plt.close(fig)
    plot_n += 1
    print(f"Plot {plot_n}: {prefix}dynamical_stream.png (svg) has been generated successfully !\n")
    # dynamical gene velocity plot
    high_quality_genes = adata.var[adata.var['fit_likelihood'] > 0.1].index
    scv.tl.rank_velocity_genes(adata, groupby=annotation_column, min_corr=0.1)
    velocity_genes_df = pd.DataFrame(adata.uns['rank_velocity_genes']['names'])
    top_velocity_genes = velocity_genes_df.iloc[:3].values.flatten().tolist()
    top_velocity_genes = list(pd.unique(top_velocity_genes))
    combined_genes = set(high_quality_genes).intersection(top_velocity_genes)
    scv.pl.velocity(
        adata,
        var_names=list(combined_genes),
        color=annotation_column,
        show=False,
        save=f"{outpath}/{prefix}dynamical_genes.png",
        dpi=300
    )
    scv.pl.velocity(
        adata,
        var_names=list(combined_genes),
        color=annotation_column,
        show=False,
        save=f"{outpath}/{prefix}dynamical_genes.svg"
    )
    plot_n += 1
    print(f"Plot {plot_n}: {prefix}dynamical_genes.png (svg) has been generated successfully !\n")
    # dynamical latent time plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    scv.pl.scatter(
        adata,
        basis='umap',
        color=annotation_column,
        ax=ax1,
        title='Annotation distribution Plot',
        show=False
    )
    scv.pl.scatter(
        adata,
        basis='umap',
        color='latent_time',
        cmap='gnuplot',
        ax=ax2,
        title='Latent Time Plot',
        legend_loc='right margin',
        show=False
    )
    plt.tight_layout()
    plt.subplots_adjust(right=0.95)
    plt.savefig(f"{outpath}/{prefix}latent_time_umap.png", dpi=300)
    plt.savefig(f"{outpath}/{prefix}latent_time_umap.svg", bbox_inches="tight")
    plt.close()
    plot_n += 1
    print(f"Plot {plot_n}: {prefix}latent_time_umap.png (svg) has been generated successfully !\n")
    # velocity pseudotime
    scv.tl.velocity_pseudotime(adata)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    scv.pl.scatter(
        adata,
        basis='umap',
        color=annotation_column,
        ax=ax1,
        title='Annotation distribution Plot',
        show=False
    )
    scv.pl.scatter(
        adata,
        basis='umap',
        color='velocity_pseudotime',
        cmap='gnuplot',
        ax=ax2,
        title='Velocity Pseudotime Plot',
        legend_loc='right margin',
        show=False
    )
    plt.tight_layout()
    plt.subplots_adjust(right=0.95)
    plt.savefig(f"{outpath}/{prefix}velocity_pseudotime_umap.png", dpi=300)
    plt.savefig(f"{outpath}/{prefix}velocity_pseudotime_umap.svg", bbox_inches="tight")
    plt.close()
    plot_n += 1
    print(f"Plot {plot_n}: {prefix}velocity_pseudotime_umap.png (svg) has been generated successfully !\n")

def save_adata(adata, outpath, prefix, mode, file_n):
    """
    Save the AnnData object to an h5ad file.
    """
    output_file = f"{outpath}/{prefix}{mode}_velocity_result.h5ad"
    adata.write(output_file, compression='gzip')
    file_n += 1
    print(f"File {file_n}: {output_file}  saved successfully !\n")

def main():
    parser = argparse.ArgumentParser(
    description="A script for RNA velocity analysis using scVelo."
    "Example: scRNA_velocity -I scRNA_velocity_input.h5ad -AC celltype -M dynamical -TC 'B_cell,T_cell,Macrophage' -T 40 -P test_ -O test"
    )
    parser.add_argument(
        "-I", "--input-file",
        type=str,
        required=True,
        help=(
            "Input AnnData (.h5ad) file for RNA velocity analysis. (Required)\n"
            "The AnnData object must contain the following precomputed elements:\n"
            "  - layers['spliced'] and layers['unspliced'] matrices\n"
            "  - UMAP coordinates in obsm['X_umap']\n"
            "  - Batch-corrected embedding in obsm['X_harmony'] (used for neighbor graph construction)"
        )
    )
    parser.add_argument(
        "-AC", "--annotation-column",
        type=str,
        required=True,
        help="Column name in adata.obs used for cell annotation. (Required)"
    )
    parser.add_argument(
        "-M", "--mode",
        type=str,
        choices=["stochastic", "dynamical"],
        default="dynamical",
        help="Velocity model mode used by scVelo (default: dynamical)"
    )
    parser.add_argument(
        "-TC", "--target-cell",
        type=str,
        default=None,
        help=(
            "Target cell types for subsetting (comma-separated list). "
            "Example: -t 'B_cell,T_cell,Macrophage'. "
            "If not set, all cells are used."
            )
    )
    parser.add_argument(
        "-T", "--threads",
        type=int,
        default=4,
        help="Number of threads used for parallel computation (default: 4)"
    )
    parser.add_argument(
        "-P", "--prefix",
        type=str,
        default="",
        help="Prefix for output h5ad filenames (default: none)"
    )
    parser.add_argument(
        "-O", "--outpath",
        type=str,
        default=".",
        help="Output folder for processed h5ad files (default: current directory)"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
        
    args = parser.parse_args()
    input_file = args.input_file
    annotation_column = args.annotation_column
    mode = args.mode
    target_cell = args.target_cell.split(',') if args.target_cell is not None else None
    threads = args.threads
    prefix = args.prefix
    outpath = args.outpath
    file_n = 0
    plot_n = 0

    print("\n====================== scRNA Velocity Analysis ====================\n")
    print("Current time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform:", platform.platform())
    print("Working directory:", os.getcwd())
    print("Output directory:", outpath)
    print("Copyright © 2026 OYJH. All Rights Reserved.")
    print("\n===================================================================\n")

    os.makedirs(outpath, exist_ok=True)
    print("Start running ......\n")
    adata = data_load(input_file, annotation_column, target_cell)
    adata = preprocess_data(adata, target_cell)
    if mode == "stochastic":
        adata = run_velocity_stochastic(adata, threads)
        visualize_stochastic(adata, annotation_column, outpath, prefix, plot_n)
    if mode == "dynamical":
        adata = run_velocity_dynamical(adata, threads)
        visualize_dynamical(adata, annotation_column, outpath, prefix, plot_n)
    save_adata(adata, outpath, prefix, mode, file_n)

if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print("\nProgram finished\nRunning time:%.2fs"%(end - start))
    print("\n============================= ALL DONE ============================\n")