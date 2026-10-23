#!/usr/bin/env python3
import sys,argparse,warnings,platform,os
from datetime import datetime
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def is_raw_counts(adata_obj, sample=1000, float_tol=1e-6):
    """
    Heuristic check whether an AnnData's X looks like raw integer counts.

    For sparse X, .data may be stored as int or float; for dense X the dtype
    reflects what was written. We sample a subset of cells, then test whether
    every value is within ``float_tol`` of an integer — this tolerates the
    common case of integer counts stored as float32 / float64.
    """
    X = adata_obj.X
    n_cells = X.shape[0]
    if n_cells == 0:
        return {"is_integer": False, "zero_fraction": 1.0,
                "max_value": 0.0, "min_value": 0.0}

    idx = np.random.choice(n_cells, min(sample, n_cells), replace=False)

    if sp.issparse(X):
        submatrix = X[idx, :].toarray()
    else:
        submatrix = np.asarray(X[idx, :])

    is_integer = bool(np.all(np.abs(submatrix - np.round(submatrix)) < float_tol))

    return {
        "is_integer": is_integer,
        "zero_fraction": float(np.sum(submatrix == 0) / submatrix.size),
        "max_value": float(submatrix.max()),
        "min_value": float(submatrix.min()),
    }


def merge_duplicate_genes(X, genes):
    """
    Collapse duplicate gene names by summing expression across rows.

    Returns ``(X_new, unique_genes)``. Uses a sparse COO + ``sum_duplicates``
    path so we never materialise a dense matrix in memory, and casts the
    accumulator to int64 when the input is integer to avoid int32 overflow
    when many duplicates sum together.
    """
    unique_genes, gene_index = np.unique(np.asarray(genes), return_inverse=True)
    n_genes_new = len(unique_genes)

    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    coo = X.tocoo()
    new_col = gene_index[coo.col]
    X_new = sp.coo_matrix(
        (coo.data, (coo.row, new_col)),
        shape=(X.shape[0], n_genes_new),
    )
    X_new.sum_duplicates()

    # Promote integer dtypes to int64 so the column-sum can't overflow.
    if np.issubdtype(X.dtype, np.integer) and X_new.dtype != np.int64:
        X_new = X_new.astype(np.int64)

    return X_new.tocsr(), unique_genes


def main():
    parser = argparse.ArgumentParser(
        description="Process h5ad files by specifying the gene name column (-V) "
                    "and ensuring the expression matrix is raw counts "
                    "(prefer adata.raw.X if available)."
    )

    parser.add_argument(
        "-I", "--input_dir",
        help="Directory containing h5ad files (mutually exclusive with -F)"
    )
    parser.add_argument(
        "-F", "--file",
        help="Process a single h5ad file (mutually exclusive with -I)"
    )
    parser.add_argument(
        "-O", "--outpath",
        type=str,
        default=".",
        help="Output folder for processed h5ad files (default: current directory)"
    )
    parser.add_argument(
        "-P", "--prefix",
        type=str,
        default="",
        help="Prefix for output h5ad filenames (default: none)"
    )
    parser.add_argument(
        "-V", "--var_index",
        help="Column in adata.var to use as index"
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if (args.input_dir is None) == (args.file is None):
        parser.error("Please provide either -I (input_dir) or -F (single file), not both/neither.")

    if args.file is not None:
        if not os.path.isfile(args.file):
            parser.error(f"Input file not found: {args.file}")
    else:
        if not os.path.isdir(args.input_dir):
            parser.error(f"Input directory not found: {args.input_dir}")

    input_dir = args.input_dir
    outpath = args.outpath
    prefix = args.prefix
    var_index = args.var_index

    print("\n========================= scRNA H5ad Process ======================\n")
    print("Current time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform:", platform.platform())
    print("Working directory:", os.getcwd())
    print("Output directory:", outpath)
    print("Copyright © 2025 OYJH. All Rights Reserved.")
    print("\n===================================================================\n")

    os.makedirs(outpath, exist_ok=True)

    # Build the list of (fname, file_path) pairs to process
    if args.file is not None:
        files_to_process = [(os.path.basename(args.file), os.path.abspath(args.file))]
    else:
        files_to_process = [
            (fname, os.path.join(input_dir, fname))
            for fname in os.listdir(input_dir)
            if fname.lower().endswith(".h5ad")
        ]

    if not files_to_process:
        print(f"No .h5ad files found in {input_dir}. Nothing to do.")
        return

    for fname, file_path in files_to_process:
        print(f"\n>>> Reading {file_path}\n")
        adata = sc.read_h5ad(file_path)
        print("\nBasic information:\n")
        print(adata)
        print("\nMetadata (obs), top 3 rows:\n")
        print(adata.obs.head(n=3))
        print("\nGene info (var), top 3 rows:\n")
        print(adata.var.head(n=3))

        if var_index is None:
            print(
                "\nWarning: You did not specify -V (the column in adata.var to use as index). "
                "The original index of adata.var will be used instead.\n"
            )
        elif var_index not in adata.var.columns:
            raise ValueError(f"\nColumn {var_index} not found in adata.var for file {fname}\n")

        # Choose the count source. var_for_genes must track the source of X so
        # that gene names stay aligned with the columns we sum over.
        X_info = is_raw_counts(adata, sample=1000)
        print(f"[counts check] adata.X  -> is_integer={X_info['is_integer']}, "
              f"zero_fraction={X_info['zero_fraction']:.4f}, "
              f"min={X_info['min_value']:.4g}, max={X_info['max_value']:.4g}")
        if X_info["is_integer"]:
            print(f"[counts source] adata.X selected as count matrix "
                  f"(shape={adata.X.shape}, dtype={adata.X.dtype}, "
                  f"sparse={sp.issparse(adata.X)}).\n")
            X = adata.X
            var_for_genes = adata.var
        elif adata.raw is not None:
            raw_info = is_raw_counts(adata.raw, sample=1000)
            print(f"[counts check] adata.raw.X -> is_integer={raw_info['is_integer']}, "
                  f"zero_fraction={raw_info['zero_fraction']:.4f}, "
                  f"min={raw_info['min_value']:.4g}, max={raw_info['max_value']:.4g}")
            if raw_info["is_integer"]:
                print(f"[counts source] adata.raw.X selected as count matrix "
                      f"(shape={adata.raw.X.shape}, dtype={adata.raw.X.dtype}, "
                      f"sparse={sp.issparse(adata.raw.X)}). "
                      "adata.X is being replaced.\n")
                X = adata.raw.X
                var_for_genes = adata.raw.var
            else:
                raise ValueError(
                    f"\nFile {fname}: neither adata.X nor adata.raw.X looks like "
                    "integer counts. No raw counts found. Aborting.\n"
                )
        else:
            raise ValueError(
                f"\nFile {fname}: adata.X does not look like integer counts and "
                "adata.raw is None. No raw counts found. Aborting.\n"
            )

        if var_index is not None:
            var_for_genes = var_for_genes.copy()
            var_for_genes[var_index] = var_for_genes[var_index].astype(str)
            var_for_genes = var_for_genes.set_index(var_index)
            genes = var_for_genes.index.astype(str).values

            n_dup = int(len(genes) - len(np.unique(genes)))
            print(f"Input genes: {len(genes)} | Unique genes: {len(np.unique(genes))} "
                  f"| Duplicates to merge: {n_dup}")
            X_new, unique_genes = merge_duplicate_genes(X, genes)

            obs = adata.obs.copy()
            adata_new = ad.AnnData(
                X=X_new,
                obs=obs,
                var=pd.DataFrame(index=pd.Index(unique_genes, name=None)),
            )
            adata_new.layers["counts"] = X_new
        else:
            obs = adata.obs.copy()
            var = var_for_genes.copy()
            adata_new = ad.AnnData(X=X, obs=obs, var=var)
            adata_new.layers["counts"] = X

        # Save processed h5ad
        out_name = f"{prefix}{os.path.splitext(fname)[0]}.h5ad"
        out_path = os.path.join(outpath, out_name)
        adata_new.write_h5ad(out_path, compression="gzip")
        print(f"\nFile : {out_name} saved successfully !\n")


if __name__ == "__main__":
    main()