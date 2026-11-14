#!/oldhome/ouyjh/miniforge3/envs/pySCENIC/bin/python
import argparse, os, sys, subprocess, platform, time, pandas as pd
from datetime import datetime

REF_FILES = {
    "mouse": {
        "tfs": "mm_mgi_tfs.txt",
        "rankings": "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
        "annotations": "motifs-v9-nr.mgi-m0.001-o0.0.tbl",
    },
    "human": {
        "tfs": "hs_hgnc_tfs.txt",
        "rankings": "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
        "annotations": "motifs-v9-nr.hgnc-m0.001-o0.0.tbl",
    },
}

ALL_STEPS = ["grn", "ctx", "aucell"]


def find_pyscenic():
    # The CLI sits next to the interpreter that runs this script.
    cand = os.path.join(os.path.dirname(sys.executable), "pyscenic")
    if os.path.exists(cand):
        return cand
    raise ValueError(
        f"pyscenic CLI not found at {cand}.\n"
        f"This script must be run with the python of the pySCENIC environment; "
        f"the current interpreter is {sys.executable}."
    )


def check_ref_files(ref_dir, species):
    need = REF_FILES[species]
    paths = {}
    for key, fname in need.items():
        p = os.path.join(ref_dir, fname)
        if not os.path.exists(p):
            raise ValueError(f"Reference file not found: {p}")
        paths[key] = p
    return paths


def count_regulons(path):
    # `pyscenic ctx` writes a two-line multi-index header before the real header.
    if not os.path.exists(path):
        return 0
    try:
        return pd.read_csv(path, skiprows=[1, 2], index_col=[0, 1]).shape[0]
    except Exception:
        with open(path) as fh:
            return max(sum(1 for line in fh if line.strip()) - 3, 0)


def check_loom(path):
    import loompy
    with loompy.connect(path) as ds:
        if "Gene" not in ds.ra:
            raise ValueError(f"Loom has no row attribute 'Gene': {path}")
        if "CellID" not in ds.ca:
            raise ValueError(f"Loom has no column attribute 'CellID': {path}")
        return ds.shape


def run_cmd(cmd, step):
    print(f"\n>>> [{step}] {' '.join(cmd)}\n")
    t0 = time.time()
    ret = subprocess.run(cmd)
    if ret.returncode != 0:
        raise ValueError(f"pyscenic {step} failed with exit code {ret.returncode}")
    print(f">>> [{step}] finished in {time.time() - t0:.1f} s")


def step_grn(cli, loom, tfs, adj, threads, method, seed, force):
    if os.path.exists(adj) and not force:
        print(f">>> [grn] skip, exists: {adj}")
        return
    cmd = [cli, "grn", "-o", adj, "--num_workers", str(threads),
           "--method", method, "--seed", str(seed), loom, tfs]
    run_cmd(cmd, "grn")


def step_ctx(cli, adj, rankings, annotations, loom, regulons, threads, mask, force):
    if os.path.exists(regulons) and not force:
        print(f">>> [ctx] skip, exists: {regulons}")
        return
    cmd = [cli, "ctx", adj, rankings,
           "--annotations_fname", annotations,
           "--expression_mtx_fname", loom,
           "--num_workers", str(threads),
           "-o", regulons]
    if mask:
        cmd.append("--mask_dropouts")
    run_cmd(cmd, "ctx")


def step_aucell(cli, loom, regulons, out_loom, threads, force):
    if os.path.exists(out_loom) and not force:
        print(f">>> [aucell] skip, exists: {out_loom}")
        return
    cmd = [cli, "aucell", "--num_workers", str(threads),
           "-o", out_loom, loom, regulons]
    run_cmd(cmd, "aucell")


def main():
    parser = argparse.ArgumentParser(
        description="Run the pySCENIC pipeline (grn -> ctx -> aucell) on a loom file. "
        "Example: scRNA_pyscenic.py -I 01_pyscenic/01.expr_norm.loom -S mouse -RD ref_data -T 40 -P 01. -O 01_pyscenic"
    )
    parser.add_argument("-I", "--input-loom", required=True,
                        help="Input loom file built by scRNA_pyscenic_prepare.r. (Required)")
    parser.add_argument("-S", "--species", required=True, choices=["mouse", "human"],
                        help="Species; selects the TF list, motif rankings and annotation table. (Required)")
    parser.add_argument("-RD", "--ref-dir", default="ref_data",
                        help="Directory holding the three reference files. (Optional, default: ./ref_data)")
    parser.add_argument("-GM", "--grn-method", default="grnboost2", choices=["grnboost2", "genie3"],
                        help="GRN inference method. (Optional, default: grnboost2)")
    parser.add_argument("-ST", "--steps", default="all",
                        help="Comma-separated steps to run from grn,ctx,aucell, or 'all'. (Optional, default: all)")
    parser.add_argument("-T", "--threads", type=int, default=12,
                        help="Number of workers passed to pyscenic. (Optional, default: 12)")
    parser.add_argument("-SE", "--seed", type=int, default=777,
                        help="Random seed for the GRN step. (Optional, default: 777)")
    parser.add_argument("--mask-dropouts", dest="mask_dropouts",
                        action="store_true", default=True,
                        help="Pass --mask_dropouts to ctx. (Default: enabled; this flag is redundant)")
    parser.add_argument("--no-mask-dropouts", dest="mask_dropouts", action="store_false",
                        help="Do not pass --mask_dropouts to ctx.")
    parser.add_argument("-F", "--force", action="store_true", default=False,
                        help="Re-run steps even if their output file already exists. (Default: off)")
    parser.add_argument("-P", "--prefix", default="",
                        help="Prefix for output files. (Optional, default: none)")
    parser.add_argument("-O", "--outpath", default=".",
                        help="Output directory. (Optional, default: current directory)")
    args = parser.parse_args()

    input_loom = args.input_loom
    species = args.species
    ref_dir = args.ref_dir
    grn_method = args.grn_method
    threads = args.threads
    seed = args.seed
    mask_dropouts = args.mask_dropouts
    force = args.force
    prefix = args.prefix
    outpath = args.outpath

    if args.steps.strip().lower() == "all":
        steps = list(ALL_STEPS)
    else:
        steps = [s.strip().lower() for s in args.steps.split(",") if s.strip()]
        bad = [s for s in steps if s not in ALL_STEPS]
        if bad:
            raise ValueError(f"Unknown step(s): {bad}; allowed: {ALL_STEPS} or 'all'")
    if threads < 1:
        raise ValueError("-T/--threads must be >= 1")
    if not os.path.exists(input_loom):
        raise ValueError(f"Input loom not found: {input_loom}")

    print("\n======================= scRNA pySCENIC ========================\n")
    print("Current time:", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    print("Operating system:", platform.system(), platform.release())
    print("Platform:", platform.platform())
    print("Working directory:", os.getcwd())
    print("Output directory:", outpath)
    print("Species:", species)
    print("Steps:", ",".join(steps))
    print("Copyright © 2026 OYJH. All Rights Reserved.")
    print("\n===============================================================\n")

    os.makedirs(outpath, exist_ok=True)

    cli = find_pyscenic()
    print(">>> pyscenic CLI:", cli)
    loom_shape = check_loom(input_loom)
    print(f">>> Input loom: {loom_shape[0]} genes x {loom_shape[1]} cells")
    refs = check_ref_files(ref_dir, species)
    print(">>> Reference directory:", ref_dir)

    adj = os.path.join(outpath, prefix + "adj.csv")
    regulons = os.path.join(outpath, prefix + "regulons.csv")
    out_loom = os.path.join(outpath, prefix + "SCENIC.loom")

    if "grn" in steps:
        existed = os.path.exists(adj)
        step_grn(cli, input_loom, refs["tfs"], adj, threads, grn_method, seed, force)
        if not existed or force:
            print(f"\nFile 1 : {adj} saved successfully !\n")
    else:
        if not os.path.exists(adj) and "ctx" in steps:
            raise ValueError(f"ctx needs {adj}; run the grn step first.")

    if "ctx" in steps:
        existed = os.path.exists(regulons)
        step_ctx(cli, adj, refs["rankings"], refs["annotations"], input_loom,
                 regulons, threads, mask_dropouts, force)
        n_reg = count_regulons(regulons)
        if not existed or force:
            print(f"\nFile 2 : {regulons} saved successfully ! ({n_reg} regulon rows)\n")
        else:
            print(f"\n>>> [ctx] reusing existing {regulons} ({n_reg} regulon rows)\n")
        if n_reg == 0:
            print("!!! No regulon passed the motif enrichment; aucell will have nothing to score.")

    if "aucell" in steps:
        if not os.path.exists(regulons):
            raise ValueError(f"aucell needs {regulons}; run the ctx step first.")
        existed = os.path.exists(out_loom)
        step_aucell(cli, input_loom, regulons, out_loom, threads, force)
        if not existed or force:
            print(f"\nFile 3 : {out_loom} saved successfully !\n")

    print(">>> Done. Results in", outpath)
    print("\n=========================All Done !============================\n")


if __name__ == "__main__":
    main()
