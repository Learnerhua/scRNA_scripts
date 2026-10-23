#!/oldhome/ouyjh/miniforge3/envs/old_base/bin/python
"""
Single-cell compositional (differential abundance) analysis with scCODA.

scCODA identifies cell-type composition shifts between conditions with a
Bayesian Dirichlet-multinomial model (Buettner, Ostner et al., Nat Commun 2021).
Because per-sample cell counts are compositional, a rise in one cell type
mechanically depresses the others -- a plain per-cell differential abundance
test would call most of those shifts significant. scCODA models that constraint,
picks one reference cell type, and uses a spike-and-slab prior to select credible
effects at a target FDR.

Accepted inputs (auto-detected by default):
  h5ad   : an AnnData .h5ad file produced from the annotated object. Only
           adata.obs is read (backed mode), so the expression matrix is never
           loaded and even a large file is fast.
  counts : one row per sample, one column per cell type, plus covariate column(s)
  obs    : one row per CELL, with a sample column and a cell-type column;
           counts are aggregated with a crosstab

Example usage:
  scRNA_composition_sccoda.py -I 05.seu_annotated.h5ad -S sample -CT celltype \\
      -GC group -O 08_composition -P 08.
  scRNA_composition_sccoda.py -I counts.tsv -GC group -O 08_composition -P 08.
  scRNA_composition_sccoda.py -I obs.tsv -T obs -S sample -CT celltype \\
      -GC group -O 08_composition -P 08.

Note: a Seurat .rds cannot be read directly by Python (S4 object); convert it to
.h5ad first (e.g. with sceasy / zellkonverter / SeuratDisk in R).
"""
import sys, argparse, warnings, os, logging, time, pickle, re, shutil
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('fontTools').setLevel(logging.WARNING)
# Quiet TensorFlow / TFP C++ chatter (scCODA samples with TFP).
os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")


def read_table(path):
    """Read a CSV/TSV by extension, falling back to delimiter sniffing."""
    sep = "\t" if path.lower().endswith((".tsv", ".txt", ".tab")) else ","
    try:
        df = pd.read_csv(path, sep=sep)
    except Exception:
        df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] == 1:                      # wrong delimiter -> retry
        df = pd.read_csv(path, sep=None, engine="python")
    return df


def build_counts(df, args):
    """
    Return (counts_df, obs_df) where counts_df is samples x cell types and
    obs_df holds one row per sample (index = sample names) with the covariate.

    counts mode: df is already samples x cell types (+ covariates); a dedicated
                 sample-name column, if present, becomes the index.
    obs mode   : df is one row per cell; counts are built with crosstab.
    """
    if args.input_type == "obs":
        for col in (args.sample_col, args.celltype_col):
            if col not in df.columns:
                sys.exit("Error: column '{}' not found. Available: {}".format(
                    col, list(df.columns)))
        sample = df[args.sample_col].astype(str)
        celltype = df[args.celltype_col].astype(str)
        counts = pd.crosstab(sample, celltype)

        obs = pd.DataFrame(index=counts.index)
        if args.group_col and args.group_col in df.columns:
            g = df[[args.sample_col, args.group_col]].drop_duplicates()
            g[args.sample_col] = g[args.sample_col].astype(str)
            # A sample must belong to exactly one group.
            if g[args.sample_col].duplicated().any():
                bad = g[args.sample_col][g[args.sample_col].duplicated()].unique()
                sys.exit("Error: sample(s) map to more than one {} : {}".format(
                    args.group_col, list(bad)))
            obs[args.group_col] = g.set_index(args.sample_col)[args.group_col]
        return counts, obs

    # ---- counts mode -------------------------------------------------------
    work = df.copy()
    if args.sample_col and args.sample_col in work.columns:
        work = work.set_index(args.sample_col)
    work.index = work.index.astype(str)

    cov_cols = [c for c in ([args.group_col] + args.extra_covariates)
                if c and c in work.columns]
    cell_cols = [c for c in work.columns if c not in cov_cols]

    # Everything non-covariate is a cell type (scCODA's own convention), so
    # refuse anything that is obviously not a count column.
    non_numeric = [c for c in cell_cols
                   if not pd.api.types.is_numeric_dtype(work[c])]
    if non_numeric:
        sys.exit("Error: non-numeric column(s) would be treated as CELL TYPES: "
                 "{}\n  -> rename/remove them, or pass them as covariates via "
                 "-GC / -EC.".format(non_numeric))

    counts = work[cell_cols].astype(int)
    obs = work[cov_cols] if cov_cols else None
    return counts, obs


def make_covariate(counts, obs, args):
    """Assemble the covariate DataFrame scCODA needs (one row per sample)."""
    if obs is not None and args.group_col in getattr(obs, "columns", []):
        cov = obs[[args.group_col]].copy()
    elif args.group_col and args.group_col in counts.columns:
        cov = counts[[args.group_col]].copy()
        counts = counts.drop(columns=[args.group_col])
    else:
        sys.exit("Error: covariate column '{}' not found.{}".format(
            args.group_col,
            " Pass -GC with the column to compare groups on."
            if args.group_col else " Pass -GC (e.g. -GC group)."))

    cov_df = cov.copy()
    for c in args.extra_covariates:
        if c and c in counts.columns:
            cov_df[c] = counts[c]
    cov_df.index = cov_df.index.astype(str)
    if not counts.index.equals(cov_df.index):
        counts = counts.loc[cov_df.index]
    return counts, cov_df


def build_formula(args):
    """Patsy formula: explicit -F wins, else wrap the group column if -GR given."""
    if args.formula:
        return args.formula
    col = args.group_col
    if args.group_ref:
        return "C({}, Treatment('{}'))".format(col, args.group_ref)
    return col


def read_h5ad_obs(path):
    """
    Read ONLY the .obs table of an h5ad file.

    backed="r" keeps the expression matrix on disk, so even a multi-GB h5ad
    costs just the metadata in memory. Falls back to a full read when the file
    is stored in a way backed mode cannot open.
    """
    import anndata as ad
    try:
        adata = ad.read_h5ad(path, backed="r")
        n_obs, n_vars = adata.n_obs, adata.n_vars
        obs = adata.obs.copy()
        try:
            adata.file.close()
        except Exception:
            pass
    except Exception:
        adata = ad.read_h5ad(path)
        n_obs, n_vars = adata.n_obs, adata.n_vars
        obs = adata.obs.copy()
    return obs, n_obs, n_vars


def load_input(path, args):
    """
    Return (table, kind, note).

    h5ad files and per-cell obs tables both come back as kind='obs' (one row per
    cell); sample x cell-type count tables come back as kind='counts'.
    """
    ext = os.path.splitext(path)[1].lower()
    kind = args.input_type
    note = ""

    if kind == "auto" and ext == ".h5ad":
        kind = "h5ad"

    if kind == "h5ad":
        obs, n_obs, n_vars = read_h5ad_obs(path)
        note = "h5ad: {} cells x {} genes (expression matrix NOT loaded)".format(
            n_obs, n_vars)
        return obs, "obs", note

    df = read_table(path)
    if kind == "auto":
        has_sample = args.sample_col in df.columns
        looks_cell_level = has_sample and len(df) > 5 * max(
            1, df[args.sample_col].nunique())
        kind = "obs" if looks_cell_level else "counts"
        note = "auto-detected '{}' from columns".format(kind)
    return df, kind, note


def save_fig(path_png, path_pdf, dpi):
    plt.savefig(path_png, dpi=dpi, bbox_inches="tight")
    plt.savefig(path_pdf, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="""
        Compositional (differential abundance) analysis of single-cell data with
        scCODA. Accepts an AnnData .h5ad file (only .obs is read), a sample x
        cell-type count table, or a per-cell obs table.
        Example:
          scRNA_composition_sccoda.py -I 05.seu_annotated.h5ad -S sample \\
              -CT celltype -GC group -O 08_composition -P 08.
        """
    )

    # ---- Required ----------------------------------------------------------
    parser.add_argument("-I", "--input", required=True, type=str,
                        help="Input file: an AnnData .h5ad, a sample x cell-type "
                             "count table, or a per-cell obs table (see -T). "
                             "CSV/TSV delimiter auto-detected. (Required)")

    # ---- Input layout ------------------------------------------------------
    parser.add_argument("-T", "--input_type", type=str, default="auto",
                        choices=["auto", "h5ad", "counts", "obs"],
                        help="Input layout: 'h5ad' (AnnData file, only .obs is "
                             "read), 'counts' (samples x cell types), 'obs' (one row "
                             "per cell) or 'auto' to detect from the file/columns "
                             "(default: auto)")
    parser.add_argument("-S", "--sample_col", type=str, default="sample",
                        help="Sample column in adata.obs / the obs table "
                             "(default: sample)")
    parser.add_argument("-CT", "--celltype_col", type=str, default="celltype",
                        help="Cell-type column in adata.obs / the obs table "
                             "(default: celltype)")

    # ---- Model -------------------------------------------------------------
    parser.add_argument("-GC", "--group_col", type=str, default="group",
                        help="Covariate column defining the groups to compare "
                             "(default: group)")
    parser.add_argument("-EC", "--extra_covariates", type=str, default="",
                        help="Comma-separated additional covariate columns "
                             "(default: none)")
    parser.add_argument("-F", "--formula", type=str, default="",
                        help="Patsy formula, overrides -GC/-GR. Wrap odd column "
                             "names: Q(\"my col\"). (default: built from -GC)")
    parser.add_argument("-GR", "--group_ref", type=str, default="",
                        help="Reference level of the group column; scCODA reports "
                             "every OTHER level against it. Default: the first "
                             "level present.")
    parser.add_argument("-R", "--reference_cell_type", type=str, default="automatic",
                        help="Reference cell type: a cell-type name, a column "
                             "index, or 'automatic'. The default 'automatic' picks "
                             "the cell type with the lowest relative-abundance "
                             "dispersion among those present in at least "
                             "--ref_absence_threshold (95 percent) of samples. "
                             "For small sample sizes (N < 8) it is RECOMMENDED to "
                             "override this with a name such as 'T' or 'B' to "
                             "avoid automatic selecting a low-abundance, low-power "
                             "cell type. See the 'Choosing the reference cell type' "
                             "section in the documentation for the full decision "
                             "framework. Examples: 'T', 'B', 'Schwann', 'automatic'.")
    parser.add_argument("-RA", "--ref_absence_threshold", type=float, default=0.05,
                        help="Absence threshold used by the automatic reference "
                             "choice (default: 0.05)")

    # ---- Sampling ----------------------------------------------------------
    parser.add_argument("-M", "--method", type=str, default="hmc",
                        choices=["hmc", "nuts", "hmc_da"],
                        help="Sampler: hmc, nuts or hmc_da (HMC + a second HMC per "
                             "cell type for the differential-abundance form) "
                             "(default: hmc)")
    parser.add_argument("-N", "--num_results", type=int, default=20000,
                        help="MCMC chain length (default: 20000)")
    parser.add_argument("-B", "--num_burnin", type=int, default=5000,
                        help="Burn-in iterations (default: 5000)")
    parser.add_argument("-L", "--num_leapfrog_steps", type=int, default=10,
                        help="HMC leapfrog steps (default: 10; nuts ignores it)")
    parser.add_argument("-TD", "--max_tree_depth", type=int, default=10,
                        help="NUTS max tree depth (default: 10; hmc ignores it)")
    parser.add_argument("-SS", "--step_size", type=float, default=0.01,
                        help="Sampler step size (default: 0.01)")
    parser.add_argument("-SD", "--seed", type=int, default=42,
                        help="Random seed for the sampler (default: 42)")
    parser.add_argument("--deterministic", dest="deterministic",
                        action="store_true",
                        help="Make the MCMC run bit-for-bit reproducible. "
                             "scCODA never passes a seed to tfp.mcmc.sample_chain "
                             "and TensorFlow's float kernels are not deterministic "
                             "by default, so identical settings otherwise give "
                             "slightly different posteriors. This flag (a) sets "
                             "TF_DETERMINISTIC_OPS=1 before TensorFlow is imported "
                             "and (b) injects -SD into tfp.mcmc.sample_chain. "
                             "Costs some speed; use it when you need to reproduce "
                             "an exact run.")
    parser.add_argument("-NT", "--num_threads", type=int, default=0,
                        help="CPU threads for TensorFlow / BLAS (default: 0 = let "
                             "TensorFlow decide). Must be set BEFORE TensorFlow is "
                             "imported, so it is applied via OMP_NUM_THREADS / "
                             "TF_NUM_INTRAOP_THREADS / TF_NUM_INTEROP_THREADS. For "
                             "small problems fewer threads are often faster than "
                             "many -- try 1, 4, 8 and 48 on your machine.")
    parser.add_argument("--verbose_sampler", dest="verbose_sampler",
                        action="store_true", default=False,
                        help="Print the sampler's per-step progress (default: off; "
                             "at 20000+ steps this floods the log)")

    # ---- Results -----------------------------------------------------------
    parser.add_argument("--fdr", type=float, default=0.05,
                        help="Target FDR for calling credible effects "
                             "(default: 0.05; the scCODA tutorial also uses 0.1-0.4)")
    parser.add_argument("--fdr_sweep", type=str, default="0.05,0.1,0.2,0.4",
                        help="Comma-separated FDR values to tabulate side by side "
                             "(default: 0.05,0.1,0.2,0.4); empty string to skip")

    # ---- Output ------------------------------------------------------------
    parser.add_argument("-P", "--prefix", type=str, default="",
                        help="Prefix for output files (default: none)")
    parser.add_argument("-O", "--outpath", type=str, default=".",
                        help="Output directory (default: current directory)")
    parser.add_argument("-W", "--width", type=float, default=8.0,
                        help="Figure width in inches (default: 8)")
    parser.add_argument("-H", "--height", type=float, default=6.0,
                        help="Figure height in inches (default: 6)")
    parser.add_argument("-dpi", "--dpi", type=int, default=300,
                        help="DPI for saved figures (default: 300)")
    parser.add_argument("--no-plot", dest="no_plot", action="store_true",
                        default=False,
                        help="Skip all figures, keep CSV/TXT/pickle outputs")
    parser.add_argument("-RS", "--rscript", type=str, default="",
                        help="Path to the Rscript binary used for the composition "
                             "sankey (Plot 1). Default: auto-detect -- 'Rscript' "
                             "on PATH, else alongside this Python interpreter "
                             "(i.e. the same conda env). Needs R + ggalluvial.")
    parser.add_argument("-ET", "--end_celltype", type=str, default="",
                        help="Cell type to place at the END of the ordering in the "
                             "composition sankey (optional). The ordering itself is "
                             "by global abundance, descending. (Named -ET here "
                             "because -EC is taken by --extra_covariates.)")
    parser.add_argument("--force_same_dir", dest="force_same_dir",
                        action="store_true", default=False,
                        help="Allow writing outputs into the input file's directory "
                             "(refused by default: source data dirs are read-only)")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.extra_covariates = [c.strip() for c in args.extra_covariates.split(",")
                             if c.strip()]

    input_file = args.input
    outpath = args.outpath
    prefix = args.prefix
    dpi = args.dpi

    # ------------------------------------------------------------------ #
    # Setup environment                                                    #
    # ------------------------------------------------------------------ #
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans",
                                       "Noto Sans CJK SC"]
    mpl.rcParams["axes.unicode_minus"] = False

    # Thread count must be pinned BEFORE TensorFlow is imported (scCODA pulls it
    # in). For tiny models like 7 samples x 18 cell types, oversubscription can
    # be SLOWER than a few threads -- measure rather than assume.
    if args.num_threads and args.num_threads > 0:
        n_threads = str(int(args.num_threads))
        for var in ("OMP_NUM_THREADS", "TF_NUM_INTRAOP_THREADS",
                    "TF_NUM_INTEROP_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
            os.environ[var] = n_threads
        print(">>> Thread limit: {} (OMP/TF intra+inter/BLAS)\n".format(n_threads))
    else:
        print(">>> Thread limit: not set (TensorFlow decides)\n")

    # Deterministic kernels must also be requested BEFORE TensorFlow is imported.
    # Without this, TF's CPU float kernels (oneDNN/Eigen) may reduce in a
    # different order on each run, so even a fully seeded sampler drifts.
    if args.deterministic:
        os.environ["TF_DETERMINISTIC_OPS"] = "1"
        print(">>> Determinism: TF_DETERMINISTIC_OPS=1\n")

    from sccoda.util import cell_composition_data as dat
    from sccoda.util import comp_ana as mod
    from sccoda.util import data_visualization as viz

    if not os.path.exists(input_file):
        sys.exit("Error: input file not found: {}".format(input_file))

    # Guard: never write into the directory that holds the source data.
    in_dir = os.path.dirname(os.path.abspath(input_file)) or os.getcwd()
    out_dir = os.path.abspath(outpath)
    if os.path.normpath(in_dir) == os.path.normpath(out_dir) and not args.force_same_dir:
        sys.exit("Error: refusing to write outputs into the input directory:\n"
                 "  {}\n  Use -O <separate analysis dir>, or pass --force_same_dir "
                 "to override.".format(in_dir))
    os.makedirs(outpath, exist_ok=True)

    log_file = os.path.join(outpath, "{}Composition_sccoda_summary.txt".format(prefix))
    logf = open(log_file, "w", encoding="utf-8")

    def write_summary(*parts):
        """Write one entry to the summary (strings only, like the R scripts)."""
        for p in parts:
            if not isinstance(p, str):
                raise TypeError("write_summary() accepts strings only")
        logf.write("".join(parts))

    # ------------------------------------------------------------------ #
    # Run info                                                             #
    # ------------------------------------------------------------------ #
    banner = "\n==================== scRNA Composition (scCODA) ====================\n"
    print(banner)
    write_summary(banner)
    # Pre-declare so the Reference celltype line can always be written, even when
    # the value is only known later (after CompositionalAnalysis resolves an
    # 'automatic' request to a concrete cell-type name).
    actual_reference = args.reference_cell_type
    for line in [
        "Current time      : " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Working directory : " + os.getcwd(),
        "Input file        : " + input_file,
        "Input type (req.) : " + args.input_type,
        "Output directory  : " + outpath,
        "Output prefix     : " + prefix,
        "Group column      : " + args.group_col,
        "Reference celltype: " + args.reference_cell_type
            + ("" if args.reference_cell_type != "automatic"
               else "  (will be resolved after model construction)"),
        "Formula           : " + build_formula(args),
        "Sampler           : {}  (N={}, burnin={}, seed={})".format(
            args.method, args.num_results, args.num_burnin, args.seed),
        "Target FDR        : {}".format(args.fdr),
        "",
        "Copyright (c) 2026 OYJH. All Rights Reserved.",
        "",
    ]:
        print("  " + line)
        write_summary(line + "\n")
    print("=" * 68 + "\n")
    write_summary("=" * 68 + "\n")

    run_start = time.time()

    # ------------------------------------------------------------------ #
    # Load data                                                            #
    # ------------------------------------------------------------------ #
    print(">>> Loading data\n")
    table, input_type, note = load_input(input_file, args)
    print("  Input type      : {}".format(input_type))
    if note:
        print("  Detail          : {}".format(note))
    print("  Table shape     : {} rows x {} columns\n".format(*table.shape))
    args.input_type = input_type

    counts, obs = build_counts(table, args)
    counts, cov = make_covariate(counts, obs, args)

    print(">>> Count table: {} samples x {} cell types\n".format(*counts.shape))
    print(counts.to_string(), "\n")

    total = int(counts.values.sum())
    # NOTE: keep this a pandas Series (indexed by sample name). Taking
    # counts.values.sum(axis=1) would give a positional numpy array, and the
    # loop below indexes it with the sample NAME -> IndexError.
    frac = counts.sum(axis=1) / total
    print("  Cells per sample :")
    for s in counts.index:
        print("    {:<20} {:>7d}  ({:.1f}%)".format(
            s, int(counts.loc[s].sum()), 100 * frac[s]))
    print("")

    # Empty samples cannot be modelled.
    empty = [s for s in counts.index if counts.loc[s].sum() == 0]
    if empty:
        sys.exit("Error: sample(s) with zero cells: {}".format(empty))

    write_summary("\n---------- Count table ({} samples x {} celltypes) ----------\n".format(
        *counts.shape))
    write_summary(counts.to_string() + "\n")
    write_summary("\nGroup membership:\n" + cov.to_string() + "\n")

    # ------------------------------------------------------------------ #
    # Build the scCODA data object                                         #
    # ------------------------------------------------------------------ #
    print(">>> Building scCODA data object\n")
    table = counts.copy()
    for c in cov.columns:
        table[c] = cov[c]
    data = dat.from_pandas(table, covariate_columns=list(cov.columns))
    print(data)
    print("")
    write_summary("\n---------- scCODA data object ----------\n{}\n".format(data))

    group_levels = list(pd.unique(cov[args.group_col].astype(str)))
    print("  Group level(s): {}\n".format(group_levels))
    if len(group_levels) < 2:
        sys.exit("Error: column '{}' has only one level: {}. Use -GC to pick the "
                 "column that separates the groups.".format(
                     args.group_col, group_levels))

    # Control vs case, derived from -GR (never hard-code group names).  Done
    # once here so every plot uses the same order.
    _ugroups = sorted(group_levels)
    if args.group_ref and args.group_ref in _ugroups:
        ctrl_lab = args.group_ref
        case_lab = [g for g in _ugroups if g != ctrl_lab][0]
    else:
        ctrl_lab, case_lab = _ugroups[0], _ugroups[-1]
    group_order = [ctrl_lab, case_lab]
    print("  Control group: {!r}   Case group: {!r}\n".format(ctrl_lab, case_lab))

    formula = build_formula(args)
    print(">>> Reference cell type: {}\n".format(args.reference_cell_type))
    if args.reference_cell_type == "automatic":
        print("  ('automatic' = lowest relative-abundance dispersion among cell "
              "types present in >= {:.0%} of samples)\n".format(
                  1 - args.ref_absence_threshold))

    # ------------------------------------------------------------------ #
    # Plot 0: relative-abundance dispersion (shows why the reference        #
    #         cell type is a sensible choice)                              #
    # ------------------------------------------------------------------ #
    # Keep scCODA's default plot (defaults for size/colours/sparse labels),
    # but tweak: (1) axis labels 2x + bold, (2) drop legend title,
    # (3) label every dot, not only the "abundant" ones.
    plot_start = time.time()
    if not args.no_plot:
        print(">>> Plot 0: relative abundance dispersion\n")
        try:
            import numpy as _np
            # label_cell_types=False so scCODA does not draw anything; we add
            # the labels ourselves below for every dot.
            ax0 = viz.rel_abundance_dispersion_plot(
                data, abundant_threshold=1 - args.ref_absence_threshold,
                figsize=(args.width, args.height), dpi=dpi,
                label_cell_types=False)
            # (1) axis labels: enlarged (1.4x) + bold
            big = {"fontsize": 1.4 * plt.rcParams["font.size"], "fontweight": "bold"}
            ax0.set_xlabel("Presence", **big)
            ax0.set_ylabel("Total dispersion", **big)
            # (3) label every dot, not just the abundant ones
            rel_abun = data.X / _np.sum(data.X, axis=1, keepdims=True)
            percent_zero = _np.sum(data.X == 0, axis=0) / data.X.shape[0]
            cell_type_disp = _np.var(rel_abun, axis=0) / _np.mean(rel_abun, axis=0)
            presence = 1 - percent_zero
            for ct, x, y in zip(data.var.index, presence, cell_type_disp):
                ax0.text(x + 0.01, y, str(ct), fontsize=8, alpha=0.85)
            # (4) dashed vertical line at the abundance threshold so the red/grey
            # boundary is visible (default threshold = 0.95 = 1 - RA).
            ax0.axvline(x=1 - args.ref_absence_threshold,
                        color="grey", linestyle="--", linewidth=1, alpha=0.7)
            # (2) drop the legend title (scCODA's default legend has it)
            leg = ax0.get_legend()
            if leg is not None:
                leg.set_title(None)
            save_fig(os.path.join(outpath, "{}composition_dispersion.png".format(prefix)),
                     os.path.join(outpath, "{}composition_dispersion.pdf".format(prefix)),
                     dpi)
            print("Plot 0: {}composition_dispersion.png(pdf) generated "
                  "successfully!\n".format(prefix))
        except Exception as exc:
            print("Plot 0 skipped ({}: {})\n".format(type(exc).__name__, exc))

    # ------------------------------------------------------------------ #
    # Plot 1: composition sankey (group x cell type, alluvial)              #
    # ------------------------------------------------------------------ #
    # Rendered by the companion R helper sankey_composition.r (same directory) with
    # ggalluvial, so the alluvial look stays faithful to the reference
    # implementation.  Requires R + ggalluvial; the figure is skipped with a
    # clear message if either is unavailable.
    if not args.no_plot:
        print(">>> Plot 1: composition sankey (alluvial)\n")
        try:
            import subprocess
            import tempfile

            # group x cell-type proportions, computed from the scCODA object
            # (rows = samples, columns = cell types, .obs carries the group).
            _counts_by_sample = pd.DataFrame(
                data.X, index=list(data.obs_names), columns=list(data.var_names))
            _grp_of_sample = data.obs[args.group_col].astype(str)
            _by_group = _counts_by_sample.groupby(_grp_of_sample.values).sum()
            _prop = _by_group.div(_by_group.sum(axis=1), axis=0)

            # Cell-type order: global abundance, descending (matches the
            # reference implementation); an optional -ET type is pushed last.
            _totals = _by_group.sum(axis=0)
            ct_order = [str(c) for c in _totals.sort_values(ascending=False).index]
            if args.end_celltype:
                if args.end_celltype in ct_order:
                    ct_order = ([c for c in ct_order if c != args.end_celltype]
                                + [args.end_celltype])
                else:
                    print("  warning: -ET {!r} not found among cell types; "
                          "ignored".format(args.end_celltype))
            print("  cell-type order (global abundance, desc): {}{}".format(
                ", ".join(ct_order[:5]),
                " ... (+{} more)".format(len(ct_order) - 5) if len(ct_order) > 5 else ""))

            _prop = _prop.reindex(columns=ct_order)
            _prop = _prop.reindex(index=group_order)
            _prop.index.name = args.group_col
            _prop = _prop.reset_index()

            r_helper = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "sankey_composition.r")
            if not os.path.exists(r_helper):
                raise FileNotFoundError("R helper not found: {}".format(r_helper))

            # Locate Rscript: explicit -RS wins, then PATH, then the sibling of
            # this Python interpreter (same conda env).
            rscript = args.rscript or shutil.which("Rscript")
            if not rscript:
                _cand = os.path.join(os.path.dirname(sys.executable), "Rscript")
                rscript = _cand if os.path.exists(_cand) else None
            if not rscript:
                raise FileNotFoundError(
                    "Rscript not found; install R or pass -RS /path/to/Rscript")

            tmpdir = tempfile.mkdtemp(prefix="sankey_")
            csv_path = os.path.join(tmpdir, "prop.csv")
            _prop.to_csv(csv_path, index=False)

            png_path = os.path.join(outpath, "{}composition_sankey.png".format(prefix))
            pdf_path = os.path.join(outpath, "{}composition_sankey.pdf".format(prefix))
            # The cell-type order travels via the CSV column order.
            cmd = [rscript, r_helper, csv_path, png_path, pdf_path,
                   str(args.width), str(args.height), ",".join(group_order)]
            proc = subprocess.run(cmd, capture_output=True, text=True)
            if proc.returncode != 0:
                raise RuntimeError(
                    "Rscript failed ({}): {}".format(
                        proc.returncode,
                        (proc.stderr or proc.stdout).strip().splitlines()[-1]
                        if (proc.stderr or proc.stdout).strip() else "no output"))

            shutil.rmtree(tmpdir, ignore_errors=True)
            print("Plot 1: {}composition_sankey.png(pdf) generated "
                  "successfully!\n".format(prefix))
        except Exception as exc:
            print("Plot 1 skipped ({}: {})\n".format(type(exc).__name__, exc))

    # ------------------------------------------------------------------ #
    # Plot 2: per-cell-type relative-abundance boxplots                    #
    # ------------------------------------------------------------------ #
    # Custom boxplot (replaces scCODA's default viz.boxplots).  The user
    # controls point size, box colour, figure width, and group order here.
    if not args.no_plot:
        print(">>> Plot 2: relative-abundance boxplots\n")

        # ctrl_lab / case_lab / group_order come from the hoisted
        # control-vs-case block above (derived from -GR).
        palette = {ctrl_lab: "#90BFF9", case_lab: "#FFA0A0"}  # user-fixed colours
        print("  ctrl = {!r} (color #90BFF9); case = {!r} (color #FFA0A0)"
              .format(ctrl_lab, case_lab))

        # -- long-format table: sample x celltype x group -> celltype x rel --
        # Re-read the per-cell obs (cheap: backed mode, .obs only -- no expression matrix).
        import anndata as _ad_for_box
        _ad_bp = _ad_for_box.read_h5ad(input_file, backed="r")
        _ct = pd.crosstab(_ad_bp.obs[args.sample_col], _ad_bp.obs[args.celltype_col])
        try: _ad_bp.file.close()
        except Exception: pass
        rel = (_ct.div(_ct.sum(axis=1), axis=0)
                 .reset_index()
                 .melt(id_vars=args.sample_col, var_name=args.celltype_col,
                       value_name="rel"))
        # attach the group column (one value per sample -- broadcast)
        # Re-read obs again to fetch the group column (lightweight, backed mode)
        _ad_g = _ad_for_box.read_h5ad(input_file, backed="r")
        _grp_map = _ad_g.obs.drop_duplicates(args.sample_col)[[args.sample_col, args.group_col]]
        try: _ad_g.file.close()
        except Exception: pass
        rel = rel.merge(_grp_map, on=args.sample_col)
        rel[args.group_col] = pd.Categorical(rel[args.group_col],
                                             categories=group_order, ordered=True)

        # -- auto figure width: 1 inch per cell type (>=20 inch) --
        n_types = rel[args.celltype_col].nunique()
        fig_w = max(20, 1.0 * n_types)
        fig_h = 6
        # honour user-supplied -W/-H if they are larger (user explicitly wanted bigger)
        fig_w = max(fig_w, args.width)
        fig_h = max(fig_h, args.height)

        import seaborn as sns
        sns.set_theme(style="white", context="notebook")
        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
        sns.boxplot(data=rel, x=args.celltype_col, y="rel", hue=args.group_col,
                    hue_order=group_order, palette=palette,
                    width=0.7, fliersize=0, ax=ax)
        # Use a black-only palette so the strip-plot dots are always solid
        # black regardless of the boxplot palette; hue + dodge still gives
        # the per-group offset.
        sns.stripplot(data=rel, x=args.celltype_col, y="rel", hue=args.group_col,
                      hue_order=group_order, palette=["black", "black"],
                      size=4.0, dodge=True, jitter=0.25,
                      edgecolor="none", alpha=0.85, ax=ax)
        # remove duplicate strip-legend entries; remove the legend title too
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[:2], labels[:2], title=None, frameon=False)
        # axis labels: enlarged (1.4x) + bold
        big = {"fontsize": 1.4 * plt.rcParams["font.size"], "fontweight": "bold"}
        ax.set_ylabel("Relative abundance", **big)
        ax.set_xlabel("Cell type",          **big)
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()

        png = os.path.join(outpath, "{}composition_boxplot.png".format(prefix))
        pdf = os.path.join(outpath, "{}composition_boxplot.pdf".format(prefix))
        fig.savefig(png, dpi=dpi, bbox_inches="tight")
        fig.savefig(pdf, bbox_inches="tight")
        plt.close(fig)
        print("Plot 2: {}composition_boxplot.png(pdf) generated "
              "successfully!\n".format(prefix))
    print("  Plotting took {:.1f} s\n".format(time.time() - plot_start))

    # ------------------------------------------------------------------ #
    # Model setup + inference                                              #
    # ------------------------------------------------------------------ #
    print(">>> Model setup and inference\n")
    print("  Formula        : {}".format(formula))
    print("  Reference type : {}".format(args.reference_cell_type))
    print("  Sampler        : {}\n".format(args.method))

    # Seeding MUST happen before CompositionalAnalysis is constructed: its
    # __init__ draws the MCMC starting values with tf.random.normal(...)
    # (scCODA_model.py, 'init_b_offset').  Setting the seed afterwards leaves
    # the initial chain state random, which is why identical runs used to drift.
    import tensorflow as tf
    # Best-effort: also pin threads through the TF API (works only if TF has not
    # created its thread pools yet -- the env vars above are the reliable path).
    if args.num_threads and args.num_threads > 0:
        try:
            tf.config.threading.set_intra_op_parallelism_threads(args.num_threads)
            tf.config.threading.set_inter_op_parallelism_threads(args.num_threads)
            print("  TF threads pinned via tf.config.threading: {}\n".format(
                args.num_threads))
        except Exception as exc:
            print("  (tf.config.threading not applied: {})\n".format(exc))
    tf.random.set_seed(args.seed)
    np.random.seed(args.seed)

    if args.deterministic:
        # scCODA calls tfp.mcmc.sample_chain WITHOUT a seed (tfp 0.25 does not
        # even accept seed= on HamiltonianMonteCarlo). tf.random.set_seed alone
        # does not pin the sampler, so inject the seed here.  Combined with
        # TF_DETERMINISTIC_OPS=1 this makes a run bit-for-bit reproducible.
        import tensorflow_probability as tfp
        if not getattr(tfp.mcmc.sample_chain, "_sccoda_seeded", False):
            _orig_sample_chain = tfp.mcmc.sample_chain

            def _seeded_sample_chain(*a, **kw):
                kw.setdefault("seed", [int(args.seed), int(args.seed)])
                return _orig_sample_chain(*a, **kw)

            _seeded_sample_chain._sccoda_seeded = True
            tfp.mcmc.sample_chain = _seeded_sample_chain
        print("  determinism: tfp.mcmc.sample_chain seeded with {}\n".format(args.seed))

    model = mod.CompositionalAnalysis(
        data, formula=formula,
        reference_cell_type=args.reference_cell_type,
        automatic_reference_absence_threshold=args.ref_absence_threshold)

    # After model construction, scCODA has resolved the 'automatic' choice to a
    # concrete cell-type name (or its index). Capture it so the summary file
    # reports what was ACTUALLY used, not just the user's request.
    if args.reference_cell_type == "automatic":
        try:
            chosen_name = str(data.var_names[model.reference_cell_type])
        except Exception:
            chosen_name = "index={}".format(model.reference_cell_type)
        actual_reference = chosen_name
        # Cross-check against any cell-type name the user passed
        if data.var_names[model.reference_cell_type] not in data.var_names:
            actual_reference = str(model.reference_cell_type)
        msg = "  (automatic resolved to: {})\n".format(actual_reference)
        print(msg)
        write_summary(msg)
    else:
        actual_reference = args.reference_cell_type

    sample_start = time.time()
    verbose = {"verbose": bool(args.verbose_sampler)}
    if args.method == "hmc":
        result = model.sample_hmc(num_results=args.num_results,
                                  num_burnin=args.num_burnin,
                                  num_leapfrog_steps=args.num_leapfrog_steps,
                                  step_size=args.step_size, **verbose)
    elif args.method == "hmc_da":
        result = model.sample_hmc_da(num_results=args.num_results,
                                     num_burnin=args.num_burnin,
                                     num_leapfrog_steps=args.num_leapfrog_steps,
                                     step_size=args.step_size, **verbose)
    else:
        result = model.sample_nuts(num_results=args.num_results,
                                   num_burnin=args.num_burnin,
                                   max_tree_depth=args.max_tree_depth,
                                   step_size=args.step_size, **verbose)
    print("\n  Sampling took {:.1f} s\n".format(time.time() - sample_start))

    # ------------------------------------------------------------------ #
    # Results                                                              #
    # ------------------------------------------------------------------ #
    print(">>> Summary (target FDR = {})\n".format(args.fdr))
    result.set_fdr(est_fdr=args.fdr)
    # NOTE: CAResult.summary() PRINTS to console and returns None -- never wrap
    # it in print() (that appends a stray "None") and never str() it to a file.
    result.summary()
    print("")

    # ---- effect / intercept tables ---------------------------------------
    intercept_df, effect_df = result.summary_prepare(est_fdr=args.fdr)
    intercept_csv = os.path.join(outpath, "{}composition_intercepts.csv".format(prefix))
    effect_csv = os.path.join(outpath, "{}composition_effects.csv".format(prefix))
    intercept_df.to_csv(intercept_csv)
    effect_df.to_csv(effect_csv)
    print("File 1: {} saved".format(os.path.basename(intercept_csv)))
    print("File 2: {} saved\n".format(os.path.basename(effect_csv)))

    # The summary file gets the tables themselves, since summary() only prints.
    write_summary("\n---------- Intercepts (FDR = {}) ----------\n".format(args.fdr))
    write_summary(intercept_df.to_string() + "\n")
    write_summary("\n---------- Effects (FDR = {}) ----------\n".format(args.fdr))
    write_summary(effect_df.to_string() + "\n")

    # ---- credible effects at the target FDR ------------------------------
    credible = result.credible_effects(est_fdr=args.fdr)
    print(">>> Credible effects (FDR = {})\n".format(args.fdr))
    print(credible.to_string())
    print("")
    credible_csv = os.path.join(outpath, "{}composition_credible_effects.csv".format(prefix))
    credible.rename("credible").to_frame().to_csv(credible_csv)
    print("File 3: {} saved\n".format(os.path.basename(credible_csv)))

    write_summary("\n---------- Credible effects (FDR = {}) ----------\n".format(args.fdr))
    write_summary(credible.to_string() + "\n")

    # ---- FDR sweep --------------------------------------------------------
    if args.fdr_sweep.strip():
        fdrs = [float(x) for x in args.fdr_sweep.split(",") if x.strip()]
        sweep = pd.DataFrame(index=credible.index)
        for f in fdrs:
            try:
                sweep["FDR={}".format(f)] = result.credible_effects(est_fdr=f)
            except Exception as exc:
                print("  FDR sweep: FDR={} skipped ({})".format(f, exc))
        sweep_csv = os.path.join(outpath, "{}composition_fdr_sweep.csv".format(prefix))
        sweep.to_csv(sweep_csv)
        print(">>> FDR sweep\n")
        print(sweep.to_string())
        print("\nFile 4: {} saved\n".format(os.path.basename(sweep_csv)))
        write_summary("\n---------- FDR sweep ----------\n" + sweep.to_string() + "\n")

    # ---- Plot 3: effect sizes (log2 fold change) -------------------------
    if not args.no_plot:
        print(">>> Plot 3: effect sizes\n")
        try:
            eff = effect_df.copy()
            # Locate columns by keyword so this survives small version changes.
            lfc_col = next((c for c in eff.columns if "log2" in str(c).lower()),
                           None)
            if lfc_col is None:
                lfc_col = next((c for c in eff.columns
                                if "final parameter" in str(c).lower()), None)
            lo_col = next((c for c in eff.columns
                           if str(c).upper().startswith("HDI") and "lower" in str(c).lower()),
                          None)
            hi_col = next((c for c in eff.columns
                           if str(c).upper().startswith("HDI") and "upper" in str(c).lower()),
                          None)
            if lfc_col is None:
                raise ValueError("no effect-size column found in {}".format(list(eff.columns)))

            # Build a plain "cell type -> credible?" map once, so the plotting
            # loop never depends on MultiIndex gymnastics.
            def cell_of(idx):
                return str(idx[-1]) if isinstance(idx, tuple) else str(idx)
            cred_map = {}
            for idx, val in credible.items():
                cred_map[cell_of(idx)] = bool(val)

            is_multi = isinstance(eff.index, pd.MultiIndex) and eff.index.nlevels >= 2
            groups = eff.groupby(level=0) if is_multi else [(None, eff)]
            for key, sub in groups:
                sub = sub.copy()
                sub.index = [cell_of(i) for i in sub.index]
                sub = sub.sort_values(lfc_col)

                y = np.arange(len(sub))
                fig, ax = plt.subplots(
                    figsize=(args.width, max(3.0, 0.32 * len(sub) + 1.5)), dpi=dpi)
                colors = ["#B22222" if cred_map.get(lab, False) else "#999999"
                          for lab in sub.index]

                if lo_col and hi_col:
                    xerr = np.vstack([(sub[lfc_col].values - sub[lo_col].values),
                                      (sub[hi_col].values - sub[lfc_col].values)])
                    ax.errorbar(sub[lfc_col].values, y, xerr=np.abs(xerr), fmt="none",
                                ecolor="#444444", elinewidth=0.8, capsize=2, zorder=1)
                ax.scatter(sub[lfc_col].values, y, c=colors, s=42, zorder=2)
                ax.axvline(0, color="black", linewidth=0.8, linestyle="--", zorder=0)
                ax.set_yticks(y)
                # bold y-axis tick labels (cell type names)
                ax.set_yticklabels(sub.index, fontweight="bold")
                # bold + 1.4x x-axis label
                big = {"fontsize": 1.4 * plt.rcParams["font.size"], "fontweight": "bold"}
                ax.set_xlabel(lfc_col, **big)
                ax.set_ylabel("")
                # drop the title (per user request)
                ax.set_title("")
                for sp in ("top", "right"):
                    ax.spines[sp].set_visible(False)
                fig.tight_layout()
                # Fixed filename: drop the suffix that encodes the covariate
                # level.  Multiple covariate levels are still saved separately
                # in script-level logic, but for the common single-level case
                # this gives a stable name.
                fig.savefig(os.path.join(outpath,
                            "{}composition_effects.png".format(prefix)), dpi=dpi)
                fig.savefig(os.path.join(outpath,
                            "{}composition_effects.pdf".format(prefix)),
                            bbox_inches="tight")
                plt.close(fig)
            print("Plot 3: {}composition_effects*.png(pdf) generated "
                  "successfully!\n".format(prefix))
        except Exception as exc:
            print("Plot 3 skipped ({}: {})\n".format(type(exc).__name__, exc))

    # ---- Save the sampler result -----------------------------------------
    pkl_file = os.path.join(outpath, "{}composition_result.pkl".format(prefix))
    result.save(pkl_file)
    print("File 5: {} saved".format(os.path.basename(pkl_file)))
    print("  (reload with:  import pickle; pickle.load(open(path,'rb')))\n")

    elapsed = time.time() - run_start
    print("Total time   : {:.1f} s\n".format(elapsed))
    write_summary("\nTotal time : {:.1f} s\n".format(elapsed))

    tail = "\n============================ All Done ============================\n"
    print(tail)
    write_summary(tail)
    logf.close()


if __name__ == "__main__":
    main()
