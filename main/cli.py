__author__ = "Maulik Nariya" 
__date__ = "November 2025"
__copyright__ = "MIT license"

"""
Sanity single-cell RNA-seq wrapper — CLI

This script defines command-line parsing for the Sanity wrapper.
It detects input mode, input/output files and directories, and minimal preprocessing parameters.
"""

import argparse
import multiprocessing as mp
import os

SUPPORTED_MODES = ("cellranger", "anndata", "matrix", "velocyto")
SUPPORTED_MATRIX_FORMATS = ("tsv", "csv", "mtx")


def _add_common_args(p: argparse.ArgumentParser) -> None:
    grp_io = p.add_argument_group("I/O & orchestration")
    grp_io.add_argument("--output-dir", required=True, help="Output directory.")
    grp_io.add_argument("--output-prefix", default="sanity", help="Output file prefix.")
    grp_io.add_argument("--tmp-dir", default=None, help="Temporary directory (defaults to output/tmp).")
    grp_io.add_argument("--force-overwrite", action="store_true", help="Overwrite existing outputs.")
    grp_io.add_argument("--dry-run", action="store_true", help="Parse and validate, then exit.")
    grp_io.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO", help="Logging level.")

    grp_pre = p.add_argument_group("Preprocessing (Scanpy filters)")
    grp_pre.add_argument("--no-preprocess", action="store_true", help="Skip preprocessing filters.")
    grp_pre.add_argument("--min-counts", type=int, default=500, help="Minimum total UMI counts per cell.")
    grp_pre.add_argument("--min-genes", type=int, default=200, help="Minimum genes detected per cell.")
    grp_pre.add_argument("--min-cells", type=int, default=3, help="Minimum number of cells a gene must appear in.")

    grp_sanity = p.add_argument_group("Sanity execution")
    grp_sanity.add_argument("--sanity-bin", default="Sanity", help="Path to Sanity executable (default assumes it's in $PATH).")
    grp_sanity.add_argument("--systemd", action="store_true", help="Launch via systemd-run.")
    grp_sanity.add_argument("--threads", type=int, default=max(1, mp.cpu_count() // 2))
    grp_sanity.add_argument("--systemd-slice", default=None)
    grp_sanity.add_argument("--systemd-nice", type=int, default=None)
    grp_sanity.add_argument("--vmin", type=float, default=0.001, help="Minimal value of variance in log transcription quotient (default: 0.001).")
    grp_sanity.add_argument("--vmax", type=float, default=50.0, help="Maximal value of variance in log transcription quotient (default: 50).")
    grp_sanity.add_argument("--nbin", type=int, default=160, help="Number of bins for the variance in log transcription quotient (default: 160).")
    grp_sanity.add_argument("--no-norm", action="store_true", help="Skip Sanity's internal cell size normalization (--no_cell_size_normalization).")
    grp_sanity.add_argument("--max-v", action="store_true", help="Use Sanity's maximum-likelihood variance mode (--get_output_for_maxlik_variance).")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="sanity_wrapper",
        description="Wrapper for Sanity single-cell RNA-seq tool — CLI only.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--mode", required=True, choices=SUPPORTED_MODES)

    grp_cellranger = p.add_argument_group("Mode: cellranger")
    grp_cellranger.add_argument("--cellranger-dir", help="Path to 'filtered_feature_bc_matrix' directory.")

    grp_matrix = p.add_argument_group("Mode: matrix")
    grp_matrix.add_argument("--matrix-file", help="Counts matrix file (tsv/csv/mtx).")
    grp_matrix.add_argument("--matrix-format", choices=SUPPORTED_MATRIX_FORMATS, default="tsv")
    grp_matrix.add_argument("--gene-col", default=None)
    grp_matrix.add_argument("--index-col", default=None)
    grp_matrix.add_argument("--transpose", action="store_true")

    grp_velo = p.add_argument_group("Mode: velocyto")
    grp_velo.add_argument("--loom-file", help="Velocyto .loom file.")

    grp_anndata = p.add_argument_group("Mode: anndata")
    grp_anndata.add_argument("--anndata-file", help="Input AnnData (.h5ad) file.")

    _add_common_args(p)
    return p


def _validate_args(args: argparse.Namespace) -> None:
    if args.output_dir:
        args.output_dir = os.path.abspath(os.path.expanduser(args.output_dir))
    if args.tmp_dir is None and args.output_dir:
        args.tmp_dir = os.path.join(args.output_dir, "tmp")

    if args.mode == "cellranger":
        if not args.cellranger_dir:
            raise SystemExit("--cellranger-dir required for cellranger mode")
        if not os.path.exists(args.cellranger_dir):
            raise SystemExit(f"Cell Ranger directory not found: {args.cellranger_dir}")

    elif args.mode == "anndata":
        if not args.anndata_file:
            raise SystemExit("--anndata-file required for anndata mode")
        if not os.path.exists(args.anndata_file):
            raise SystemExit(f"AnnData file not found: {args.anndata_file}")

    elif args.mode == "matrix":
        if not args.matrix_file:
            raise SystemExit("--matrix-file required for matrix mode")
        if not os.path.exists(args.matrix_file):
            raise SystemExit(f"Matrix file not found: {args.matrix_file}")
        if args.index_col and not args.gene_col:
            args.gene_col = args.index_col

    elif args.mode == "velocyto":
        if not args.loom_file:
            raise SystemExit("--loom-file required for velocyto mode")
        if not os.path.exists(args.loom_file):
            raise SystemExit(f"Loom file not found: {args.loom_file}")
    

    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")


def build_parser() -> argparse.ArgumentParser:
    return _build_parser()


def parse_args(argv=None):
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_args(args)
    return args


if __name__ == "__main__":
    import json
    ns = parse_args()
    print(json.dumps(vars(ns), indent=2, sort_keys=True))
