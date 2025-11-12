__author__ = "Maulik Nariya" 
__date__ = "November 2025"
__copyright__ = "MIT license"

"""
Sanity single-cell RNA-seq wrapper — Main runner

This script coordinates the full workflow:
  1. Parse CLI arguments (via cli.py)
  2. Load input data (cellranger, matrix, or velocyto)
  3. Apply Scanpy preprocessing filters
  4. Export gene×cell TSV files for Sanity
  5. Run Sanity using systemd-run or directly
  6. Load LTQs and LTQ errors into AnnData layers
  7. Save final .h5ad file
"""

import os
import subprocess
import scanpy as sc
from anndata import AnnData
import numpy as np
import pandas as pd

from cli import parse_args
import preprocessing as prep


def run_sanity_tsv(tsv_file: str, args, label: str = "matrix"):
    out_dir = os.path.dirname(tsv_file)
    base = os.path.splitext(os.path.basename(tsv_file))[0]
    out_ltq = os.path.join(out_dir, "log_transcription_quotients.txt")
    out_err = os.path.join(out_dir, "ltq_error_bars.txt")

    cmd = [
        args.sanity_bin,
        "-f", tsv_file,
        "-d", out_dir,
        "-n", str(args.threads),
        "-vmin", str(args.vmin),
        "-vmax", str(args.vmax),
        "-nbin", str(args.nbin),
        "-max_v", str(args.max_v),
        "-no_norm", str(args.no_norm)
        
    ]
    # Optional Sanity boolean flags
    if args.no_norm:
        cmd += ["-no_norm", "true"]
    if args.max_v:
        cmd += ["-max_v", "true"]

    if args.systemd:
        systemd_cmd = ["systemd-run", "--scope", "--collect", "--quiet"]
        if args.systemd_slice:
            systemd_cmd += ["--slice", args.systemd_slice]
        if args.systemd_nice is not None:
            systemd_cmd += ["--nice", str(args.systemd_nice)]
        systemd_cmd += cmd
        cmd = systemd_cmd

    print(f"[INFO] Running Sanity ({label})...")
    print(" ".join(cmd))
    
    subprocess.run(cmd, check=True)
    # Rename Sanity outputs based on mode and label
    prefix = args.output_prefix
    out_dir = args.output_dir
    
    # Default Sanity output filenames
    ltq_src = os.path.join(out_dir, "log_transcription_quotients.txt")
    err_src = os.path.join(out_dir, "ltq_error_bars.txt")
    
    # Determine destination filenames
    if args.mode == "velocyto":
        # label can be matrix / spliced / unspliced
        if label == "matrix":
            ltq_dst = os.path.join(out_dir, f"{prefix}_matrix_LTQ.txt")
            err_dst = os.path.join(out_dir, f"{prefix}_matrix_LTQerr.txt")
        elif label == "spliced":
            ltq_dst = os.path.join(out_dir, f"{prefix}_spliced_LTQ.txt")
            err_dst = os.path.join(out_dir, f"{prefix}_spliced_LTQerr.txt")
        elif label == "unspliced":
            ltq_dst = os.path.join(out_dir, f"{prefix}_unspliced_LTQ.txt")
            err_dst = os.path.join(out_dir, f"{prefix}_unspliced_LTQerr.txt")
        else:
            ltq_dst = os.path.join(out_dir, f"{prefix}_{label}_LTQ.txt")
            err_dst = os.path.join(out_dir, f"{prefix}_{label}_LTQerr.txt")
    else:
        # non-velocyto modes
        ltq_dst = os.path.join(out_dir, f"{prefix}_LTQ.txt")
        err_dst = os.path.join(out_dir, f"{prefix}_LTQerr.txt")
    
    # Rename Sanity outputs
    os.rename(ltq_src, ltq_dst)
    os.rename(err_src, err_dst)

    try:
        os.remove(tsv_file)
        print(f"[INFO] Removed intermediate file: {tsv_file}")
    except Exception as e:
        print(f"[WARN] Could not remove intermediate file {tsv_file}: {e}")
    
    return ltq_dst, err_dst


def read_ltq_file(path: str) -> np.ndarray:
    """Read Sanity LTQ or LTQerr TSV file into numpy array."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    return df.values


def attach_sanity_outputs(adata: AnnData, ltq_path: str, err_path: str, label: str):
    """Attach LTQ and LTQerr outputs to AnnData.layers."""
    ltq_df = pd.read_csv(ltq_path, sep="\t", index_col=0)
    ltq_err_df = pd.read_csv(err_path, sep="\t", index_col=0)
    
    # Align to master gene list (adata.var_names)
    ltq_df = ltq_df.reindex(adata.var_names).fillna(0)
    ltq_err_df = ltq_err_df.reindex(adata.var_names).fillna(0)
    
    # Transpose to cells × genes
    ltq = ltq_df.values.T
    ltq_err = ltq_err_df.values.T
    
    adata.layers[f"{label}_LTQ"] = ltq
    adata.layers[f"{label}_LTQerr"] = ltq_err


def main():
    args = parse_args()

    # Step 1: Load data
    print('[INFO] Loading input data (this may take a few minutes depending on file size)...')
    if args.mode == "cellranger":
        adata = prep.load_cellranger(args.cellranger_dir)
    elif args.mode == "anndata":
        adata = prep.load_anndata(args.anndata_file)
    elif args.mode == "matrix":
        adata = prep.load_matrix_tabular(args.matrix_file, args.matrix_format, args.gene_col, args.transpose)
    elif args.mode == "velocyto":
        adata = prep.load_velocyto_loom(args.loom_file)
    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

    # Step 2: Preprocessing
    print("[INFO] Preprocessing: filtering cells and genes with low counts...")
    adata = prep.apply_filters(adata, args)
    
    # Store raw matrix layer (only for cellranger and simple matrix modes)
    if args.mode in {"cellranger", "matrix", "anndata"}:
        adata.layers["matrix"] = adata.X.copy()

    # Step 3: Export TSV(s)
    tsv_files = prep.export_tsvs(adata, args)

    # Step 4: Run Sanity
    sanity_outputs = {}
    for label, tsv_path in tsv_files.items():
        ltq_path, ltq_err = run_sanity_tsv(tsv_path, args, label)
        sanity_outputs[label] = (ltq_path, ltq_err)

    # Step 5: Integrate results
    for label, (ltq_path, ltq_err) in sanity_outputs.items():
        attach_sanity_outputs(adata, ltq_path, ltq_err, label)

    # Step 6: Save AnnData
    out_h5ad = os.path.join(args.output_dir, f"{args.output_prefix}.h5ad")
    print(f"[INFO] Saving AnnData to {out_h5ad}")
    adata.write_h5ad(out_h5ad)
    print(f"[INFO] Execution complete. Output is saved as '{out_h5ad}'.")

    # Cleanup: remove Sanity output text files
    to_delete = []
    if args.mode == "velocyto":
        to_delete += [
            os.path.join(args.output_dir, f"{args.output_prefix}_matrix_LTQ.txt"),
            os.path.join(args.output_dir, f"{args.output_prefix}_matrix_LTQerr.txt"),
            os.path.join(args.output_dir, f"{args.output_prefix}_spliced_LTQ.txt"),
            os.path.join(args.output_dir, f"{args.output_prefix}_spliced_LTQerr.txt"),
            os.path.join(args.output_dir, f"{args.output_prefix}_unspliced_LTQ.txt"),
            os.path.join(args.output_dir, f"{args.output_prefix}_unspliced_LTQerr.txt"),
        ]
    else:
        to_delete += [
            os.path.join(args.output_dir, f"{args.output_prefix}_LTQ.txt"),
            os.path.join(args.output_dir, f"{args.output_prefix}_LTQerr.txt"),
        ]
    
    for f in to_delete:
        if os.path.exists(f):
            try:
                os.remove(f)
                print(f"[INFO] Removed Sanity output: {f}")
            except Exception as e:
                print(f"[WARN] Could not remove {f}: {e}")


if __name__ == "__main__":
    main()
