__author__ = "Maulik Nariya" 
__date__ = "November 2025"
__copyright__ = "MIT license"

"""
Sanity single-cell RNA-seq wrapper — Preprocessing

This module provides functions to:
  • Load input data into AnnData (Cell Ranger, TSV/CSV/TXT, Velocyto)
  • Accept both compressed (.gz) and uncompressed files for tabular input
  • Apply Scanpy-based filters (min_counts, min_genes, min_cells)
  • Export preprocessed data as uncompressed gene×cell TSV files for Sanity.
"""

import gzip
import os
import sys
from typing import Optional, Iterable, Mapping

import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import scanpy as sc
    from anndata import AnnData
except ImportError:
    sys.stderr.write("[FATAL] scanpy/anndata are required. Install with: pip install scanpy anndata\n")
    raise


def ensure_dirs(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def open_maybe_gzip(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def write_gene_by_cell_tsv(
    X: sp.spmatrix | np.ndarray,
    gene_names: Iterable[str],
    cell_names: Iterable[str],
    out_path: str,
    chunk_rows: int = 2000,
) -> None:
    gene_names = np.asarray(list(gene_names))
    cell_names = np.asarray(list(cell_names))

    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)

    with open(out_path, "wt") as fh:
        fh.write("gene\t" + "\t".join(map(str, cell_names)) + "\n")
        n_genes = X.shape[0]
        for start in range(0, n_genes, chunk_rows):
            end = min(start + chunk_rows, n_genes)
            block = X[start:end].toarray() if sp.issparse(X) else X[start:end]
            for i, row in enumerate(block):
                row_str = "\t".join(map(str, row.astype(int)))
                fh.write(f"{gene_names[start + i]}\t{row_str}\n")


# ====================
# Loaders per input mode
# ====================

def load_cellranger(cellranger_dir: str) -> AnnData:
    print(f"[INFO] Loading Cell Ranger data from {cellranger_dir} ...")
    adata = sc.read_10x_mtx(cellranger_dir, var_names="gene_ids", make_unique=True)
    adata.var_names_make_unique()
    return adata

def load_anndata(path: str) -> AnnData:
    print(f"[INFO] Loading AnnData file {path} ...")
    adata = sc.read_h5ad(path)
    adata.var_names_make_unique()
    return adata


def load_matrix_tabular(path: str, fmt: Optional[str], gene_col: Optional[str], transpose: bool) -> AnnData:
    print(f"[INFO] Loading matrix file {path} ...")
    fmt = fmt or os.path.splitext(path)[1].lstrip('.')
    if fmt in {"tsv", "tsv.gz", "txt", "txt.gz"}:
        sep = "\t"
    elif fmt in {"csv", "csv.gz"}:
        sep = ","
    else:
        raise ValueError(f"Unsupported file format: {fmt}")

    open_func = gzip.open if path.endswith(".gz") else open
    with open_func(path, 'rt') as f:
        df = pd.read_csv(f, sep=sep, header=0)

    if gene_col is None:
        gene_col = df.columns[0]

    if gene_col in df.columns:
        df = df.set_index(gene_col)

    if transpose:
        df = df.T

    X = sp.csr_matrix(df.values)
    adata = AnnData(X=X)
    adata.var_names = df.index.astype(str)
    adata.obs_names = df.columns.astype(str)
    return adata


def load_velocyto_loom(path: str) -> AnnData:
    print(f"[INFO] Loading Velocyto loom file {path} ...")
    adata = sc.read_loom(path, sparse=True, validate=False)
    adata.var_names_make_unique()
    return adata


# ==========================
# Preprocessing (Scanpy pp)
# ==========================

def apply_filters(adata: AnnData, args) -> AnnData:
    if args.no_preprocess:
        return adata

    print(f"[INFO] Initial shape before filtering: {adata.shape}")
    
    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)

    # Apply filtering to all layers (important for velocyto)
    obs_mask = adata.obs_names.isin(adata.obs_names)
    var_mask = adata.var_names.isin(adata.var_names)
    obs_idx = np.where(obs_mask)[0]
    var_idx = np.where(var_mask)[0]

    for layer in list(adata.layers.keys()):
        layer_data = adata.layers[layer]
        if sp.issparse(layer_data):
            adata.layers[layer] = layer_data[obs_idx[:, None], var_idx]
        else:
            adata.layers[layer] = layer_data[obs_idx[:, None], var_idx]

    print(f"[INFO] Shape after filtering: {adata.shape}")
    return adata
    

# ==============================
# Export TSVs for each run mode
# ==============================


def export_tsvs(adata: AnnData, args) -> Mapping[str, str]:
    ensure_dirs(args.output_dir)
    prefix = os.path.join(args.output_dir, args.output_prefix)
    written = {}

    # Always write the main counts matrix
    out_counts = prefix + ".matrix.tsv"
    write_gene_by_cell_tsv(adata.X.T, adata.var_names, adata.obs_names, out_counts)
    written["matrix"] = out_counts

    # For velocyto mode, also write spliced/unspliced layers if present
    if args.mode == "velocyto":
        for layer in ("spliced", "unspliced"):
            if layer in adata.layers:
                outp = f"{prefix}.{layer}.tsv"
                write_gene_by_cell_tsv(adata.layers[layer].T, adata.var_names, adata.obs_names, outp)
                written[layer] = outp

    return written
