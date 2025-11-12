# pySanity
pySanity is a python wrapper that integrates Sanity into the scanpy ecosystem. In addition to input formats accepted by Sanity natively, the cellranger output `mtx` and the `tsv` files, it also accepts Anndata `.h5ad` file and the velocyto `loom` file. This can be specified using the `--mode` option. **Note** that in the `--matrix` mode it expects a cell by gene, as opposed to gene by cell matrix as expected by Sanity and in `--anndata` mode it expects to have raw counts in `adata.X`. By default, before running Sanity it peforms preprocessing to remove cells and genes with low counts, this can be changed using the `--no-preprocess` option. It saves the output in an Anndata file (in `adata.layers['matrix_LTQ']` and `adata.layers['matrix'_LTQerr]`).

### Dependencies
Python ≥ 3.9  
Sanity  
scanpy  
