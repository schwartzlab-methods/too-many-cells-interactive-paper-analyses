"""
Utilites to prepare read count matrix for analysis
"""


import scanpy as sc
from scipy.sparse import issparse


### Checks for read count data
def is_counts(mtx):
    if issparse(mtx):
        x = mtx.data
    else:
        x = mtx
    return ((x - x.astype(int)) == 0).all()


def normalize(adata):
    # Check that raw counts are present in the input matrix
    if not is_counts(adata.layers["raw"]):
        raise ValueError("Input matrix does not contain raw counts")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['lognorm'] = adata.X.copy()


def reduce_dimensions(adata, leiden_res=0.1, id=""):
    """
    Perform dimensionality reduction and clustering.
    """
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=leiden_res, key_added=f"leiden{id}")
    adata.obs[[f"UMAP1{id}", f"UMAP2{id}"]] = adata.obsm["X_umap"][:, :2]


def get_qc_metrics(adata, mito_prefix="MT-", ribo_prefix=["RPL", "RPS"]):
    """
    Annotate AnnData object with quality control metrics for single-cell
    sequencing data.
    """
    if not is_counts(adata.layers["raw"]):
        raise ValueError("Input matrix does not contain raw counts")
    # Calculate expression of mitochodnrial and ribosomal genes
    adata.obs["pct_mito"] = score_expression(adata, gene_prefix=mito_prefix) * 100
    adata.obs["pct_ribo"] = score_expression(adata, gene_prefix=ribo_prefix) * 100
    # Calculate total UMI counts and detected genes per cell
    adata.obs["total_umi"] = adata.layers["raw"].sum(axis=1).astype(int)
    adata.obs["n_genes"] = (adata.layers["raw"] > 0).sum(axis=1)
    # Calculate total UMI counts and detected cells per gene
    adata.var["total_umi"] = adata.layers["raw"].sum(axis=0).A1.astype(int)
    adata.var["n_cells"] = (adata.layers["raw"] > 0).sum(axis=0).A1


### Scoring functions
def load_genes(genes_path):
    with open(genes_path, "r", encoding="utf-8") as f:
        gene_list = f.read().splitlines()
    return gene_list


def score_cell_cycle(adata, s_genes, g2m_genes, layer="lognorm"):
    s = s_genes
    g2m = g2m_genes
    adata.X = adata.layers["lognorm"].copy()
    sc.tl.score_genes_cell_cycle(adata, s_genes=s, g2m_genes=g2m)


def score_expression(
    adata, gene_prefix=None, gene_list=None, layer="raw", average=False, ratio=True
):
    # Handle gene prefix options
    if not gene_prefix:
        gene_prefix = ""
    elif isinstance(gene_prefix, list):
        gene_prefix = tuple(gene_prefix)
    # Check if genes are present in the AnnData object and subset to matches
    prefix_vars = adata.var_names.str.startswith(gene_prefix)
    found_vars = True if gene_list is None else adata.var_names.isin(gene_list)
    gene_idx = prefix_vars & found_vars
    if gene_idx.sum() == 0:
        raise ValueError("No genes found with the specified prefix(es)")
    # Calculate gene expression scores
    x = adata.layers[layer]
    if average:
        return x[:, gene_idx].mean(axis=1)
    else:
        geneset_counts = x[:, gene_idx].sum(axis=1)
    if ratio:
        return geneset_counts / x.sum(axis=1)
    else:
        return geneset_counts
