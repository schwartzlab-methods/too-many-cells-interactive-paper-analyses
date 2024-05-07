"""
Helper functions for TooManyCells Interactive demonstration
"""

import gzip
import os
from pathlib import Path
from typing import Dict

import altair as alt
import anndata as ad
import gseapy as gp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, rgb2hex, to_hex
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from scipy.stats import kruskal, mannwhitneyu
from scipy.stats.mstats import gmean
from statsmodels.stats.multitest import multipletests

from . import altair_themes

# Set up Altair theme
alt.themes.register("publishTheme", altair_themes.publishTheme)
alt.themes.enable("publishTheme")
alt.data_transformers.disable_max_rows()


def find_parent_dir(dirpath, filename):
    for root, _, files in os.walk(dirpath, topdown=False):
        if any(filename in f for f in files):
            return root
    return None


# Operations on sparse matrices
def _apply_by_row(func, input_mtx, *args, **kwargs):
    ncells = input_mtx.shape[0]
    output_mtx = input_mtx.copy()
    for cell in range(ncells):
        cellreads = input_mtx.getrow(cell)
        new_cellreads = func(cellreads.data, *args, **kwargs)
        output_mtx.data[output_mtx.indptr[cell] : output_mtx.indptr[cell + 1]] = (
            new_cellreads
        )
    return output_mtx


# Normalize data
def standardize(mtx):
    def row_standardize(reads):
        std = np.std(reads)
        norm = (reads - np.mean(reads)) / std
        return norm

    matrixnorm = _apply_by_row(row_standardize, mtx)
    return matrixnorm


def normalize(adata, target_sum=1e4):
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()


# Dimensionality reduction
def reduce_dimensions(adata):
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, random_state=1)
    sc.tl.umap(adata, random_state=1)
    X_umap = adata.obsm["X_umap"].copy()
    adata.obs[["UMAP1", "UMAP2"]] = X_umap


# Gene signature scoring functions
def signature_score(adata, geneset_path, geneset_id, inplace=True):
    """
    Calculate gene signature scores from average expression of weighted genes.
    """
    try:
        geneset = pd.read_csv(geneset_path)
    except AttributeError:
        "Cannot open path to geneset.csv file"
    # Scaling allows for read counts to be centered around the mean, so that
    # final scores can be maximized regardless of expression directionality
    try:
        raw = adata.layers["raw"].todense()
    except AttributeError:
        "Cannot find raw layer in AnnData object"
    stdev = np.std(raw, axis=0)
    znorm = raw - np.mean(raw, axis=0) / stdev
    reads = (
        pd.DataFrame(znorm, index=adata.obs_names, columns=adata.var_names)
        .fillna(0)
        .replace([np.inf, -np.inf], 0)
    )
    # Weights correspond to upregulated/downregulated genes
    gene_up = [gene for gene in geneset[geneset["weight"] >= 0]["gene"]]
    gene_down = [gene for gene in geneset[geneset["weight"] < 0]["gene"]]
    adata_up = reads[reads.columns.intersection(gene_up)]
    adata_down = reads[reads.columns.intersection(gene_down)] * -1
    adata_all = pd.concat([adata_up, adata_down], axis=1)
    gene_signature_score = np.round(adata_all.mean(axis=1), 3)
    if inplace == True:
        adata.obs.loc[:, f"{geneset_id}_score"] = gene_signature_score
    else:
        return gene_signature_score


def txt_to_list(txt_path):
    with open(txt_path, "r") as f:
        return f.read().splitlines()


# GSEA functions
def get_log2fc(adata, node_info, node_a=None, node_b=None, col_name="deg_groups"):
    """
    Get log2(B/A) fold change from TooManyCells clusters.
    """
    adata.obs[col_name] = np.nan
    if not node_a and not node_b:
        raise ValueError("Please provide node IDs for comparison.")
    if node_a:
        a = "A"
        adata.obs.loc[node_idx(adata.obs, node_a, node_info), col_name] = a
    else:
        a = "rest"
    if node_b:
        b = "B"
        adata.obs.loc[node_idx(adata.obs, node_b, node_info), col_name] = b
    else:
        b = "all"
    sc.tl.rank_genes_groups(
        adata,
        groupby=col_name,
        groups=[b],
        reference=a,
        layer="lognorm",
        method="wilcoxon",
    )
    log2fc = sc.get.rank_genes_groups_df(adata, group="B")
    log2fc["pvals_adj<=0.05"] = log2fc["pvals_adj"] <= 0.05
    log2fc.sort_values("logfoldchanges", ascending=False, inplace=True)
    return log2fc


def get_gsea_preranked(log2fc, geneset_dict, prerank_params):
    results = []
    for geneset_name, geneset in geneset_dict.items():
        res = gp.prerank(rnk=log2fc, gene_sets=geneset, **prerank_params).res2d
        res["Gene set"] = geneset_name
        results.append(res)
    gsea_df = pd.concat(results, axis=0)
    gsea_df.sort_values("NES", ascending=False, inplace=True)
    gsea_df["qval<=0.05"] = gsea_df["FDR q-val"] <= 0.05
    return gsea_df


def preranked(log2fc, gene_dict: Dict, min_size=5):
    pre_res = gp.prerank(
        rnk=log2fc,
        gene_sets=gene_dict,
        threads=4,
        min_size=min_size,
        max_size=1000,
        permutation_num=100,
        outdir=None,
        seed=1,
    )
    return pre_res


def onevsrest(
    adata: ad.AnnData, layer: str, cluster_path: os.PathLike, save_dir: os.PathLike
):
    def colwise_mean(mtx: csr_matrix) -> np.ndarray:
        sums = mtx.sum(axis=0).A1
        nonzero = mtx.getnnz(axis=0)
        mean = sums / nonzero
        return mean

    # Create dictionary of node IDs and corresponding barcodes
    clusters_df = pd.read_csv(cluster_path, index_col="cell")
    clusters_dict = clusters_df["path"].to_dict()
    node_dict = {}
    for barcode, path in clusters_dict.items():
        path_list = [int(x) for x in path.split("/")]
        for node in path_list:
            if node not in node_dict:
                node_dict[node] = [barcode]
            else:
                node_dict[node].append(barcode)
    # Calculate log2FC for each node
    df_list = []
    for node_id, barcodes in node_dict.items():
        a = adata[barcodes].layers[layer]
        b = adata[~adata.obs_names.isin(barcodes), :].layers[layer]
        ratio = colwise_mean(a) / colwise_mean(b)
        vars = adata.var_names[~np.isnan(ratio)]
        log2fc = np.log2(ratio[~np.isnan(ratio)])
        log2fc_df = pd.DataFrame(log2fc, index=vars, columns=["log2fc"]).dropna()
        log2fc_df.sort_values(by="log2fc", ascending=False, inplace=True)
        log2fc_path = os.path.join(save_dir, f"log2fc_{node_id}_vs_rest.csv")
        log2fc_df.to_csv(log2fc_path)
        log2fc_df["node_id"] = node_id
        df_list.append(log2fc_df)
    return pd.concat(df_list)


def mannwhitney(obs: pd.DataFrame, obs_key: str) -> pd.DataFrame:
    score_dict = {}
    score_list = []
    for sample, scores in obs.groupby("sample"):
        score_dict[sample] = scores[obs_key].to_list()
        score_list.append(scores[obs_key].to_list())
    cell_line = []
    u_score = []
    pvals = []
    pvals_adj = []
    for cell in ["lncap", "skmel28", "pc9", "mdamb231", "dnd41_long", "dnd41_short"]:
        if cell.startswith("dnd41"):
            ctrl = "dnd41_control"
        else:
            ctrl = f"{cell}_control"
        res = mannwhitneyu(x=score_dict[ctrl], y=score_dict[cell])
        cell_line.append(cell)
        u_score.append(res[0])
        pvals.append(res[1])
    pvals_adj = multipletests(pvals, alpha=0.05, method="fdr_bh")[1]
    df = pd.DataFrame(
        list(zip(cell_line, u_score, pvals, pvals_adj)),
        columns=["Cell line", "U score", "p-value", "p-value (adj)"],
    )
    h = kruskal(*score_list)
    df["Kruskal-Wallis H"] = h[0]
    df["Kruskal-Wallis p-value"] = h[1]
    return df


def rank_product(foldchange: pd.DataFrame, n_permutations: int = 100):
    """
    Calculate the rank product of the foldchange values.
    """
    rank = foldchange.rank(axis=0, ascending=False)
    rank["geo_mean"] = gmean(rank, axis=1, nan_policy="omit")
    rank.sort_values("geo_mean", inplace=True)
    rank["rank"] = range(1, len(rank) + 1)
    # Calculate the significance for each rank product value by permuting the
    # foldchange values randomly, calculating the rank product, and comparing
    # how likely the given value or better is observed.
    shuffle = rank.drop(["geo_mean", "rank"], axis=1).copy()
    c = pd.DataFrame(index=shuffle.index)
    for k in range(n_permutations):
        for col in shuffle.columns:
            shuffle[col] = shuffle[col].sample(frac=1).values
        geo_mean = gmean(shuffle, axis=1)
        c[k] = geo_mean <= rank["geo_mean"]
    erp = c.sum(axis=1) / n_permutations
    pfp = erp / rank["rank"]
    rank["pVal"] = pfp
    return rank


# Saving MEX files and corresponding label metadata
def adata_to_10x(adata, key, mex_path, attr="layers"):
    Path(mex_path).mkdir(parents=True, exist_ok=True)
    # Transpose read counts to matrix market format, such that (rows, columns)
    # correspond to (features, cells)
    with gzip.open(os.path.join(mex_path, "matrix.mtx.gz"), "wb") as mtx:
        mmwrite(mtx, getattr(adata, attr)[key].T, symmetry="symmetric")
    # Barcode observations
    barcodes_path = os.path.join(mex_path, "barcodes.tsv.gz")
    adata.obs.to_csv(barcodes_path, columns=[], sep="\t", header=False)
    # Gene features
    features_path = os.path.join(mex_path, "features.tsv.gz")
    var_df = pd.DataFrame(
        {
            "gene_id": adata.var.reset_index()["index"],
            "gene_symbol": adata.var.reset_index()["index"],
            "type": "Gene Expression",
        }
    )
    var_df.to_csv(features_path, sep="\t", header=False, index=False)


def export_labels(folder, adata_obs, obs_key):
    labels = adata_obs[obs_key].reset_index()
    labels.columns = ["item", "label"]
    labels.to_csv(os.path.join(folder, f"labels_{obs_key}.csv"), index=None)


def write_adata_label(adata, obs_key, data_dir):
    # Write labels file for TooManyCells
    labels_target = os.path.join(data_dir, f"labels_{obs_key}.csv")
    adata_obs = pd.DataFrame(adata.obs.reset_index())
    labels_df = adata_obs[["index", obs_key]]
    labels_df.columns = ["item", "label"]
    labels_df.to_csv(labels_target, header=True, index=False)


def convert_to_cat(df, col, order):
    df[col] = df[col].astype(str).astype("category").cat.reorder_categories(order)


# Handling TooManyCells hierarchical node information
def node_means(obs, obs_key, node_info, cluster_key="tmc"):
    subtree_dict = dict(zip(node_info["node"], node_info["subtree"].str.split("/")))
    vals = {}
    for node_id, subtree in subtree_dict.items():
        vals[node_id] = obs[obs_key][obs[cluster_key].isin(subtree)].mean()
    vals_df = pd.DataFrame().from_dict(vals, orient="index")
    vals_df = vals_df.reset_index().rename(columns={"index": "node_id", 0: obs_key})
    return vals_df


def node_idx(obs, node, node_info, cluster_key="tmc"):
    subtree_dict = dict(zip(node_info["node"], node_info["subtree"].str.split("/")))
    idx = obs[cluster_key].isin(subtree_dict[node])
    return idx


# Plotting functions
def save_altair(plot, plot_id, results_dir, img_formats=["svg", "html", "png"]):
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    for img_format in img_formats:
        plot_path = f"{results_dir}/{plot_id}.{img_format}"
        plot.save(plot_path)


def plot_umap(obs, c, ctitle, cmap, embedding="UMAP"):
    umap = (
        alt.Chart(obs)
        .mark_circle(size=3, opacity=0.8)
        .encode(
            alt.X(f"{embedding}1:Q"),
            alt.Y(f"{embedding}2:Q"),
            alt.Color(c).title(ctitle).scale(**cmap),
        )
    )
    return umap


def interpolate_colors(x, palette="tab20", hex_list=None):
    if not isinstance(x, int):
        x = len(x)
    if hex_list is None:
        # Palette can either be a continuous or continuous color map
        cmap = plt.get_cmap(palette)
        if isinstance(cmap, ListedColormap):
            hex_list = [to_hex(c) for c in cmap.colors]
        else:
            cmap = plt.get_cmap(palette, x)
            hex_list = [rgb2hex(cmap(i)) for i in range(cmap.N)]
    if x <= len(hex_list):
        colors = hex_list[:x]
        return colors
    else:
        cmap = LinearSegmentedColormap.from_list("cmap", hex_list)
        extended_scheme = cmap(np.linspace(0, 1, x))
        colors = [to_hex(c, keep_alpha=False) for c in extended_scheme]
        return colors
