"""
Address reviewer comments and prepare rebuttal.
"""

import os
from pathlib import Path

import altair as alt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from gseapy import parser, prerank
from matplotlib.colors import rgb2hex

from scripts.configs import *
from src import altair_themes
from src.helpers import (
    node_idx,
    node_means,
    save_altair,
    txt_to_list,
    write_adata_label,
)

# Set up Altair theme
alt.themes.register("publishTheme", altair_themes.publishTheme)
alt.themes.enable("publishTheme")
alt.data_transformers.disable_max_rows()


# Define helper functions
def convert_to_cat(df, col, order):
    df[col] = df[col].astype(str).astype("category").cat.reorder_categories(order)


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


def get_log2fc(adata, node_a, node_b, node_info, col_name="deg_groups"):
    """Get log2(B/A) fold change from TooManyCells clusters."""
    adata.obs[col_name] = np.nan
    adata.obs.loc[node_idx(adata.obs, node_a, node_info), col_name] = "A"  # reference
    adata.obs.loc[node_idx(adata.obs, node_b, node_info), col_name] = "B"  # numerator
    sc.tl.rank_genes_groups(
        adata,
        groupby=col_name,
        groups=["B"],
        reference="A",
        layer="lognorm",
        method="wilcoxon",
    )
    log2fc = sc.get.rank_genes_groups_df(adata, group="B")
    log2fc["pvals_adj<=0.05"] = log2fc["pvals_adj"] <= 0.05
    # Sanity check to ensure log2fc is correct
    vals = adata[:, log2fc["names"][0]].layers["lognorm"]
    vals_high = vals[adata.obs[col_name] == "B"].data.mean()
    vals_low = vals[adata.obs[col_name] == "A"].data.mean()
    if vals_high > vals_low:
        return log2fc
    else:
        raise ValueError("Log2FC direction is not correct.")


def get_gsea_preranked(log2fc, geneset_dict, prerank_params):
    results = []
    for geneset_name, geneset in geneset_dict.items():
        res = prerank(rnk=log2fc, gene_sets=geneset, **prerank_params).res2d
        res["Gene set"] = geneset_name
        results.append(res)
    gsea_df = pd.concat(results, axis=0)
    gsea_df["qval<=0.05"] = gsea_df["FDR q-val"] <= 0.05
    return gsea_df


# Load raw read counts and normalize according to standard Scanpy workflow
results_dir = os.path.join(config["results"], "rebuttal")
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(os.path.join(config["scrnaseq_data"], "adata.h5ad"))
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()
sc.pp.scale(adata)

# Perform dimensionality reduction using PCA and UMAP
sc.tl.pca(adata)
sc.pp.neighbors(adata, random_state=1)
sc.tl.umap(adata, random_state=1)
X_umap = adata.obsm["X_umap"].copy()
adata.obs[["UMAP1", "UMAP2"]] = X_umap

# Run Harmony batch correction
adata.obs[["cell", "tx"]] = adata.obs["sample"].str.split("_", expand=True)
sc.external.pp.harmony_integrate(adata, key="cell", basis="X_pca")
sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony")
sc.tl.umap(adata, neighbors_key="harmony")
adata.obs[["Harmony UMAP1", "Harmony UMAP2"]] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap"] = X_umap.copy()

# Cluster using uncorrected neighbours and plot UMAP
res = 0.1
sc.tl.leiden(adata, neighbors_key="neighbors", resolution=res)
clusters = adata.obs["leiden"].sort_values().unique().astype(str)
convert_to_cat(adata.obs, "leiden", clusters)
set1 = [rgb2hex(c) for c in plt.get_cmap("Set1").colors]
leiden_cmap = {
    "domain": clusters,
    "range": set1[: len(clusters)],
}
adata.uns["leiden_colors"] = leiden_cmap["domain"]
umap = plot_umap(adata.obs, "leiden", "Leiden", leiden_cmap)
save_altair(umap, f"umap_leiden_{res}", results_dir)

# Cluster using Harmony-corrected neighbours and plot UMAP
lh = "leiden_harmony"
sc.tl.leiden(adata, neighbors_key="harmony", key_added=lh, resolution=res)
clusters_harmony = adata.obs[lh].sort_values().unique().astype(str)
convert_to_cat(adata.obs, lh, clusters_harmony)
leiden_harmony_cmap = {
    "domain": clusters_harmony,
    "range": set1[: len(clusters_harmony)],
}
adata.uns[f"{lh}_colors"] = leiden_harmony_cmap["domain"]
umap_harmony = plot_umap(adata.obs, lh, "Leiden", leiden_harmony_cmap, "Harmony UMAP")
save_altair(umap_harmony, f"umap_{lh}_{res}", results_dir)

# Calculate cell cycle phase on lognorm transformed counts
adata.X = adata.layers["lognorm"].copy()
s_phase = txt_to_list(config["s_phase"])
g2m_phase = txt_to_list(config["g2m_phase"])
sc.tl.score_genes_cell_cycle(adata, s_genes=s_phase, g2m_genes=g2m_phase)
phases = ["G1", "S", "G2M"]
convert_to_cat(adata.obs, "phase", phases)
phase_cmap = {
    "domain": phases,
    "range": set1[: len(phases)],
}
adata.uns["phase_colors"] = phase_cmap["range"]

# Format observation names
obs = adata.obs.copy()
obs["tx"] = obs["tx"].fillna(value="treated")
obs["cell"] = obs["cell"].replace(
    {
        "dnd41": "DND-41",
        "pc9": "PC9",
        "lncap": "LNCaP",
        "mdamb231": "MDA-MB-231",
        "skmel28": "SK-MEL-28",
    }
)
obs.loc[(obs["sample"] == "dnd41_long"), "tx"] = "treated (long)"
obs.loc[(obs["sample"] == "dnd41_short"), "tx"] = "treated (short)"
obs["name"] = obs[["cell", "tx"]].astype(str).agg(" ".join, axis=1)

# Overlay TooManyCells cluster information and calculate diapause means per node
clusters_path = os.path.join(config["pruned"], "clusters.csv")
clusters = pd.read_csv(clusters_path, index_col="cell")
obs["TMC"] = clusters["cluster"].astype(int).astype(str)
node_info = pd.read_csv(os.path.join(config["pruned"], "node_info.csv"))
diapause_df = node_means(obs, "diapause_score", node_info)
diapause_path = os.path.join(results_dir, "annotate_diapause.csv")
diapause_df.to_csv(diapause_path, index=False, header=True)

# Bin diapause score
obs["diapause_binned"] = pd.cut(obs["diapause_score"], bins=12)
bins = obs["diapause_binned"].cat.categories.astype(str)
convert_to_cat(obs, "diapause_binned", bins)
cmap = plt.get_cmap("YlGnBu", len(bins) - 1)
hexmap = [rgb2hex(cmap(i)) for i in range(cmap.N)][::-1] + ["#ffffff"]
adata.uns["diapause_binned_colors"] = hexmap

# Add colour scheme for sample conditions
paired = [rgb2hex(c) for c in plt.get_cmap("Paired").colors][:10]
paired.insert(2, "#003f94")
name_cmap = {
    "domain": [
        "DND-41 control",
        "DND-41 treated (short)",
        "DND-41 treated (long)",
        "LNCaP control",
        "LNCaP treated",
        "MDA-MB-231 control",
        "MDA-MB-231 treated",
        "PC9 control",
        "PC9 treated",
        "SK-MEL-28 control",
        "SK-MEL-28 treated",
    ],
    "range": paired,
}
convert_to_cat(obs, "name", name_cmap["domain"])
adata.uns["name_colors"] = name_cmap["range"]

# Export PCAs for use with TooManyCells
pca = pd.DataFrame(
    adata.obsm["X_pca"],
    index=adata.obs_names,
    columns=[f"PC{i}" for i in range(1, 51)],
)
pca.T.to_csv(os.path.join(results_dir, "pca/matrix.csv"))

pca_harmony = pd.DataFrame(
    adata.obsm["X_pca_harmony"],
    index=adata.obs_names,
    columns=[f"PC{i}" for i in range(1, 51)],
)
pca_harmony.T.to_csv(os.path.join(results_dir, "pca_harmony/matrix.csv"))

# Load gene sets and analysis parameters for GSEA
geneset_dict = {}
geneset_names = {"h": "Hallmark", "c2": "C2", "c6": "C6"}
for geneset_file in os.listdir(config["gmt_msigdb"]):
    geneset_path = os.path.join(config["gmt_msigdb"], geneset_file)
    geneset_name = geneset_names[geneset_file.split(".")[0]]
    geneset_dict[geneset_name] = parser.read_gmt(geneset_path)
prerank_params = {
    "n_perms": 100,
    "min_size": 3,
    "max_size": 1000,
    "ascending": False,
    "threads": 100,
    "seed": 1,
    "outdir": None,
}

# Run differential expression analysis for node 48 (ID2 high) vs. 45 (ID2 low)
adata.obs = obs.copy()
deg_id2 = get_log2fc(adata, 45, 48, node_info)
deg_id2.to_csv(os.path.join(results_dir, "log2fc_id2_48v45.csv"), index=False)
cols = ["names", "logfoldchanges"]
gsea_id2 = get_gsea_preranked(deg_id2[cols], geneset_dict, prerank_params)
gsea_id2.to_csv(os.path.join(results_dir, "gsea_id2_48v45.csv"), index=False)

# Run DEG on node 36 vs. 128
deg_breast = get_log2fc(adata, 128, 36, node_info)
deg_breast.to_csv(os.path.join(results_dir, "log2fc_breast_36v128.csv"), index=False)
deg_breast = deg_breast[deg_breast["pvals_adj"] <= 0.05]
gsea_breast = get_gsea_preranked(deg_breast[cols], geneset_dict, prerank_params)
gsea_breast.to_csv(os.path.join(results_dir, "gsea_breast_128v36.csv"), index=False)

# Run DEG on Harmony cluster 7 vs rest
sc.tl.rank_genes_groups(
    adata,
    groupby=lh,
    groups=["7"],
    reference="rest",
    layer="lognorm",
    method="wilcoxon",
)
deg_leiden = sc.get.rank_genes_groups_df(adata, group="7")
deg_leiden["pvals_adj<=0.05"] = deg_leiden["pvals_adj"] <= 0.05
deg_leiden.to_csv(os.path.join(results_dir, f"log2fc_{lh}_7vrest.csv"), index=False)
gsea_leiden = get_gsea_preranked(deg_leiden[cols], geneset_dict, prerank_params)
gsea_leiden.to_csv(os.path.join(results_dir, f"gsea_{lh}_7vrest.csv"), index=False)

# Save results to file
adata.write_h5ad(os.path.join(results_dir, "adata.h5ad"))
for key in ["phase", "name", "leiden", "leiden_harmony", "diapause_binned"]:
    write_adata_label(adata, key, results_dir)
