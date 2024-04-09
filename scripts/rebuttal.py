"""
Address reviewer comments and prepare rebuttal.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from matplotlib.colors import rgb2hex

from scripts.configs import *
from src.helpers import node_means, txt_to_list, write_adata_label

# Load raw read counts and normalize according to standard Scanpy workflow
adata = sc.read_h5ad(f"{config['scrnaseq_data']}/adata.h5ad")
adata.obs["total_umi"] = adata.layers["raw"].sum(axis=1).astype(int)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
adata.layers["lognorm"] = adata.X.copy()
sc.pp.scale(adata)

# Perform dimensionality reduction using PCA and UMAP
sc.tl.pca(adata)
sc.pp.neighbors(adata, random_state=1)
sc.tl.umap(adata, random_state=1)
adata.obs[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"].copy()

# Calculate cell cycle phase on lognorm transformed counts
s_phase = txt_to_list(config["s_phase"])
g2m_phase = txt_to_list(config["g2m_phase"])
sc.tl.score_genes_cell_cycle(adata, s_genes=s_phase, g2m_genes=g2m_phase)

# Format observations
obs = adata.obs.copy()
obs[["cell", "tx"]] = obs["sample"].str.split("_", expand=True)
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

# Overlay TooManyCells cluster information and calculate overlay values per node
clusters = pd.read_csv(f"{config['pruned']}/clusters.csv", index_col="cell")
obs["TMC"] = clusters["cluster"].astype(int).astype(str)
node_info = pd.read_csv(f"{config['pruned']}/node_info.csv")
umi_df = node_means(obs, "total_umi", node_info)
diapause_df = node_means(obs, "diapause_score", node_info)

# Bin diapause score and set categorical colour scheme for CELLXGENE
bins = [-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1]
obs["diapause_binned"] = pd.cut(obs["diapause_score"], bins=bins)
bin_order = obs["diapause_binned"].sort_values().unique().astype(str)
obs["diapause_binned"] = (
    obs["diapause_binned"]
    .astype(str)
    .astype("category")
    .cat.reorder_categories(bin_order)
)
adata.uns["diapause_binned_colors"] = [
    "#034e7b",
    "#0570b0",
    "#3690c0",
    "#74a9cf",
    "#a6bddb",
    "#d0d1e6",
    "#f1eef6",
    "#ef8a62",
]
len(obs["diapause_binned"].cat.categories) == len(adata.uns["diapause_binned_colors"])
dict(zip(obs["diapause_binned"].cat.categories, adata.uns["diapause_binned_colors"]))

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
obs["name"] = (
    obs["name"]
    .astype(str)
    .astype("category")
    .cat.reorder_categories(name_cmap["domain"])
)
adata.uns["name_colors"] = name_cmap["range"]

# Add colour scheme for cell cycle phase
phase_cmap = {
    "domain": ["G1", "S", "G2M"],
    "range": [rgb2hex(c) for c in plt.get_cmap("Set1").colors][:3],
}
obs["phase"] = (
    obs["phase"].astype("category").cat.reorder_categories(phase_cmap["domain"])
)
adata.uns["phase_colors"] = phase_cmap["range"]

# Run Harmony
adata.obs = obs.copy()
sc.external.pp.harmony_integrate(adata, key="cell")
sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony")
sc.tl.umap(adata, neighbors_key="harmony")
sc.tl.leiden(adata, neighbors_key="harmony", key_added="leiden_harmony", resolution=0.1)
leiden_harmony_cmap = {
    "domain": [rgb2hex(c) for c in plt.get_cmap("Set1").colors][:8],
    "range": adata.obs["leiden_harmony"].sort_values().unique().astype(str),
}
adata.uns["leiden_harmony_colors"] = leiden_harmony_cmap["domain"]
adata.obs["leiden_harmony"] = adata.obs["leiden_harmony"].astype(str).astype
sc.tl.leiden(adata, neighbors_key="neighbors", resolution=0.1)
adata.obs[["Harmony UMAP1", "Harmony UMAP2"]] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap"] = adata.obs[["UMAP1", "UMAP2"]].to_numpy().copy()

# Save results to file
results_dir = f"{config['results']}/rebuttal"
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata.write_h5ad(f"{results_dir}/adata.h5ad")
for key in ["phase", "name", "leiden", "leiden_harmony", "diapause_binned"]:
    write_adata_label(adata, key, results_dir)
obs.to_csv(f"{results_dir}/obs.csv")
umi_df.to_csv(f"{results_dir}/annotate_umi.csv", index=False, header=True)
diapause_df.to_csv(f"{results_dir}/annotate_diapause.csv", index=False, header=True)

# TO-DO: Export Harmony PCA for use with TooManyCells
pca_df = pd.DataFrame(
    adata.obsm["X_pca_harmony"],
    index=adata.obs_names,
    columns=[f"PC{i}" for i in range(1, 51)],
)
pca_path = os.path.join(results_dir, "pca_harmony/matrix.csv")
pca_df.T.to_csv(pca_path)
