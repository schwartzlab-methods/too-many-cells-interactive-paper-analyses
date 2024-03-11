"""
Address reviewer comments and prepare rebuttal.
"""

from pathlib import Path

import pandas as pd
import scanpy as sc

from scripts.configs import *
from src.helpers import write_adata_label


def txt_to_list(txt_path):
    with open(txt_path, "r") as f:
        return f.read().splitlines()


def node_means(obs, obs_key, cluster_key="TMC"):
    node_info = pd.read_csv(f"{config['pruned']}/node_info.csv")
    subtree_dict = dict(zip(node_info["node"], node_info["subtree"].str.split("/")))
    vals = {}
    for node_id, subtree in subtree_dict.items():
        vals[node_id] = obs[obs_key][obs[cluster_key].isin(subtree)].mean()
    vals_df = pd.DataFrame().from_dict(vals, orient="index")
    vals_df = vals_df.reset_index().rename(columns={"index": "node_id", 0: obs_key})
    return vals_df


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
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata.obs[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"].copy()

# Calculate cell cycle phase
s_phase = txt_to_list(config["s_phase"])
g2m_phase = txt_to_list(config["g2m_phase"])
sc.tl.score_genes_cell_cycle(adata, s_genes=s_phase, g2m_genes=g2m_phase)

# Run Harmony
adata.obs[["cell", "tx"]] = adata.obs["sample"].str.split("_", expand=True)
sc.external.pp.harmony_integrate(adata, key="cell")
sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony")
sc.tl.umap(adata, neighbors_key="harmony")
adata.obs[["Harmony UMAP1", "Harmony UMAP2"]] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap_harmony"] = adata.obsm["X_umap"].copy()
adata.obsm["X_umap"] = adata.obs[["UMAP1", "UMAP2"]].to_numpy().copy()

# Format observations
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

# Overlay TooManyCells cluster information and calculate overlay values per node
clusters = pd.read_csv(f"{config['pruned']}/clusters.csv", index_col="cell")
adata.obs["TMC"] = clusters["cluster"].astype(int).astype(str)
umi_df = node_means(adata.obs, obs_key="total_umi")
diapause_df = node_means(adata.obs, obs_key="diapause_score")

# Save results to file
results_dir = f"{config['results']}/rebuttal"
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata.obs = obs.copy()
adata.write_h5ad(f"{results_dir}/adata.h5ad")
for key in ["phase", "name"]:
    write_adata_label(adata, key, results_dir)
obs.to_csv(f"{results_dir}/obs.csv")
umi_df.to_csv(f"{results_dir}/annotate_umi.csv", index=False, header=True)
diapause_df.to_csv(f"{results_dir}/annotate_diapause.csv", index=False, header=True)
