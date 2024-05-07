"""
Specify category colors in H5AD file for use in CellxGene
"""

import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from gigascience import config
from src.helpers import (
    convert_to_cat,
    interpolate_colors,
    node_means,
    normalize,
    plot_umap,
    reduce_dimensions,
    save_altair,
    txt_to_list,
    write_adata_label,
)

# Load raw read counts and normalize according to standard Scanpy workflow
results_dir = config["cellxgene"]
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(config["raw_adata"])
obs = adata.obs
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
normalize(adata)
reduce_dimensions(adata)

# Cluster using uncorrected neighbours and plot UMAP
sc.tl.leiden(adata, resolution=0.1)
clusters = adata.obs["leiden"].unique().sort_values().astype(str)
convert_to_cat(adata.obs, col="leiden", order=clusters)
leiden_cmap = {
    "domain": clusters,
    "range": interpolate_colors(clusters, palette="tab10"),
}
adata.uns["leiden_colors"] = leiden_cmap["range"]
umap = plot_umap(adata.obs, c="leiden", ctitle="Leiden", cmap=leiden_cmap)
save_altair(umap, plot_id=f"umap_leiden", results_dir=results_dir)

# Calculate cell cycle phase on lognorm transformed counts
adata.X = adata.layers["lognorm"].copy()
s = txt_to_list(config["s_phase"])
g2m = txt_to_list(config["g2m_phase"])
sc.tl.score_genes_cell_cycle(adata, s_genes=s, g2m_genes=g2m)
phases = ["G1", "G2M", "S"]
convert_to_cat(adata.obs, col="phase", order=phases)
adata.uns["phase_colors"] = interpolate_colors(phases, palette="Set1")

# Overlay TooManyCells cluster information and calculate diapause means per node
tmc_clusters = pd.read_csv(config["tmc_clusters"], index_col="cell")
obs["tmc"] = tmc_clusters["cluster"].astype(int).astype(str)
node_info = pd.read_csv(config["tmc_node_info"])
diapause_df = node_means(obs, obs_key="diapause_score", node_info=node_info)
diapause_path = os.path.join(results_dir, "tmci_diapause.csv")
diapause_df.to_csv(diapause_path, index=False, header=True)

# Bin diapause score
bins = list(range(0, 12))
obs["diapause_binned"] = pd.cut(obs["diapause_score"], len(bins), labels=bins)
hexmap = interpolate_colors(len(bins) - 1, palette="YlGnBu")[::-1]
adata.uns["diapause_binned_colors"] = hexmap + ["#ff0000"]

# Add colour scheme for sample conditions
names = [
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
]
paired = interpolate_colors(len(names) - 1, palette="Paired")
paired.insert(2, "#003f94")
convert_to_cat(obs, col="name", order=names)
adata.uns["name_colors"] = paired

# Add Harmony corrected UMAP and Leiden clustering
harmony_filename = os.path.join(config["harmony"], "umap_harmony.csv")
harmony = pd.read_csv(harmony_filename, index_col=0)
adata.obsm["X_umap_harmony"] = harmony[["UMAP1", "UMAP2"]].values
adata.obs["leiden_harmony"] = harmony["leiden"].copy()
harmony_clusters = adata.obs["leiden_harmony"].sort_values().unique().astype(str)
convert_to_cat(adata.obs, "leiden_harmony", harmony_clusters)
adata.uns["leiden_harmony_colors"] = interpolate_colors(harmony_clusters, "tab10")

# Add regression corrected UMAP and Leiden clustering
regress_filename = os.path.join(config["regress"], "umap_regress.csv")
regress = pd.read_csv(regress_filename, index_col=0)
adata.obsm["X_umap_regress"] = regress[["UMAP1", "UMAP2"]].values
adata.obs["leiden_regress"] = regress["leiden"].copy()
regress_clusters = adata.obs["leiden_regress"].sort_values().unique().astype(str)
convert_to_cat(adata.obs, "leiden_regress", regress_clusters)
adata.uns["leiden_regress_colors"] = interpolate_colors(regress_clusters, "tab10")

# Save results to file
adata.obs = obs.copy()
adata.write_h5ad(os.path.join(results_dir, "adata.h5ad"))
for key in [
    "phase",
    "name",
    "leiden",
    "leiden_harmony",
    "leiden_regress",
    "diapause_binned",
    "diapause_score",
]:
    write_adata_label(adata, key, results_dir)
