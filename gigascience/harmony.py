"""
Implement batch correction of read count data using Harmony
"""

import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from gigascience import config
from src.helpers import normalize

results_dir = config["harmony"]
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(config["raw_adata"])
adata.obs[["cell", "tx"]] = adata.obs["sample"].str.split("_", expand=True)
batch_key = "cell"

# Normalize data and run Harmony batch correction
normalize(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)
sc.external.pp.harmony_integrate(adata, key=batch_key, basis="X_pca")
sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony")
sc.tl.umap(adata, neighbors_key="harmony")
sc.tl.leiden(adata, neighbors_key="harmony", resolution=0.1)

# Export PCAs for use with TooManyCells
pca_harmony = pd.DataFrame(
    adata.obsm["X_pca_harmony"],
    index=adata.obs_names,
    columns=[f"PC{i}" for i in range(1, 51)],
)
pca_harmony.T.to_csv(os.path.join(results_dir, "pca_harmony.csv"))

# Export UMAPs and Leiden clusters
adata.obs[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"].copy()
filename = os.path.join(results_dir, "umap_harmony.csv")
adata.obs[["UMAP1", "UMAP2", "leiden"]].to_csv(filename)
