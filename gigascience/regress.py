"""
Implement batch correction of read count data using linear regression
"""

import os
from pathlib import Path

import scanpy as sc
from scipy.sparse import csr_matrix

from gigascience import config
from src.helpers import normalize, reduce_dimensions

# Load raw read counts
results_dir = config["regress"]
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(config["raw_adata"])
adata.obs[["cell", "tx"]] = adata.obs["sample"].str.split("_", expand=True)

# Normalize according to standard Scanpy workflow and regress out batch effects
# between different experiments and save to file
normalize(adata)
sc.pp.regress_out(adata, keys=["cell"], n_jobs=100)
reduce_dimensions(adata)
sc.tl.leiden(adata, resolution=0.1)
adata.layers["regress_cell"] = csr_matrix(adata.X, copy=True, dtype="float32")
adata.X = adata.layers["raw"].copy()
adata.write(os.path.join(results_dir, "adata.h5ad"))

# Export UMAPs and Leiden clusters
adata.obs[["UMAP1", "UMAP2"]] = adata.obsm["X_umap"].copy()
filename = os.path.join(results_dir, "umap_regress.csv")
adata.obs[["UMAP1", "UMAP2", "leiden"]].to_csv(filename)
