"""
Integrate batch correction and TooManyCells module with TMCI workflow
"""

import os

import anndata as ad
import scanpy as sc
from toomanycells import TooManyCells as tmc

from gigascience import config
from src.prep import normalize

# Load raw counts data and normalize
raw = sc.read_h5ad(config["raw_adata"])
adatas = {}
adatas["raw"] = ad.AnnData(raw.layers["raw"], obs=raw.obs, var=raw.var)
normalize(raw)
batch_key = "cell"

# Generate Harmony batch corrected embedding
hdata = raw.copy()
sc.pp.scale(hdata)
sc.tl.pca(hdata)
sc.external.pp.harmony_integrate(hdata, key=batch_key, basis="X_pca")
adatas["harmony_pca"] = ad.AnnData(
    hdata.obsm["X_pca_harmony"],
    obs=hdata.obs,
    var=["PC" + str(i) for i in range(1, 51)],
)

# Generate MNN batch corrected embedding
mdata = raw.copy()
sc.pp.highly_variable_genes(mdata, batch_key=batch_key)
var_idx = mdata.var.highly_variable_nbatches > 1
var_genes = var_idx[var_idx].index
subset = {k: mdata[idx] for k, idx in mdata.obs.groupby(batch_key).groups.items()}
mdata = sc.external.pp.mnn_correct(*subset, batch_key=batch_key)
adatas["mnn"] = ad.AnnData(mdata[0][:, var_genes], obs=mdata.obs, var=var_genes)

# Generate Combat batch corrected embedding
cdata = raw.copy()
sc.pp.combat(cdata, key=batch_key, covariates="tx")

# Run TooManyCells on raw counts data
adata = sc.read_h5ad(os.path.join(config["regress"], "adata.h5ad"))
results_dir = config["toomanycells"]
tmc_raw = tmc(adata, os.path.join(results_dir, "output_raw"))
tmc_raw.run_spectral_clustering()
tmc_raw.store_outputs()
tmc_raw.A.obs[["sp_cluster", "sp_path"]]

# Regress out cell types and take absolute value of counts
adata.X = adata.layers["regress_cell"].copy()
tmc_regress = tmc(adata, os.path.join(results_dir, "output_regress_v16"))
tmc_regress.run_spectral_clustering()
tmc_regress.store_outputs()
tmc_regress.A.obs[["sp_cluster", "sp_path"]]
