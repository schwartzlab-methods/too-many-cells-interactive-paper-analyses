"""
Download single-cell RNA sequencing data from cancer drug-tolerant persister
cell line experiments and aggregate into one file for downstream use with
TooManyCells and differential expresssion analysis.
"""

import os
from pathlib import Path

import pandas as pd
from anndata import concat

from analysis import config
from src.downloads import download_mex_v2, download_mex_v3, download_tar
from src.helpers import adata_to_10x, export_labels

# Load read-count data from GEO repository links, as specified in the config
# file
adatas = {}
for data_type, url_dict in config["data"].items():
    if data_type == "v2":
        for sample, urls in url_dict.items():
            adatas[sample] = download_mex_v2(**urls)
    if data_type == "v3":
        for sample, urls in url_dict.items():
            adatas[sample] = download_mex_v3(**urls)
    if data_type == "tar":
        for sample, urls in url_dict.items():
            adatas[sample] = download_tar(**urls)

# Filter out DND-41 cells not found within the doublet-removal whitelist
dnd = {
    "1": "DND-41 control",
    "2": "DND-41 treated short",
    "3": "DND-41 treated long",
}
whitelist = pd.read_csv(config["whitelist"], header=None)
whitelist[["bc", "id"]] = whitelist[0].str.split("-", expand=True)
for id, sample in dnd.items():
    barcodes = [bc + "-1" for bc in whitelist[whitelist["id"] == id]["bc"].tolist()]
    adatas[sample] = adatas[sample][barcodes].copy()

# Concatenate data from each experiment and format metadata
adata = concat(adatas, label="sample")
obs = adata.obs
# Recode barcode IDs to indicate sample origin
bc = [barcode[:17] for barcode in obs.index]
adata.obs_name = [b + str(i) for b, i in zip(bc, obs["sample"].factorize()[0])]
adata.obs_names_make_unique()

# Annotate cell line name and treatment duration
cols = ["cell", "tx", "duration"]
obs[cols] = obs["sample"].str.split(" ", expand=True)
obs.loc[obs["cell"].isin(["LNCaP", "PC9", "SK-MEL-28"]), "duration"] = "short"
obs.loc[obs["cell"].isin(["MDA-MB-231"]), "duration"] = "long"
obs["condition"] = obs[cols].astype(str).agg("_".join, axis=1).str.lower()

# Preprocess data and save AnnData object to file
adata.layers["raw"] = adata.X.copy()
Path(config["download"]).mkdir(parents=True, exist_ok=True)
adata.write_h5ad(os.path.join(config["download"], "adata.h5ad"))
# Save MEX format files and labels for use with TooManyCells
mex_path = os.path.join(config["download"], "mex_raw")
adata_to_10x(adata, "raw", mex_path)
export_labels(mex_path, adata.obs, "condition")
