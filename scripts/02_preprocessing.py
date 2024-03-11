"""Aggregate datasets and perform gene signature scoring."""

import os
import sys

import anndata as ad
import pandas as pd
import scanpy as sc
import yaml
from scipy.sparse import csr_matrix

from src.helpers import adata_to_10x, export_labels, signature_score


def main():
    with open("./analysis/config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # Read counts in the data folder are read into the mtx_dict dictionary
    # and annotated with the sample ID prior to export in 10X format
    mtx_dict = {}
    for mtx_file in os.listdir(config["raw_data"]):
        mtx_celltype = mtx_file.split(".")[0]
        mtx_path = os.path.join(config["raw_data"], mtx_file)
        # Matrix files can be either txt or 10X format
        if mtx_path.endswith(".txt.gz"):
            mtx_dict[mtx_celltype] = ad.read_text(mtx_path).transpose()
        else:
            mtx_dict[mtx_celltype] = sc.read_10x_mtx(mtx_path)
        # Add observation metadata to AnnData object
        obs = pd.DataFrame(index=mtx_dict[mtx_celltype].obs_names)
        obs["sample"] = mtx_celltype
        mtx_dict[mtx_celltype].obs = obs

    # DND-41 cell line data has a whitelist to remove doublets
    whitelist = pd.read_csv("./src/whitelist_dnd41.csv", header=None)
    # Long-term treated cells have suffix '-3' instead of '-1'
    long = mtx_dict["dnd41_long"].copy()
    long.obs_names = [x.replace("-1", "-3") for x in long.obs_names]
    mtx_dict["dnd41_long"] = long[long.obs_names.intersection(whitelist[0])]
    # Short-term treated cells have suffix '-2' instead of '-1'
    short = mtx_dict["dnd41_short"].copy()
    short.obs_names = [x.replace("-1", "-2") for x in short.obs_names]
    mtx_dict["dnd41_short"] = short[short.obs_names.intersection(whitelist[0])]
    # Control untreated cells have suffix of '-1'
    ctrl = mtx_dict["dnd41_control"].copy()
    mtx_dict["dnd41_control"] = ctrl[ctrl.obs_names.intersection(whitelist[0])]

    # Calculate diapause expression scores for each cell using average of
    # z-scaled gene counts and export to file for use with TooManyCells
    adata = ad.concat(mtx_dict, index_unique="-", label="sample")
    adata.layers["raw"] = adata.X
    diapause_path = "./analysis/diapause_genes.csv"
    signature_score(adata, geneset_path=diapause_path, geneset_id="diapause")
    adata.obs["diapause_score"].to_csv(
        os.path.join(config["scrnaseq_data"], "mtx_diapause.csv")
    )

    # Diapause scores need to be convered into 10x matrix market format to use
    # with TooManyCells interactive tool
    diapause_scores = adata.obs[["diapause_score"]]
    diapause_adata = sc.AnnData(
        X=diapause_scores,
        obs=diapause_scores.index.to_frame(),
        var=diapause_scores.columns.to_frame(),
    )
    diapause_adata.X = csr_matrix(diapause_adata.X)
    adata_to_10x(
        diapause_adata, obs_key="mtx_diapause", data_dir=config["scrnaseq_data"]
    )

    # Export treatment labels to use with TooManyCells differential expression
    # analysis
    adata.obs["treatment"] = [
        "control" if x.endswith("control") else "treated" for x in adata.obs["sample"]
    ]
    export_labels(config["scrnaseq_data"], adata.obs, "treatment")
    # Create labels for GSEA within short-term tx cell lines
    adata.obs["tx"] = [
        "control" if x.endswith("control") else "short_tx" for x in adata.obs["sample"]
    ]
    adata.obs.loc[adata.obs["sample"] == "dnd41_long", "tx"] = "long_tx"
    adata.obs.loc[adata.obs["sample"] == "mdamb231", "tx"] = "long_tx"
    adata.obs["short_tx"] = adata.obs["tx"]
    short_cells = ["dnd41_control", "lncap_control", "pc9_control", "skmel28_control"]
    adata.obs.loc[adata.obs["sample"].isin(short_cells), "short_tx"] = "short_control"
    export_labels(config["scrnaseq_data"], adata.obs, "short_tx")
    # Create labels for GSEA within long-term tx cell lines
    adata.obs["long_tx"] = adata.obs["tx"]
    long_cells = ["dnd41_control", "mdamb231_control"]
    adata.obs.loc[adata.obs["sample"].isin(long_cells), "long_tx"] = "long_control"
    export_labels(config["scrnaseq_data"], adata.obs, "long_tx")
    del adata.obs["tx"]

    # Export AnnData object to file for use with TooManyCells
    adata_to_10x(adata, obs_key="mtx_raw", data_dir=config["scrnaseq_data"])
    adata.write_h5ad(config["adata"])


if __name__ == "__main__":
    sys.exit(main())
