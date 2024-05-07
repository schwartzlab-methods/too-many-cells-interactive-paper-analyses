"""
Address reviewer comments and run GSEA on batch-corrected data
"""

import os
from pathlib import Path

import pandas as pd
import scanpy as sc
from gseapy import parser

from gigascience import config
from src.helpers import get_gsea_preranked, get_log2fc, node_idx

# Load raw read counts and annotate samples
results_dir = config["gsea"]
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata = sc.read_h5ad(os.path.join(config["cellxgene"], "adata.h5ad"))
cols = ["names", "logfoldchanges"]
obs = adata.obs

# Add TMC cluster and node_info annotation
clusters = pd.read_csv(config["tmc_clusters"])
node_info = pd.read_csv(config["tmc_node_info"])
obs["tmc"] = clusters["cluster"].astype(int)

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
# Run GSEA on node 48 (ID2 high) vs. 45 (ID2 low)
deg_id2 = get_log2fc(adata, 45, 48, node_info)
deg_id2.to_csv(os.path.join(results_dir, "log2fc_id2.csv"), index=False)
gsea_id2 = get_gsea_preranked(deg_id2[cols], geneset_dict, prerank_params)
gsea_id2.to_csv(os.path.join(results_dir, "gsea_id.csv"), index=False)

# Run GSEA on node 126 vs. 3 for treated breast cancer cell line
subset = adata[adata.obs["name"] == "MDA-MB-231 treated"].copy()
deg_breast = get_log2fc(subset, 3, 126, node_info)
deg_breast.to_csv(os.path.join(results_dir, "log2fc_breast.csv"), index=False)
gsea_breast = get_gsea_preranked(deg_breast[cols], geneset_dict, prerank_params)
gsea_breast.to_csv(os.path.join(results_dir, "gsea_breast.csv"), index=False)

# Run GSEA on TMC PCA Harmony cluster 33448 vs. rest
adata.obs["deg_harmony"] = "0"
adata.obs.loc[node_idx(adata.obs, 33448, node_info, "TMC_harmony"), "deg_harmony"] = (
    "33448"
)
sc.tl.rank_genes_groups(
    adata,
    groupby="deg_harmony",
    groups=["33448"],
    reference="rest",
    layer="lognorm",
    method="wilcoxon",
)
deg_tmc = get_log2fc(adata, "rest", "33448", node_info)
