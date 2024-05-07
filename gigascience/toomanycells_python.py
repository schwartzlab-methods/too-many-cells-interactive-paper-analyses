"""
Perform TooManyCells clustering on batch-corrected read counts data
"""

import os

import scanpy as sc
from toomanycells import TooManyCells as tmc

from gigascience import config

# Run TooManyCells on raw counts data
adata = sc.read_h5ad(os.path.join(config["regress"], "adata.h5ad"))
results_dir = config["toomanycells"]
adata.X = adata.layers["regress_cell"].copy()
tmc_regress = tmc(adata, os.path.join(results_dir, "output_regress_v16"))
tmc_regress.run_spectral_clustering()
tmc_regress.store_outputs()
