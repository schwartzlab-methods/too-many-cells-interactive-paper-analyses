"""
Single-cell sequencing data aggregation package for Python
"""

import gzip
import os
from pathlib import Path

import pandas as pd
from scipy.io import mmwrite


def adata_to_10x(adata, mex_path):
    Path(mex_path).mkdir(parents=True, exist_ok=True)
    # Transpose read counts to matrix market format, such that (rows, columns)
    # correspond to (features, cells)
    with gzip.open(os.path.join(mex_path, "matrix.mtx.gz"), "wb") as mtx:
        mmwrite(mtx, adata.X.T, symmetry="symmetric")
    # Barcode observations
    barcodes_path = os.path.join(mex_path, "barcodes.tsv.gz")
    adata.obs.index.to_csv(barcodes_path, sep="\t", header=False)
    # Gene features
    features_path = os.path.join(mex_path, "features.tsv.gz")
    var_df = pd.DataFrame(
        {
            "gene_id": adata.var.reset_index()["index"],
            "gene_symbol": adata.var.reset_index()["index"],
            "type": "Gene Expression",
        }
    )
    var_df.to_csv(features_path, sep="\t", header=False, index=False)
