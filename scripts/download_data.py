"""
Download scRNA-seq data sets and save to file as an AnnData object.
"""

import gzip
import os
import shutil
import tarfile
import tempfile
from pathlib import Path
from urllib import request

import anndata as ad
import pandas as pd
import scanpy as sc

from scripts.config import *


def download_mex_v3(barcodes, features, matrix):
    with tempfile.TemporaryDirectory() as temp_dir:
        with open(temp_dir + "/barcodes.tsv.gz", "wb") as b:
            b.write(request.urlopen(barcodes).read())
        with open(temp_dir + "/features.tsv.gz", "wb") as f:
            f.write(request.urlopen(features).read())
        with open(temp_dir + "/matrix.mtx.gz", "wb") as m:
            m.write(request.urlopen(matrix).read())
        adata = sc.read_10x_mtx(temp_dir)
    return adata


def download_mex_v2(barcodes, genes, matrix):
    with tempfile.TemporaryDirectory() as temp_dir:
        with open(temp_dir + "/barcodes.tsv.gz", "wb") as b:
            b.write(request.urlopen(barcodes).read())
        with open(temp_dir + "/genes.tsv.gz", "wb") as f:
            f.write(request.urlopen(genes).read())
        with open(temp_dir + "/matrix.mtx.gz", "wb") as m:
            m.write(request.urlopen(matrix).read())
        # Create another directory with extracted files
        with tempfile.TemporaryDirectory() as temp_dir2:
            for file in os.listdir(temp_dir):
                if file.endswith(".gz"):
                    with gzip.open(temp_dir + "/" + file, "rb") as zipped:
                        with open(temp_dir2 + "/" + file[:-3], "wb") as unzipped:
                            shutil.copyfileobj(zipped, unzipped)
            adata = sc.read_10x_mtx(temp_dir2)
    return adata


def download_tar(url, suffix):
    def find_parent_dir(dirpath, filename):
        for root, _, files in os.walk(dirpath, topdown=False):
            if any(filename in f for f in files):
                return root
        return None

    with open(tempfile.NamedTemporaryFile(suffix=suffix).name, "wb") as f:
        f.write(request.urlopen(url).read())
    with tempfile.TemporaryDirectory() as temp_dir:
        with tarfile.open(f.name) as tar:
            tar.extractall(temp_dir)
        mtx_dir = find_parent_dir(temp_dir, "matrix.mtx.gz")
        adata = sc.read_10x_mtx(mtx_dir)
    return adata


def adata_to_10x(adata, out_dir):
    Path(out_dir).mkdir(parents=True, exist_ok=True)


# Download read count data
adatas = {}
geo_mex = config["mex"]
for sample, urls in geo_mex.items():
    if "features" in urls:
        # Parse file download URLs from config
        barcodes = urls["url"] + urls["barcodes"]
        features = urls["url"] + urls["features"]
        matrix = urls["url"] + urls["matrix"]
        # Download and save data to AnnData object
        if features.endswith("features.tsv.gz"):
            raw = download_mex_v3(barcodes, features, matrix)
        elif features.endswith("genes.tsv.gz"):
            raw = download_mex_v2(barcodes, features, matrix)
    elif "tar" in urls:
        tar = urls["url"] + urls["tar"]
        raw = download_tar(tar, ".tar.gz")
    # Annotate observations with sample name
    raw.obs[["cell_line", "treatment"]] = sample.split("_")
    adatas[sample] = raw
adata = ad.concat(adatas)

# Subset to only include cells from the DND-41 whitelist
wl = pd.read_csv(config["dnd41_whitelist"], header=None)[0]
dnd41_wl = {
    "dnd41_control": wl[wl.str.endswith("1")],
    "dnd41_short": wl[wl.str.endswith("2")].replace("2", "1"),
    "dnd41_long": wl[wl.str.endswith("3")].replace("3", "1"),
}
for sample, cells in dnd41_wl.items():
    wl_idx = adatas[sample].obs_names.isin(cells)
    adatas[sample] = adatas[sample][wl_idx]

# Write AnnData to H5AD
adata.obs_names_make_unique()
results_dir = os.path.join(config["results"], Path(__file__).stem)
Path(results_dir).mkdir(parents=True, exist_ok=True)
adata.write_h5ad(f"{results_dir}/raw.h5ad")

# Write AnnData in 10X format for downstream use with TooManyCells
