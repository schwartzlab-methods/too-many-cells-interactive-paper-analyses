"""
Download scRNA-seq data sets and save to file as an AnnData object.
"""

import tarfile
import tempfile
from urllib import request

import scanpy as sc

from scripts.config import *


def tar_to_adata(filename):
    with tempfile.TemporaryDirectory() as temp_mex:
        with tarfile.open(filename) as tar:
            tar.extractall(path=temp_mex)
            mtx_dir = find_parent_dir(temp_mex, "matrix.mtx")
            adata = sc.read_10x_mtx(mtx_dir)
    return adata


def download_url(url, suffix=".tar.gz"):
    with open(tempfile.NamedTemporaryFile(suffix=suffix).name, "wb") as f:
        opener = urllib.request.build_opener()
        opener.addheaders = [("User-Agent", "Mozilla/5.0")]
        urllib.request.install_opener(opener)
        with urllib.request.urlopen(url) as handle:
            f.write(handle.read())
            return f


adata_dict = {}
for celltype, url in config["pbmc_urls"].items():
    with open(download_url(url).name) as f:
        adata_dict[celltype] = tar_to_adata(f.name)
