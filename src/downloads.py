"""
Functions for downloading and reading single-cell cell
"""

import gzip
import os
import shutil
import tarfile
from tempfile import NamedTemporaryFile, TemporaryDirectory
from urllib.request import urlopen

from scanpy import read_10x_mtx

from src.helpers import find_parent_dir


def download_tar(url, read_func=read_10x_mtx, suffix=".tar.gz", **kwargs):
    """
    Download compressed 10X Matrix Exchange directory to AnnData.
    """
    with NamedTemporaryFile(suffix=suffix) as f:
        f.write(urlopen(url).read())
        with TemporaryDirectory() as temp_dir:
            with tarfile.open(f.name) as tar:
                tar.extractall(temp_dir, filter=None)
            mtx_dir = find_parent_dir(temp_dir, "matrix.mtx.gz")
            adata = read_func(mtx_dir)
    return adata


def download_mex_v3(barcodes, features, matrix):
    """
    Download 10X CellRanger v3 feature-barcode matrix files to AnnData.
    """
    with TemporaryDirectory() as temp_dir:
        with open(temp_dir + "/barcodes.tsv.gz", "wb") as b:
            b.write(urlopen(barcodes).read())
        with open(temp_dir + "/features.tsv.gz", "wb") as f:
            f.write(urlopen(features).read())
        with open(temp_dir + "/matrix.mtx.gz", "wb") as m:
            m.write(urlopen(matrix).read())
        adata = read_10x_mtx(temp_dir)
    return adata


def download_mex_v2(barcodes, genes, matrix):
    """
    Download 10X CellRanger v2 gene-barcode matrix files to Anndata.
    """
    with TemporaryDirectory() as temp_dir:
        with open(temp_dir + "/barcodes.tsv.gz", "wb") as b:
            b.write(urlopen(barcodes).read())
        with open(temp_dir + "/genes.tsv.gz", "wb") as f:
            f.write(urlopen(genes).read())
        with open(temp_dir + "/matrix.mtx.gz", "wb") as m:
            m.write(urlopen(matrix).read())
        # Unzip .gz files and read into Python
        for file in os.listdir(temp_dir):
            with gzip.open(temp_dir + "/" + file, "rb") as zipped:
                with open(temp_dir + "/" + file[:-3], "wb") as unzipped:
                    shutil.copyfileobj(zipped, unzipped)
        adata = read_10x_mtx(temp_dir)
    return adata
