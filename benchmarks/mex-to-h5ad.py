from sys import argv
import os

import anndata as ad
from scanpy import read_10x_mtx


def prepare(mtx_dir_path: str, outpath: str):

    con = None

    for root, dirs, files in os.walk(mtx_dir_path):
        for dir in dirs:
            adat = read_10x_mtx(os.path.join(root, dir))
            if con:
                con = ad.concat([con, adat])
            else:
                con = adat

    con.write(outpath)


if __name__ == "__main__":
    prepare(argv[1], argv[2])
