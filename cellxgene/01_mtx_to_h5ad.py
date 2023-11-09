"""
Visualize case study scRNA-seq data using CELLxGENE Annotate for comparison.
"""

import os
import sys
import yaml

import scanpy as sc


def main():
    with open("./scripts/config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    adata = sc.read_10x_mtx(config['mtx_tfidf'])
    adata_path = os.path.join(config['scrnaseq_data'], 'cellxgene.h5ad')
    adata.write_h5ad(adata_path)
    
    
if __name__ == '__main__':
    sys.exit(main())
