"""
Add annotations to prepared CELLxGENE read count matrix. 
"""

import os
import sys
import yaml

import pandas as pd
import scanpy as sc


def main():
    with open("./scripts/config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    adata_path = os.path.join(config['scrnaseq_data'], 'cellxgene-prepared.h5ad')
    adata = sc.read_h5ad(adata_path)
    
    # Add cell line and treatment condition information
    labels_path = os.path.join(config['mtx_raw'], 'labels.csv')
    adata.obs['cell_line'] = pd.read_csv(labels_path, index_col='item')['label']
    scores_path = os.path.join(config['scrnaseq_data'], 'mtx_diapause.csv')
    
    # Add diapause score information
    adata.obs['diapause'] = pd.read_csv(scores_path, index_col=0)['diapause_score']
    adata_path = os.path.join(config['scrnaseq_data'], 'cellxgene-annotated.h5ad')
    
    # Add cluster information from pruned TooManyCells tree
    tmc_path = os.path.join(config['pruned'], 'clusters.csv')
    tmc_clusters = pd.read_csv(tmc_path, index_col='cell')['cluster'].copy()
    adata.obs['toomanycells'] = tmc_clusters.astype(str)
    adata.write_h5ad(adata_path)
    
    
if __name__ == '__main__':
    sys.exit(main())
