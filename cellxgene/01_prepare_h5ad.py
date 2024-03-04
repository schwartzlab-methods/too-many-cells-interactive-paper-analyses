"""
Visualize case study scRNA-seq data using CELLxGENE Annotate for comparison.
"""

import os
import sys
import yaml

import pandas as pd
import scanpy as sc


def main():
    with open("./scripts/config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    
    # Prepare UMAP embeddings from log-normalized read counts
    rdata = sc.read_10x_mtx(config['mtx_raw'])
    rdata.layers['raw'] = rdata.X.copy()
    sc.pp.normalize_total(rdata, target_sum=10000)
    rdata.layers['norm'] = rdata.X.copy()
    sc.pp.log1p(rdata)
    rdata.layers['lognorm'] = rdata.X.copy()
    sc.tl.pca(rdata)
    sc.pp.neighbors(rdata)
    sc.tl.leiden(rdata)
    sc.tl.umap(rdata)
    
    # Load in TF-IDF normalized read counts from TooManyCells output. Add the
    # UMAP embedding layer and observation annotations.
    adata = sc.read_10x_mtx(config['mtx_tfidf'])
    adata.obsm['X_pca'] = rdata.obsm['X_pca'].copy()
    adata.obsm['X_umap'] = rdata.obsm['X_umap'].copy()
    adata.obs['Leiden clusters'] = rdata.obs['leiden'].copy()
    
    # Add cell line and treatment condition information
    labels_path = os.path.join(config['mtx_raw'], 'labels.csv')
    adata.obs['Cell line'] = pd.read_csv(labels_path, index_col='item')['label']
    
    # Add diapause score information
    scores_path = os.path.join(config['scrnaseq_data'], 'mtx_diapause.csv')
    adata.obs['Diapause score'] = pd.read_csv(scores_path, index_col=0)['diapause_score']
    
    # Add cluster information from pruned TooManyCells tree
    tmc_path = os.path.join(config['pruned'], 'clusters.csv')
    tmc_clusters = pd.read_csv(tmc_path, index_col='cell')['cluster'].copy()
    adata.obs['TooManyCells clusters'] = tmc_clusters.astype(str)
    
    # Write to file
    adata_path = os.path.join(config['scrnaseq_data'], 'cellxgene.h5ad')
    adata.write_h5ad(adata_path)
    
    
if __name__ == '__main__':
    sys.exit(main())
