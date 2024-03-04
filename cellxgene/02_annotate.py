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
    
    # Load and compare to TMCI output
    tfidf = sc.read_10x_mtx(config['mtx_tfidf'])
    diapause = sc.read_10x_mtx(config['mtx_diapause'])
    cluster_path = os.path.join(config['pruned'], 'diapause/clusters.csv')
    df = tfidf[:, ['ID2', 'MTOR']].to_df()
    df['cluster'] = pd.read_csv(cluster_path, index_col='cell')['cluster']
    df['diapause'] = diapause.to_df()['diapause_score']
    df.groupby('cluster').mean().describe()
    df.describe()
    
    # Load CELLxGENE prepared AnnData file
    cg_path = os.path.join(config['scrnaseq_data'], 'cellxgene-prepared2.h5ad')
    cgdata = sc.read_h5ad(cg_path)
    cfdf = cgdata[:, ['ID2', 'MTOR']].to_df()
    cfdf['cluster'] = pd.read_csv(cluster_path, index_col='cell')['cluster']
    cfdf['diapause'] = diapause.to_df()['diapause_score']
    cfdf.groupby('cluster').mean().describe()
    cfdf.describe()
    
    cgp_path = os.path.join(config['scrnaseq_data'], 'cellxgene-prepared.h5ad')
    cgpdata = sc.read_h5ad(cg_path)
    
if __name__ == '__main__':
    sys.exit(main())
