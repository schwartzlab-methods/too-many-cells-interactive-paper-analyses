# TooManyCellsInteractive Paper Analysis

### Scripts
- _01_download_data.sh_: Download scRNA-seq read count data from GEO data.
- _02_preprocessing.py_: Aggregate scRNA-seq datasets across cell lines into AnnData.
- _03_toomanycells.sh_: Run TooManyCells to generate dendrograms and clusters. 
- _04_interactive.sh_: Run TMCI to generate interactive web page for dendrograms.
- _05_differential.sh_: Run differential expression analysis on scRNA-seq data.
- _06_gsea.py_: Run rank product and gene set enrichment analysis on scRNA-seq data.
- _07_plot.py_: Plot manuscript figures using Altair library.
- _config.yaml_: Configuration file containing relative paths for analysis. 

### Source
- _gene_sets_msigdb_: Contains `.gmt` files that describe genes in a given gene set, sourced from MSigDB.
- _altair_themes.py_: Defines plot configurations for Altair library. 
- _diapause_genes.csv_: List of diapause signature genes and corresponding weights, sourced from Rehman et al, Cell 2021. 
- _helpers.py_: Defines helper functions for case study analysis. 
- _whitelist_dnd41.csv_: Defines whitelist of cell barcodes for DND-41 cells after doublet removal. 
