#/bin/bash

DATA_DIR=./results/rebuttal/

# Run TooManyCells on PCA embedding
mkdir -p $DATA_DIR"pca_tmc"
too-many-cells make-tree \
--matrix-path $DATA_DIR"pca/matrix.csv" \
--labels-file $DATA_DIR"labels_name.csv" \
--normalization "NoneNorm" \
--output $DATA_DIR"pca_tmc" \
> $DATA_DIR"pca_tmc/clusters.csv"

# Run TooManyCells on Harmony-corrected PCA embedding
mkdir -p $DATA_DIR"pca_harmony_tmc"
too-many-cells make-tree \
--matrix-path $DATA_DIR"pca_harmony/matrix.csv" \
--labels-file $DATA_DIR"labels_name.csv" \
--normalization "NoneNorm" \
--output $DATA_DIR"pca_harmony_tmc" \
> $DATA_DIR"pca_harmony_tmc/clusters.csv"

# Run TooManyCells on Harmony-corrected PCA embedding wih lower modularity
mkdir -p $DATA_DIR"pca_harmony_tmc_modularity-1M"
too-many-cells make-tree \
--matrix-path $DATA_DIR"pca_harmony/matrix.csv" \
--labels-file $DATA_DIR"labels_name.csv" \
--min-modularity -1000000 \
--normalization "NoneNorm" \
--output $DATA_DIR"pca_harmony_tmc_modularity-1M" \
> $DATA_DIR"pca_harmony_tmc_modularity-1M/clusters.csv"