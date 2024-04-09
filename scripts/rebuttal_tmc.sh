#/bin/bash

DATA_DIR=./results/rebuttal/

too-many-cells make-tree \
--matrix-path $DATA_DIR"pca_harmony/matrix.csv" \
--labels-file $DATA_DIR"labels_name.csv" \
--normalization "NoneNorm" \
--output $DATA_DIR"pca_harmony_tmc" \
> $DATA_DIR"pca_harmony_tmc/clusters.csv"