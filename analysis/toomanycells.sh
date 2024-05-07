#!/bin/bash
# Run TooManyCells `make-tree` to generate input files for interactive
# visualization with TMCI

MEX_DIR=./results/download/mex_raw
OUTPUT_DIR=./results/toomanycells
FULL_TREE="$OUTPUT_DIR/01_full_tree"
PRUNED_TREE="$OUTPUT_DIR/02_pruned_tree"

# Create the base tree structure with TF-IDF normalization
mkdir -p $FULL_TREE $PRUNED_TREE
too-many-cells make-tree \
--matrix-path $MEX_DIR \
--labels-file $MEX_DIR/labels_sample.csv \
--filter-thresholds "(250, 1)" \
--normalization "TfIdfNorm" \
--output $FULL_TREE \
--matrix-output $OUTPUT_DIR/mex_tfidf \
> $FULL_TREE/clusters.csv

# Prune TF-IDF normalized tree: Min distance search of 0.019
too-many-cells make-tree \
--prior $FULL_TREE \
--labels-file $MEX_DIR/labels_sample.csv \
--min-distance-search 0.019 \
--draw-collection "NoLeaf" \
--draw-node-number \
--output $PRUNED_TREE \
> $PRUNED_TREE/clusters.csv

# Perform upper quartile normalize on raw counts
too-many-cells matrix-output \
--matrix-path $MEX_DIR \
--filter-thresholds "(250, 1)" \
--normalization "UQNorm" \
--mat-output $OUTPUT_DIR/mex_upperquartile

# Perform quantile normalization on raw counts
too-many-cells matrix-output \
--matrix-path $MEX_DIR \
--filter-thresholds "(250, 1)" \
--normalization "QuantileNorm" \
--mat-output $OUTPUT_DIR/mex_quantilenorm

# Edit features and clusters files to conform to standard 10X Matrix Exchange
# format by replacing NA column with string in features file
norms=(upperquartile quantilenorm)
for norm in ${norms[@]}; do
  gunzip $OUTPUT_DIR/mex_$norm/features.tsv.gz
  sed -i 's/NA/Gene Expression/g' $OUTPUT_DIR/mex_$norm/features.tsv
  gzip $OUTPUT_DIR/mex_$norm/features.tsv
done
exit 0