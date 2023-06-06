#!/bin/bash
DATA_DIR=./data
RESULT_DIR=./results
PRIOR_DIR="$RESULT_DIR/01_full_tree"
PRUNED_DIR="$RESULT_DIR/02_pruned_tree"

# Create the base tree structure with TF-IDF normalization
mkdir $RESULT_DIR
mkdir $PRIOR_DIR
mkdir $PRUNED_DIR
too-many-cells make-tree \
--matrix-path "$DATA_DIR/mtx_raw" \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--filter-thresholds "(250, 1)" \
--normalization "TfIdfNorm" \
--output "$RESULT_DIR/01_full_tree" \
--matrix-output "$DATA_DIR/mtx_tfidf" \
> "$RESULT_DIR/01_full_tree/clusters.csv"

# Prune TF-IDF normalized tree: Min distance search of 0.019
too-many-cells make-tree \
--prior "$RESULT_DIR/01_full_tree" \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--min-distance-search 0.019 \
--draw-collection "NoLeaf" \
--draw-node-number \
--output "$RESULT_DIR/02_pruned_tree" \
> "$RESULT_DIR/02_pruned_tree/clusters.csv"

# Overlay diapause signature score
mkdir "$RESULT_DIR/02_pruned_tree/diapause"
too-many-cells make-tree \
--prior "$RESULT_DIR/02_pruned_tree" \
--matrix-path "$DATA_DIR/mtx_diapause.csv" \
--matrix-transpose \
--draw-leaf "DrawItem (DrawThresholdContinuous [(\"diapause_score\", MadMedian 1)])" \
--draw-colors "[\"#FF0000\", \"#D3D3D3\"]" \
--output "$RESULT_DIR/02_pruned_tree/diapause" \
> "$RESULT_DIR/02_pruned_tree/diapause/clusters.csv"

# Upper quartile normalization
too-many-cells matrix-output \
--matrix-path "$DATA_DIR/mtx_raw" \
--filter-thresholds "(250, 1)" \
--normalization "UQNorm" \
--mat-output "$DATA_DIR/mtx_upperquartile"

# Quantile normalization
too-many-cells matrix-output \
--matrix-path "$DATA_DIR/mtx_raw" \
--filter-thresholds "(250, 1)" \
--normalization "QuantileNorm" \
--mat-output "$DATA_DIR/mtx_quantilenorm"

### Edit features and clusters files to become compatible with Python scripts
# Replace NA column with string in features file
gunzip "$DATA_DIR/mtx_tfidf/features.tsv.gz"
sed -i 's/NA/Gene Expression/g' "$DATA_DIR/mtx_tfidf/features.tsv"
gzip "$DATA_DIR/mtx_tfidf/features.tsv"

gunzip "$DATA_DIR/mtx_upperquartile/features.tsv.gz"
sed -i 's/NA/Gene Expression/g' "$DATA_DIR/mtx_upperquartile/features.tsv"
gzip "$DATA_DIR/mtx_upperquartile/features.tsv"

gunzip "$DATA_DIR/mtx_quantilenorm/features.tsv.gz"
sed -i 's/NA/Gene Expression/g' "$DATA_DIR/mtx_quantilenorm/features.tsv"
gzip "$DATA_DIR/mtx_quantilenorm/features.tsv"