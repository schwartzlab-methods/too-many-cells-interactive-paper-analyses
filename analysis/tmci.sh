#!/bin/bash
# Visualize TooManyCells outputs using TMCI

START_SCRIPT=/home/christie/tmci/start-and-load.sh
MEX_DIR=/home/christie/tmci-paper-analyses/results/download
TREE_DIR=/home/christie/tmci-paper-analyses/results/toomanycells

sudo $START_SCRIPT \
--matrix-dir "$TREE_DIR/mex_tfidf" \
--tree-path "$TREE_DIR/01_full_tree/cluster_tree.json" \
--label-path "$MEX_DIR/labels_sample.csv" \
--port 3000

sudo $START_SCRIPT \
--matrix-dir "$TREE_DIR/mex_tfidf" \
--tree-path "$TREE_DIR/02_pruned_tree/cluster_tree.json" \
--label-path "$MEX_DIR/labels_sample.csv" \
--port 3001