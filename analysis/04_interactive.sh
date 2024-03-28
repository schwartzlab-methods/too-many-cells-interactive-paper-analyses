#!/bin/bash

DATA_DIR=/home/christie/tmci-paper-analyses/data
RESULT_DIR=/home/christie/tmci-paper-analyses/results

# Full tree with gene expression matrix overlay
sudo /home/christie/too-many-cells-interactive/start-and-load.sh \
--matrix-dir "$DATA_DIR/mtx_tfidf" \
--tree-path "$RESULT_DIR/01_full_tree/cluster_tree.json" \
--label-path "$DATA_DIR/mtx_raw/labels.csv" \
--port 3000

sudo /home/christie/too-many-cells-interactive/start-and-load.sh \
--matrix-dir "$DATA_DIR/mtx_tfidf" \
--tree-path "$RESULT_DIR/02_pruned_tree/cluster_tree.json" \
--label-path "$DATA_DIR/mtx_raw/labels.csv" \
--port 3000

# Full tree with diapause scores overlay
sudo /home/christie/too-many-cells-interactive/start-and-load.sh \
--matrix-dir "$DATA_DIR/mtx_diapause" \
--tree-path "$RESULT_DIR/01_full_tree/cluster_tree.json" \
--label-path "$DATA_DIR/mtx_raw/labels.csv" \
--port 3001

sudo /home/christie/too-many-cells-interactive/start-and-load.sh \
--matrix-dir "$DATA_DIR/mtx_diapause" \
--tree-path "$RESULT_DIR/02_pruned_tree/cluster_tree.json" \
--label-path "$DATA_DIR/mtx_raw/labels.csv" \
--port 3001

# Run `ssh -L 3001:localhost:3001 -L 3000:localhost:3000 HOSTNAME` locally to
# bind ports and access the graphical interface on your browser