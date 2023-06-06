#!/bin/bash

DATA_DIR=./data
RESULT_DIR=./results
PRIOR_DIR="$RESULT_DIR/01_full_tree"
PRUNED_DIR="$RESULT_DIR/02_pruned_tree"
DIFF_DIR="$RESULT_DIR/03_differential"
MTX="mtx_upperquartile"


### Run differential expression analysis on bulk
# ctrl vs. resistant
mkdir $DIFF_DIR
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/labels_treatment.csv" \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"control\"], [\"treated\"])" \
> "$DIFF_DIR/all_ctrl_vs_res.csv"

# ctrl vs. short-term treatment
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRUNED_DIR \
--labels-file "$DATA_DIR/labels_short_tx.csv" \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"short_control\"], [\"short_tx\"])" \
> "$DIFF_DIR/shorttx_ctrl_vs_res.csv"

# ctrl vs. long-term treatment
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRUNED_DIR \
--labels-file "$DATA_DIR/labels_long_tx.csv" \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"long_control\"], [\"long_tx\"])" \
> "$DIFF_DIR/longtx_ctrl_vs_res.csv"

### Differential expression analysis per cell line
# NSCLC PC9 control vs. resistant
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--top-n 100000 \
--nodes "([2], [277])" \
--labels "([\"pc9_control\"], [\"pc9\"])" \
> "$DIFF_DIR/pc9_ctrl_vs_res.csv"

# Melanoma SK-MEL-28 control vs. resistant
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--top-n 100000 \
--nodes "([790], [855])" \
--labels "([\"skmel28_control\"], [\"skmel28\"])" \
> "$DIFF_DIR/skmel28_ctrl_vs_res.csv"

# Prostate LNCaP control vs. resistant
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--top-n 100000 \
--nodes "([507], [440])" \
--labels "([\"lncap_control\"], [\"lncap\"])" \
> "$DIFF_DIR/lncap_ctrl_vs_res.csv"

# T-ALL DND-41 control vs. resistant
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--top-n 100000 \
--nodes "([613], [711])" \
--labels "([\"dnd41_control\"], [\"dnd41_short\"])" \
> "$DIFF_DIR/dnd41short_ctrl_vs_res.csv"

too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--top-n 100000 \
--nodes "([613], [754])" \
--labels "([\"dnd41_control\"], [\"dnd41_long\"])" \
> "$DIFF_DIR/dnd41long_ctrl_vs_res.csv"

# Breast MDA-MB-231 control vs. resistant
too-many-cells differential \
--matrix-path "$DATA_DIR/$MTX" \
--prior $PRIOR_DIR \
--labels-file "$DATA_DIR/mtx_raw/labels.csv" \
--top-n 100000 \
--nodes "([911], [990])" \
--labels "([\"mdamb231_control\"], [\"mdamb231\"])" \
> "$DIFF_DIR/mdamb231_ctrl_vs_res.csv"