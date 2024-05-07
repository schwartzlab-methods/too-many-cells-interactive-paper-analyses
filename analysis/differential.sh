#!bin/bash/
# Run differential expression analysis between cell populations of interest
# using TooManyCells

FULL_TREE=./results/toomanycells/01_full_tree
PRUNED_TREE=./results/toomanycells/02_pruned_tree
OUTPUT_DIR=./results/differential
MEX_PATH=./results/toomanycells/mex_upperquartile
LABELS_DIR=./results/toomanycells/labels

# Control vs. any treatment
mkdir -p $OUTPUT_DIR
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $FULL_TREE \
--labels-file $LABELS_DIR/control_treated.csv \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"control\"], [\"treated\"])" \
> $OUTPUT_DIR/control_vs_treated.csv

# Control vs. short-term treatment
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/control_short.csv \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"control\"], [\"short\"])" \
> $OUTPUT_DIR/control_vs_short.csv

# Control vs. long-term treatment
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/control_long.csv \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"control\"], [\"long\"])" \
> $OUTPUT_DIR/control_vs_long.csv

# PC9 control vs. treated
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/cell_control_treated.csv \
--top-n 100000 \
--nodes "([2], [277])" \
--labels "([\"pc9_control\"], [\"pc9_treated\"])" \
> $OUTPUT_DIR/pc9_control_vs_treated.csv

# SK-MEL-28 control vs. treated
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/cell_control_treated.csv \
--top-n 100000 \
--nodes "([790], [855])" \
--labels "([\"skmel28_control\"], [\"skmel28_treated\"])" \
> $OUTPUT_DIR/skmel28_control_vs_treated.csv

# LNCaP control vs. treated
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/cell_control_treated.csv \
--top-n 100000 \
--nodes "([507], [440])" \
--labels "([\"lncap_control\"], [\"lncap_treated\"])" \
> $OUTPUT_DIR/lncap_control_vs_treated.csv

# DND-41 control vs. short treatment
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/cell_control_short.csv \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"dnd41_control\"], [\"dnd41_short\"])" \
> $OUTPUT_DIR/dnd41_control_vs_short.csv

# DND-41 control vs. long treatment
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/cell_control_long.csv \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"dnd41_control\"], [\"dnd41_long\"])" \
> $OUTPUT_DIR/dnd41_control_vs_long.csv

# MDA-MB-231 control vs. treated
too-many-cells differential \
--matrix-path $MEX_PATH \
--prior $PRUNED_TREE \
--labels-file $LABELS_DIR/cell_control_treated.csv \
--top-n 100000 \
--nodes "([0], [0])" \
--labels "([\"mda231_control\"], [\"mda231_treated\"])" \
> $OUTPUT_DIR/mda231_control_vs_treated.csv
exit 0