#/bin/bash
# Run TooManyCells `make-tree` on batch-corrected read counts.


HARMONY_DIR=./results/gigascience/harmony/
REGRESS_DIR=./results/gigascience/regress/
LABELS_DIR=./results/gigascience/cellxgene/
OUT_DIR=./results/gigascience/toomanycells/


# Make tree from regression-corrected read counts 
mkdir -p $OUT_DIR"regress"
too-many-cells make-tree \
--matrix-path $REGRESS_DIR"mex" \
--labels-file $LABELS_DIR"labels_name.csv" \
--normalization "NoneNorm" \
--output $OUT_DIR"regress" \
--shift-positive \
--dense \
> $OUT_DIR"regress/clusters.csv"

# Make tree from Harmony-corrected PCA embedding wih lower modularity
mkdir -p $OUT_DIR"pca_harmony"
too-many-cells make-tree \
--matrix-path $HARMONY_DIR"pca_harmony.csv" \
--labels-file $LABELS_DIR"labels_name.csv" \
--min-modularity -100000000000000000000000000000000000000000000000000000000000 \
--normalization "NoneNorm" \
--shift-positive \
--output $OUT_DIR"pca_harmony" \
--dense \
> $OUT_DIR"pca_harmony/clusters.csv"

exit 0