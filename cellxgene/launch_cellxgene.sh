#!/bin/bash

pip install cellxgene

python "./cellxgene/01_mtx_to_h5ad.py"

cellxgene prepare ./data/cellxgene.h5ad --output ./data/cellxgene-prepared.h5ad

python "./cellxgene/02_annotate.py"

cellxgene launch ./data/cellxgene-annotated.h5ad --open