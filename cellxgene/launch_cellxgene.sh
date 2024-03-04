#!/bin/bash

pip install cellxgene

python "./cellxgene/01_prepare_h5ad.py"

cellxgene prepare ./data/cellxgene.h5ad --output ./data/cellxgene-prepared.h5ad

cellxgene launch ./data/cellxgene-prepared.h5ad --open