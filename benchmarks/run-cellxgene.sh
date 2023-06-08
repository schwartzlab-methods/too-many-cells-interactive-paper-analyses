#! /usr/bin/env bash

set -euo pipefail

input_path=/code/matrices

python3 mex-to-h5ad.py "${input_path}" cellxgene.h5ad

cellxgene prepare cellxgene.h5ad --overwrite --output cellxgene-prepared.h5ad

cellxgene launch --host 0.0.0.0 --port 5004 cellxgene-prepared.h5ad