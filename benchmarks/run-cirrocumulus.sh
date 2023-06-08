#! /usr/bin/env bash

set -euo pipefail

input_path=/code/matrices

python3 mex-to-h5ad.py "${input_path}" cirrocumulus.h5ad

cirro launch --host 0.0.0.0 --port 5004 cirrocumulus.h5ad