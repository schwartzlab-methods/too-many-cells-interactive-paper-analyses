#! /usr/bin/env bash

set -euo pipefail

tar -cJvf benchmarks.tar.xz bootstrap-instance-cellx.sh Dockerfile.cellxgene mex-to-h5ad.py run-cellxgene.sh bench-cellxgene.sh .dockerignore

rsync benchmarks.tar.xz spot:/home/ubuntu
