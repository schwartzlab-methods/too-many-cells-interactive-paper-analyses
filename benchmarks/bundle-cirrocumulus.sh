#! /usr/bin/env bash

set -euo pipefail

tar -cJvf benchmarks.tar.xz bootstrap-instance-cirro.sh Dockerfile.cirrocumulus mex-to-h5ad.py run-cirrocumulus.sh bench-cirro.sh .dockerignore

rsync benchmarks.tar.xz spot:/home/ubuntu
