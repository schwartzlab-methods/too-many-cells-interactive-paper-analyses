"""
Provide the application's entry-point
"""

import sys

from scanpy import read_h5ad

from h5ad_to_mex.cli import _run_cli

from .. import adata_to_10x

rc = 1
try:
    h5ad_path, save_dir = _run_cli()
    adata = read_h5ad(h5ad_path)
    adata_to_10x(adata, mex_path=save_dir)
    rc = 0
except Exception as e:
    print("Error: %s" % e, file=sys.stderr)
sys.exit(rc)
