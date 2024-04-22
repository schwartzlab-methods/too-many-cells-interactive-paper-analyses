"""
Create the application's command-line interface
"""

import argparse


def _run_cli(argv=None):
    # Parse the command-line inputs to the application
    parser = argparse.ArgumentParser(
        prog="h5ad_to_mex",
        description="Convert local H5AD read counts file to 10X-format MEX",
    )
    parser.add_argument(
        "--h5ad_path",
        nargs="*",
        action="append",
        help="Path to local read counts H5AD file",
    )
    parser.add_argument(
        "--save_dir",
        nargs="*",
        action="append",
        help="Output directory for writing the MEX files",
    )
    args = parser.parse_args(argv)
    return args.h5ad_path, args.save_dir
