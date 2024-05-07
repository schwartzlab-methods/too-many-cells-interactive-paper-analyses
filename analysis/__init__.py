"""
Module for case study analysis of cancer drug-tolerant persister cell line
experiments using TooManyCells Interactive.
"""

from pathlib import Path

import yaml

# Import configs
url_path = Path(__file__).parent / "config.yaml"
with url_path.open(mode="rb") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
