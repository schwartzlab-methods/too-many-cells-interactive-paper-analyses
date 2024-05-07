"""
Submodule for additional analyses of GigaScience submission
"""

from pathlib import Path

import yaml

# Import configs
config_path = Path(__file__).parent / "config.yaml"
with config_path.open(mode="rb") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
