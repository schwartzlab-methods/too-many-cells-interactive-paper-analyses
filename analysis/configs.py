"""
Load configs from `config.yaml` file for use in other scripts.
"""

import yaml

with open("./scripts/configs.yaml", "r", encoding="utf-8") as file:
    config = yaml.load(file, Loader=yaml.FullLoader)
