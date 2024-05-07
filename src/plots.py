"""
Utilites for plotting data
"""

from pathlib import Path
import os


def format_obs(obs, species="human"):
    # Calculate QC metrics
    obs["norm_umi"] = obs["total_umi"] / obs[f"{species}_reads"]
    obs["read_depth"] = 1 / obs["norm_umi"]
    obs["hk_umi"] = obs["pct_hk"] * obs["total_umi"]
    obs["norm_hk"] = obs["hk_umi"] / obs[f"{species}_reads"]
    return obs


def plot_qc_metrics(obs):
    pass


def save_altair(plot, plot_id, results_dir, img_formats=["svg", "html", "png"]):
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    for img_format in img_formats:
        plot_path = os.path.join(results_dir, f"{plot_id}.{img_format}")
        plot.save(plot_path)