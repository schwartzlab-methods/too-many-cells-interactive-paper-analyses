"""
Register new color schemes for rebutal plots
"""

from pathlib import Path

import altair as alt
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import rgb2hex

from scripts.configs import *
from src import altair_themes

alt.themes.register("publishTheme", altair_themes.publishTheme)
alt.themes.enable("publishTheme")
alt.data_transformers.disable_max_rows()
alt.data_transformers.enable("vegafusion")


def save_altair(plot, plot_id, results_dir, img_formats=["svg", "html", "png"]):
    Path(results_dir).mkdir(parents=True, exist_ok=True)
    for img_format in img_formats:
        plot_path = f"{results_dir}/{plot_id}.{img_format}"
        plot.save(plot_path)


# Load observations file
obs = pd.read_csv(f"{config['results']}/rebuttal/obs.csv")
results_dir = f"{config['results']}/rebuttal_plots"
Path(results_dir).mkdir(parents=True, exist_ok=True)

### Reviewer 1
# Define color schemes
paired = [rgb2hex(c) for c in plt.get_cmap("Paired").colors][:10]
paired.insert(2, "#003f94")
name_cmap = {
    "domain": [
        "DND-41 control",
        "DND-41 treated (short)",
        "DND-41 treated (long)",
        "LNCaP control",
        "LNCaP treated",
        "MDA-MB-231 control",
        "MDA-MB-231 treated",
        "PC9 control",
        "PC9 treated",
        "SK-MEL-28 control",
        "SK-MEL-28 treated",
    ],
    "range": paired,
}
phase_cmap = {
    "domain": ["G1", "S", "G2M"],
    "range": [rgb2hex(c) for c in plt.get_cmap("Set1").colors][:3],
}

# Plot umaps
embeddings = ["UMAP", "Harmony UMAP"]
base = alt.Chart(obs).mark_circle(size=2)
for embedding in embeddings:
    points = base.encode(
        alt.X(f"{embedding}1").title(None), alt.Y(f"{embedding}2").title(None)
    )
    embedding_name = embedding.split(" ")[0].lower()
    sample = points.encode(
        alt.Color("name").title("Sample").scale(name_cmap).sort(name_cmap["domain"])
    )
    phase = points.encode(
        alt.Color("phase")
        .title("Cell cycle")
        .scale(phase_cmap)
        .sort(phase_cmap["domain"])
    )
    umap = alt.vconcat(sample, phase).configure_view(fill=None)
    save_altair(umap, embedding_name, results_dir)

### Reviewer 3
cellxgene_cmap_path = f"{config['scrnaseq_data']}/cellxgene_cmap.csv"
cellxgene_cmap = pd.read_csv(cellxgene_cmap_path)["RGB Hex"].to_list()
