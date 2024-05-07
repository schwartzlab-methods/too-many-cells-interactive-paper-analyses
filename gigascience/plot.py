"""
Plot figures for GigaScience reviewer comments
"""

import os

import altair as alt
import pandas as pd
import scanpy as sc
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests

from gigascience import config
from src import altair_themes
from src.plots import save_altair

alt.themes.register("publishTheme", altair_themes.publishTheme)
alt.themes.enable("publishTheme")
alt.data_transformers.disable_max_rows()


def u_test(obs, name, name_key, tx, tx_key="Treatment", value_key="ID2"):
    df = pd.DataFrame()
    for cell in name:
        x = obs[(obs[name_key] == cell) & (obs[tx_key] == tx[0])][value_key]
        y = obs[(obs[name_key] == cell) & (obs[tx_key] == tx[1])][value_key]
        h, h_pval = kruskal(x, y)
        u, u_pval = mannwhitneyu(x, y, alternative="less")
        df.loc[cell, ["H", "H_pval", "U", "U_pval"]] = h, h_pval, u, u_pval
    df["U_pval_adj"] = multipletests(df["U_pval"], alpha=0.05, method="fdr_bh")[1]
    return df


# Load data
uq = sc.read_10x_mtx(config["uq_mex"])
obs = uq.obs
obs["name"] = pd.read_csv(config["name"], index_col="item")["label"]
obs[["cell", "tx"]] = obs["name"].str.split(" ", n=1, expand=True)
short = ["DND-41", "LNCaP", "PC9", "SK-MEL-28"]
long = ["DND-41", "MDA-MB-231"]
short_mask = obs["cell"].isin(short) & obs["tx"].isin(["treated", "treated (short)"])
long_mask = obs["cell"].isin(long) & obs["tx"].isin(["treated", "treated (long)"])
obs["Treatment"] = "Control"
obs.loc[short_mask, "Treatment"] = "Short tx"
obs.loc[long_mask, "Treatment"] = "Long tx"
obs["diapause"] = pd.read_csv(config["diapause"], index_col="item")["label"]
obs["ID2"] = uq[:, "ID2"].X.A

# Plot cell counts
cell_order = ["DND-41", "LNCaP", "MDA-MB-231", "PC9", "SK-MEL-28"]
tx_order = ["Control", "Short tx", "Long tx"]
counts = (
    alt.Chart(obs)
    .mark_bar()
    .encode(
        alt.X("cell").title("Cell line").axis(labelAngle=-45).sort(cell_order),
        alt.XOffset("Treatment").sort(tx_order),
        alt.Y("count()").title("Number of cells"),
        alt.Color("Treatment").sort(tx_order).scale(scheme="yellowgreenblue"),
    )
)
save_altair(counts, "fig4a_counts", config["plots"])

# Create boxplot of ID2 expression and calculate stats
id2 = obs[obs["ID2"] > 0]
id2_plot = (
    alt.Chart(id2)
    .mark_boxplot(
        extent="min-max",
        outliers=False,
        size=5,
        median=alt.MarkConfig(color="black"),
    )
    .encode(
        alt.X("cell").title("Cell line").axis(labelAngle=-45).sort(cell_order),
        alt.XOffset("Treatment").sort(tx_order),
        alt.Y("ID2").title("ID2 expression").scale(type="symlog"),
        alt.Color("Treatment").sort(tx_order).scale(scheme="yellowgreenblue"),
    )
)
save_altair(id2_plot, "fig4d_id2", config["plots"])
id2_short = u_test(id2, short, "cell", ["Control", "Short tx"])
id2_long = u_test(id2, long, "cell", ["Control", "Long tx"])
id2_u = pd.concat([id2_short, id2_long], keys=["Short tx", "Long tx"])
id2_u.to_csv(os.path.join(config["plots"], "fig4d_id2_stats.csv"))

# Create boxplot of diapause expression
diapause = (
    alt.Chart(obs)
    .mark_boxplot(
        extent="min-max",
        outliers=False,
        size=5,
        median=alt.MarkConfig(color="black"),
    )
    .encode(
        alt.X("cell").title("Cell line").axis(labelAngle=-45).sort(cell_order),
        alt.XOffset("Treatment").sort(tx_order),
        alt.Y("diapause").title("Diapause score"),
        alt.Color("Treatment").sort(tx_order).scale(scheme="yellowgreenblue"),
    )
)
save_altair(diapause, "fig4k_diapause", config["plots"])
obs.groupby(["cell", "Treatment"])["diapause"].describe()
d_short = u_test(obs, short, "cell", ["Control", "Short tx"], value_key="diapause")
d_long = u_test(obs, long, "cell", ["Control", "Long tx"], value_key="diapause")
d_u = pd.concat([d_short, d_long], keys=["Short tx", "Long tx"])
d_u.to_csv(os.path.join(config["plots"], "fig4k_diapause_stats.csv"))
