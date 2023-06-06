"""Plot figures."""

import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from altair_saver import save

import altair as alt
from src import altair_themes
from src.helpers import standardize, signature_score, mannwhitney

alt.themes.register("publishTheme", altair_themes.publishTheme)
alt.themes.enable("publishTheme")
alt.data_transformers.disable_max_rows()


def main():
    with open("./analysis/config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # Read in TF-IDF normalized AnnData object and annotate observations
    adata = sc.read_h5ad(config["adata"])
    tfidf = sc.read_10x_mtx(config["mtx_tfidf"])
    tfidf.obs = adata[tfidf.obs.index].obs.copy()
    tfidf.write_h5ad(config["adata_tfidf"])
    uqnorm = sc.read_h5ad(config["adata_uq"])
    tfidf.layers["uqnorm"] = uqnorm.layers["uqnorm"]
    # Read in cluster annotations from pruned tree
    clusters = pd.read_csv(
        os.path.join(config["full"], "clusters.csv"), index_col="cell"
    )
    clusters_dict = clusters["path"].to_dict()
    node_dict = {}
    for barcode, path in clusters_dict.items():
        path_list = [int(x) for x in path.split("/")]
        for node in path_list:
            if node not in node_dict:
                node_dict[node] = [barcode]
            else:
                node_dict[node].append(barcode)
    main_branch_dict = {
        2: "pc9_control",
        277: "pc9",
        790: "skmel28_control",
        855: "skmel28",
        507: "lncap_control",
        440: "lncap",
        613: "dnd41_control",
        711: "dnd41_short",
        754: "dnd41_long",
        911: "mdamb231_control",
        990: "mdamb231",
    }
    tfidf.obs["main_branch"] = None
    tfidf.obs["cell"] = [x.split("_")[0] for x in tfidf.obs["sample"]]
    for node_id, cell_line in main_branch_dict.items():
        tfidf.obs.loc[
            (tfidf.obs_names.isin(node_dict[node_id]))
            & (tfidf.obs["sample"] == cell_line),
            "main_branch",
        ] = cell_line
    tfidf.obs["cell"] = tfidf.obs["cell"].str.upper()
    tfidf.obs["tx"] = tfidf.obs["treatment"].str.capitalize()
    tfidf.obs.loc[tfidf.obs.long_tx == "long_tx", "tx"] = "Long tx"
    tfidf.obs.loc[tfidf.obs.long_tx == "short_tx", "tx"] = "Short tx"
    obs = tfidf.obs.copy()

    # Fig A: Plot number of cells per celltype
    tx_order = ["Control", "Short tx", "Long tx"]
    cell_order = ["LNCAP", "SKMEL28", "PC9", "MDAMB231", "DND41"]
    cell_count_plot = (
        alt.Chart(obs)
        .mark_bar()
        .encode(
            x=alt.X(
                "cell:N",
                sort=cell_order,
                title="Cell line",
                axis=alt.Axis(labelAngle=-45),
            ),
            xOffset=alt.XOffset("tx:N", sort=tx_order),
            y=alt.Y("count()", title="Number of cells"),
            color=alt.Color("tx:N", sort=tx_order, title="Treatment"),
        )
        .properties(height=150, width=150)
    )
    for img_format in [".png", ".svg", ".html"]:
        name = os.path.join(config["plots"], f"total_cells_plot{img_format}")
        save(cell_count_plot, name)

    # Fig D: Boxplot of ID2 expression per celltype
    obs["id2"] = uqnorm[obs.index, "ID2"].to_df().values
    id2_df = obs[obs["id2"] > 0].copy()
    id2_plot = (
        alt.LayerChart(id2_df)
        .transform_aggregate(
            min="min(id2)",
            max="max(id2)",
            q1="q1(id2)",
            median="median(id2)",
            q3="q3(id2)",
            groupby=["cell", "tx"],
        )
        .encode(
            x=alt.X(
                "cell:N",
                sort=cell_order,
                title="Cell line",
                axis=alt.Axis(labelAngle=-45),
            ),
            y=alt.Y(title="ID2 expression", scale=alt.Scale(type="symlog")),
            xOffset=alt.XOffset("tx:N", sort=tx_order),
            tooltip=["cell:N", "tx:N", "min:Q", "max:Q", "mean:Q", "q1:Q", "q3:Q"],
        )
        .properties(width=200, height=150)
        .add_layers(
            alt.Chart().mark_rule().encode(x="cell:N", y="min:Q", y2="max:Q"),
            alt.Chart()
            .mark_bar(width=10)
            .encode(
                x="cell:N",
                y="q1:Q",
                y2="q3:Q",
                color=alt.Color("tx:N", title="Treatment", sort=tx_order),
            ),
            alt.Chart().mark_tick(color="black", width=10).encode(y="median:Q"),
        )
    )
    for img_format in [".png", ".svg", ".html"]:
        name = os.path.join(config["plots"], f"id2_boxplot_dropzero_median{img_format}")
        save(id2_plot, name)
    # Calculate stats and write to file
    id2_stats = mannwhitney(id2_df, "id2")
    id2_stats.to_csv(
        os.path.join(config["plots"], "id2_stats_dropzero.csv"), index=False
    )

    # Fig J: Boxplot of diapause scores per celltype from major nodes
    signature_score(adata, geneset_path=config['diapause'], geneset_id='diapause')
    obs['diapause_score'] = adata.obs['diapause_score']
    diapause_plot = alt.LayerChart(obs).transform_aggregate(
        min='min(diapause_score)',
        max='max(diapause_score)',
        q1='q1(diapause_score)',
        median='median(diapause_score)',
        q3='q3(diapause_score)',
        groupby=['cell', 'tx']
    ).encode(
        x=alt.X('cell:N', sort=cell_order, title='Cell line', axis=alt.Axis(labelAngle=-45)),
        y=alt.Y(title='Diapause score'),
        xOffset=alt.XOffset('tx:N', sort=tx_order),
        tooltip=['cell:N', 'tx:N', 'min:Q', 'max:Q', 'mean:Q', 'q1:Q', 'q3:Q']
    ).properties(width=200, height=150
    ).add_layers(
        alt.Chart().mark_rule().encode(x='cell:N', y='min:Q', y2='max:Q'),
        alt.Chart().mark_bar(width=10).encode(
            x='cell:N',
            y='q1:Q',
            y2='q3:Q',
            color=alt.Color('tx:N', title='Treatment', sort=tx_order)
        ),
        alt.Chart().mark_tick(color='black', width=10).encode(y='median:Q')
    )
    for img_format in [".png", ".svg", ".html"]:
        name = os.path.join(config["plots"], f"diapause_boxplot{img_format}")
        save(diapause_plot, name)
    diapause_stats = mannwhitney(adata.obs, "diapause_score")
    diapause_stats.to_csv(
        os.path.join(config["plots"], "diapause_stats.csv"), index=False
    )

    # Create boxplot of total read counts per cell
    obs["total_counts"] = adata[obs.index].X.sum(axis=1)
    counts_plot = (
        alt.Chart(obs)
        .mark_boxplot()
        .encode(
            x=alt.X(
                "cell:N",
                title="Cell line",
                axis=alt.Axis(labelAngle=-45),
                sort=cell_order,
            ),
            xOffset=alt.XOffset("tx:N", sort=tx_order),
            y=alt.Y("total_counts:Q", title="Total reads per cell"),
            color=alt.Color("tx:N", sort=tx_order, title="Treatment"),
            order="order:Q",
        )
        .properties(width=300)
        .configure_view(stroke=None)
    )
    for img_format in [".png", ".svg", ".html"]:
        name = os.path.join(config["plots"], f"read_counts_plot{img_format}")
        save(counts_plot, name)

    # Separate out upregulated and downregulated diapause signature genes and
    # plot to troubleshoot directionality of GSEA
    diapause_path = os.path.join(config["gmt_custom"], "diapause_genes.csv")
    geneset = pd.read_csv(diapause_path)
    genes_up = [gene for gene in geneset[geneset["weight"] > 0]["gene"]]
    genes_down = [gene for gene in geneset[geneset["weight"] < 0]["gene"]]
    znorm = standardize(adata.layers["raw"])
    reads = pd.DataFrame(
        znorm.todense(), index=adata.obs_names, columns=adata.var_names
    )
    adata_up = reads[reads.columns.intersection(genes_up)]
    adata_down = reads[reads.columns.intersection(genes_down)]
    adata.obs["score_up"] = np.round(adata_up.mean(axis=1), 3)
    adata.obs["score_down"] = np.round(adata_down.mean(axis=1), 3)
    # Plot diapause up scores
    diapause_up = (
        alt.Chart(adata.obs)
        .mark_boxplot()
        .encode(
            x=alt.X(
                "cell:N",
                title="Cell line",
                axis=alt.Axis(labelAngle=-45),
                sort=cell_order,
            ),
            xOffset=alt.XOffset("tx:N", sort=tx_order),
            y=alt.Y("score_up:Q", title="Diapause score (upregulated genes)"),
            color=alt.Color("tx:N", sort=tx_order, title="Treatment"),
            order="order:Q",
        )
        .properties(width=300)
        .configure_view(stroke=None)
    )
    for img_format in [".png", ".svg", ".html"]:
        save(
            diapause_up,
            os.path.join(config["plots"], f"diapause_up_boxplot{img_format}"),
        )
    # Plot diapause down scores
    diapause_down = (
        alt.Chart(adata.obs)
        .mark_boxplot()
        .encode(
            x=alt.X(
                "cell:N",
                title="Cell line",
                axis=alt.Axis(labelAngle=-45),
                sort=cell_order,
            ),
            xOffset=alt.XOffset("tx:N", sort=tx_order),
            y=alt.Y("score_down:Q", title="Diapause score (upregulated genes)"),
            color=alt.Color("tx:N", sort=tx_order, title="Treatment"),
            order="order:Q",
        )
        .properties(width=300)
        .configure_view(stroke=None)
    )
    for img_format in [".png", ".svg", ".html"]:
        save(
            diapause_down,
            os.path.join(config["plots"], f"diapause_down_boxplot{img_format}"),
        )


if __name__ == "__main__":
    sys.exit(main())
