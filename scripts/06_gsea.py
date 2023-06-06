"""
Perform differential expression analysis and GSEA on single-cell RNA-seq data
"""

import os
import sys

import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
import matplotlib.pyplot as plt
import altair as alt
from altair_saver import save

import src.altair_themes as altair_themes
from src.helpers import preranked, rank_product, onevsrest 

alt.themes.register("publishTheme", altair_themes.publishTheme)
alt.themes.enable("publishTheme")
alt.data_transformers.disable_max_rows()
plt.rcParams['svg.fonttype'] = 'none'


def main():
    with open("./analysis/config.yaml", "r") as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # Load in results from differential expression analysis between control
    # and resistant cells per cell line and aggregate into a dataframe
    log2fc_dict = {}
    foldchange = pd.DataFrame()
    for log2fc_file in os.listdir(config["diff"]):
        log2fc_path = os.path.join(config["diff"], log2fc_file)
        cell_name = log2fc_file.split("_")[0]
        log2fc_df = pd.read_csv(log2fc_path, index_col=0)
        log2fc_df.replace([np.inf, -np.inf], np.nan, inplace=True)
        log2fc_df.dropna(subset=["log2FC"], how="all", inplace=True)
        log2fc_df = log2fc_df.sort_values("log2FC", ascending=False)
        log2fc_dict[cell_name] = log2fc_df["log2FC"]
    foldchange = pd.concat(log2fc_dict, axis=1)

    # Using table of foldchange values per cell line, calculate the rank
    # product for downstream use with Metascape. This outputs all significantly
    # expressed genes, filtered by p-value.
    foldchange.drop(["shorttx", "longtx", "all"], axis=1, inplace=True)
    foldchange.dropna(axis=0, how="any", inplace=True)
    rank_prod = rank_product(foldchange)
    rank_prod[rank_prod["pVal"] <= 0.05].to_csv(
        os.path.join(config["gsea"], "rankprod_all.csv"), index_label="Gene"
    )
    # Short-term treatment vs. control
    shortterm = ["skmel28", "pc9", "lncap", "dnd41short"]
    rank_prod = rank_product(foldchange[shortterm])
    rank_prod[rank_prod["pVal"] <= 0.05].to_csv(
        os.path.join(config["gsea"], "rankprod_short.csv"), index_label="Gene"
    )
    # Long-term treatment vs. control
    longterm = ["mdamb231", "dnd41long"]
    rank_prod = rank_product(foldchange[longterm])
    rank_prod[rank_prod["pVal"] <= 0.05].to_csv(
        os.path.join(config["gsea"], "rankprod_long.csv"), index_label="Gene"
    )
    
    # Using the log2FC values, run GSEA preranked on each control vs. treated
    # comparison within cell lines using Hallmark, C2, and C6 gene sets from
    # MSigDB
    geneset_dict = {}
    # Populate dictionary with MSigDB gene sets
    for geneset_file in os.listdir(config["gmt_msigdb"]):
        geneset_path = os.path.join(config["gmt_msigdb"], geneset_file)
        geneset_name = geneset_file.split(".")[0]
        geneset_dict[geneset_name] = gp.parser.read_gmt(geneset_path)
    # Populate dictionary with custom diapause and TFEB gene sets from CSV
    plot_dict = {}
    for geneset_file in os.listdir(config["gmt_custom"]):
        geneset_path = os.path.join(config["gmt_custom"], geneset_file)
        geneset_df = pd.read_csv(geneset_path)
        geneset_name = geneset_file.split("_")[0]
        plot_dict[f'{geneset_name}_up'] = list(geneset_df.gene[geneset_df.weight > 0])
        plot_dict[f'{geneset_name}_down'] = list(geneset_df.gene[geneset_df.weight < 0])
    plot_dict['G2M'] = geneset_dict['c2']['FISCHER_G2_M_CELL_CYCLE']
    geneset_dict['custom'] = plot_dict
            
    # Run GSEA
    for cell_name, log2fc in log2fc_dict.items():
        gsea_dict = {}
        for geneset_name, geneset in geneset_dict.items():
            res = preranked(log2fc, geneset)
            res_df = pd.DataFrame(res.res2d)
            res_df["Geneset"] = geneset_name
            res_df.sort_values("NES", ascending=True, inplace=True)
            gsea_dict[geneset_name] = res_df
            # Plot GSEA results of Fischer G2M cell cycle geneset
            if geneset_name == "custom":
                terms = res_df.Term
                for term in terms:
                    plotname = os.path.join(
                        config["gsea"],
                        f'{cell_name}_{geneset_name}_{term}.svg',
                    )
                    gp.gseaplot(
                        rank_metric=res.ranking,
                        term=term,
                        ofname=plotname,
                        **res.results[term],
                    )
        gsea_df = pd.concat(gsea_dict, axis=0, ignore_index=True)
        gsea_path = os.path.join(config["gsea"], f"{cell_name}_gsea.csv")
        gsea_df.sort_values("NES", ascending=False).to_csv(gsea_path)
        # Create scatterplot of FDR vs. NES
        plot_df = gsea_df[gsea_df.Geneset != 'custom']
        gsea_plot = (
            alt.Chart(plot_df)
            .mark_line(
                color="black",
                point=alt.OverlayMarkDef(size=40, color='black')
            )
            .encode(
                x=alt.X("NES:Q"),
                y=alt.Y("FDR q-val:Q"),
                tooltip=["Term", "NES", "FDR q-val"],
            )
            .properties(width=100, height=100)
        )
        g2m_plot = alt.Chart(plot_df).mark_circle(
                color="red",
                size=100
            ).encode(
                x=alt.X("NES:Q"),
                y=alt.Y("FDR q-val:Q"),
                tooltip=["Term", "NES", "FDR q-val"],
            ).properties(
                width=100, height=100
            ).transform_filter(
                alt.FieldEqualPredicate(field='Term', equal='FISCHER_G2_M_CELL_CYCLE')
            )
        final_plot = alt.layer(gsea_plot, g2m_plot, data=plot_df).facet(column='Geneset')
        for img_format in ["png", "svg", "html"]:
            save(
                final_plot,
                os.path.join(config["gsea"], f"{cell_name}_fdr_nes.{img_format}"),
            )

    # Aggregate one vs. all results and export into annotation table format
    # for the TooManyCells Interactive tool
    adata_uq = sc.read_10x_mtx(config['mtx_uq'])
    adata_uq.layers['uqnorm'] = adata_uq.X
    adata = sc.read_h5ad(config['adata'])
    adata_uq.obs = adata_uq.obs.join(adata.obs)
    adata_uq.write_h5ad(config['adata_uq'])
    # Perform one vs. rest differential expression analysis
    # TO-DO: Avoid writing out .csv to file
    onevsall_fc = onevsrest( adata_uq,
        layer="uqnorm",
        cluster_path=os.path.join(config["pruned"], "clusters.csv"),
        save_dir=config["onevsall_fc"],
    )
    onevsall_dict = {}
    for node_id, df in onevsall_fc.groupby("node_id"):
        # Perform Hallmark, C2, and C6 GSEA on each one vs. rest comparison
        for geneset_name, geneset in geneset_dict.items():
            res = preranked(log2fc, geneset)
            res_df = pd.DataFrame(res.res2d)
            res_df["geneset"] = geneset_name
            res_df["node_id"] = node_id
            res_df.sort_values("NES", ascending=True, inplace=True)
            onevsall_dict[node_id] = res_df
    # Concatenate diapause geneset GSEA results and save to file
    onevsall_df = pd.concat(onevsall_dict, axis=0, ignore_index=True)
    custom_df = onevsall_df[onevsall_df.Geneset == "custom"]
    for term, df in custom_df.groupby("Term"):
        filename = os.path.join(config["onevsall_gsea"], f"overlay_{term}.csv")
        df[["node_id", "NES"]].sort_values("node_id").to_csv(filename, index=False)


if __name__ == "__main__":
    sys.exit(main())
