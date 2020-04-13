"""UMAP of Library Strategy"""
import sys

import numpy as np
import pandas as pd
import joblib
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

MAPPER = {
    "est": "EST",
    "wgs": "WGS",
    "chip": "ChIP-Seq",
    "outlier": "Outlier",
    "other": "Other",
}

COLORS = {"est": "C1", "wgs": "C2", "chip": "C4", "outlier": "C3", "other": "lightgray"}

SCATTER_STYLE = dict(s=10, edgecolors="w", linewidths=0.2,)


def main():
    style_use(snakemake.params.get("style", "sra"))

    umap = pd.read_parquet(snakemake.input.umap)
    iso = joblib.load(snakemake.input.iso)

    for strategy in ["est", "wgs", "chip"]:
        inliers, outliers = inliner_outlier_split(iso, strategy)
        umap_ = label_outliers(umap, inliers, outliers, strategy)
        grps = umap_.groupby("label")

        _, ax = plt.subplots()
        plot(grps, ax, strategy)
        plot(grps, ax, "outlier")
        plot(grps, ax, "other")

        ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
        sns.despine(ax=ax, left=True, bottom=True)
        plt.legend()
        plt.savefig(snakemake.output[strategy])


def inliner_outlier_split(model, strategy):
    features = pd.read_parquet(snakemake.input[strategy])
    return model.inliers(features).index, model.outliers(features).index


def label_outliers(umap: pd.DataFrame, inliers: pd.Index, outliers: pd.Index, strategy: str):
    inlier_mask = umap.index.isin(inliers)
    outlier_mask = umap.index.isin(outliers)

    df = umap.copy()
    df["label"] = "other"
    df.loc[inlier_mask, "label"] = strategy
    df.loc[outlier_mask, "label"] = "outlier"
    return df


def plot(grps, ax, label):
    df = grps.get_group(label)
    color = COLORS[label]
    label_ = MAPPER[label]

    if label == "other":
        zorder = 0
    elif label == "outlier":
        zorder = 3
    else:
        zorder = 2

    return ax.scatter(df.UMAP1, df.UMAP2, c=color, label=label_, zorder=zorder, **SCATTER_STYLE)


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                umap="../../output/library_strategy-wf/umap_prealn_features.parquet",
                iso="../../output/library_strategy-wf/isolation_forest_model.pkl",
                est="../../output/library_strategy-wf/est_features.parquet",
                wgs="../../output/library_strategy-wf/wgs_features.parquet",
                chip="../../output/library_strategy-wf/chip_features.parquet",
            ),
            output="",
        )

    main()
