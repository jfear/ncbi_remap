"""UMAP of Library Strategy"""
import sys

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

CATEGORIES = ["RNA-Seq", "Outlier", "Other"]
COLORS = ["C0", "C3", "lightgray"]
ZORDER = [1, 2, 0]
SCATTER_STYLE = dict(s=10, edgecolors="w", linewidths=0.2,)


def main():
    style_use(snakemake.params.get("style", "sra"))
    embeddings = wrangle_data()
    ax = plot(embeddings)
    plt.savefig(snakemake.output[0])


def wrangle_data():
    labels = (
        pd.read_parquet(snakemake.input.labels)
        .library_strategy.squeeze()
        .map(lambda x: x if x in CATEGORIES else "Other")
    )

    iso = joblib.load(snakemake.input.iso) # type: ncbi_remap.iforest.SraIsolationForest
    labels[labels.index.isin(iso.outliers_all.index)] = "Outlier"

    return pd.read_parquet(snakemake.input.umap).join(labels)


def plot(embeddings):
    for cat, color, zorder in zip(CATEGORIES, COLORS, ZORDER):
        df = embeddings.query(f"library_strategy == '{cat}'")
        plt.scatter(df.UMAP1, df.UMAP2, c=color, label=cat, zorder=zorder, **SCATTER_STYLE)

    ax = plt.gca()
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    sns.despine(ax=ax, left=True, bottom=True)
    plt.legend()
    return ax


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                umap="../../output/library_strategy-wf/umap_prealn_features_embeddings.parquet",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
                iso="../../output/library_strategy-wf/isolation_forest.pkl"
            ),
            output="",
        )

    main()
