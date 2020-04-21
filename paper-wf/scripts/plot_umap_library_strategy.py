"""Plot the UMAP panel for library strategy"""
import sys

import pandas as pd
import joblib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

CATEGORIES = ["RNA-Seq", "RNA-Seq Outlier", "EST", "WGS", "ChIP-Seq", "Other"]
COLORS = ["C0", "C3", "C1", "C2", "C4", "lightgray"]
ZORDERS = [4, 5, 3, 2, 1, 0]
MARKERS = ["o", "^", "o", "o", "o"]


def main():
    labels = get_labels()
    umap = pd.read_parquet(snakemake.input.umap).join(labels)

    style_use("sra")
    plt.rcParams["figure.figsize"] = (3.335, 2.509)

    _, ax = plt.subplots() # type: plt.Figure, plt.Axes
    grps = umap.groupby("library_strategy")
    for strategy, color, marker, zorder in zip(CATEGORIES, COLORS, MARKERS, ZORDERS):
        _df = grps.get_group(strategy)
        ax.scatter(
            _df.UMAP1,
            _df.UMAP2,
            color=color,
            marker=marker,
            label=strategy,
            zorder=zorder,
            rasterized=True,
            **snakemake.params[0]
        )
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", yticks=[-10, 0], xticks=[-10, 0, 10])
    sns.despine(ax=ax, left=True, bottom=True)
    plt.legend(loc="upper left", bbox_to_anchor=(-.1, 1))
    # ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    # ax.yaxis.set_major_locator(plt.MaxNLocator(4))

    plt.savefig(snakemake.output[0])


def get_labels():
    labels = pd.read_parquet(snakemake.input.labels).library_strategy.squeeze()

    # Focus on a subset of library strategies.
    labels[~labels.isin(CATEGORIES)] = "Other"

    # Flag outliers
    iso = joblib.load(snakemake.input.iso)  # type: SraIsolationForest
    labels[labels.index.isin(iso.outliers_all.index)] = "RNA-Seq Outlier"

    return labels


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                umap="../../output/library_strategy-wf/umap_prealn_features.parquet",
                iso="../../output/library_strategy-wf/isolation_forest.pkl",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            ),
        )

    main()
