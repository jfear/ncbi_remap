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


SELECTED_STRATEGY = ["RNA-Seq", "EST", "WGS", "ChIP-Seq", "Other"]
LIBRARY_STRATEGY_COLORS = ["C0", "C1", "C2", "C4", "lightgray"]
RNASEQ_OUTLIER_COLORS = ["lightgray", "C0", "C3"]


def main():
    umap = pd.read_parquet(snakemake.input.umap)

    style_use(snakemake.params.get("style", "sra"))
    plt.rcParams["figure.figsize"] = (8, 8)

    _, ax1 = plt.subplots()
    ax2 = inset_axes(ax1, 3, 3, bbox_to_anchor=[0, 0, 0.4, 0.4], bbox_transform=ax1.transAxes)

    umap_rnaseq_outliers(umap, ax1)
    umap_library_strategy(umap, ax2)
    add_legend(ax1)

    plt.savefig(snakemake.output[0])


def umap_rnaseq_outliers(umap: pd.DataFrame, ax: plt.Axes):
    iso = joblib.load(snakemake.input.iso)  # type: SraIsolationForest

    umap_other = umap[~umap.index.isin(iso.index)]
    ax.scatter(
        umap_other.UMAP1, umap_other.UMAP2, zorder=0, c=RNASEQ_OUTLIER_COLORS[0], rasterized=True
    )

    umap_rnaseq = umap[umap.index.isin(iso.inliers_all.index)]
    ax.scatter(
        umap_rnaseq.UMAP1, umap_rnaseq.UMAP2, zorder=1, c=RNASEQ_OUTLIER_COLORS[1], rasterized=True
    )

    umap_outlier = umap[umap.index.isin(iso.outliers_all.index)]
    ax.scatter(
        umap_outlier.UMAP1,
        umap_outlier.UMAP2,
        zorder=2,
        c=RNASEQ_OUTLIER_COLORS[2],
        rasterized=True,
    )
    sns.despine(ax=ax, left=True, bottom=True)

    return clean_umap_axis(ax)


def umap_library_strategy(umap: pd.DataFrame, ax: plt.Axes):
    labels = pd.read_parquet(snakemake.input.labels).library_strategy.squeeze()
    labels[~labels.isin(SELECTED_STRATEGY)] = "Other"  # Focus on selected strategies

    umap_ = umap.join(labels)
    grps = umap_.groupby("library_strategy")

    for strategy, color in list(zip(SELECTED_STRATEGY, LIBRARY_STRATEGY_COLORS))[::-1]:
        df = grps.get_group(strategy)
        ax.scatter(df.UMAP1, df.UMAP2, c=color, label=strategy, rasterized=True)

    # Add boarder around inset
    for loc in ["left", "top", "right", "bottom"]:
        ax.spines[loc].set_visible(True)

    return clean_umap_axis(ax, labels=False)


def clean_umap_axis(ax: plt.Axes, labels=True):
    ax.set(xticks=[], yticks=[])
    if labels:
        ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    return ax


def add_legend(ax: plt.Axes):
    strategies = SELECTED_STRATEGY.copy()
    strategies.insert(1, "RNA-Seq Outliers")

    colors = LIBRARY_STRATEGY_COLORS.copy()
    colors.insert(1, RNASEQ_OUTLIER_COLORS[2])

    legend_elements = [
        Line2D(
            [0], [0], marker="o", color=color, label=strategy, markerfacecolor=color, linestyle=""
        )
        for strategy, color in zip(strategies, colors)
    ]

    ax.legend(
        handles=legend_elements,
        loc="upper left",
        bbox_to_anchor=(0, 0.8, 0.1, 0.1),
        bbox_transform=ax.transAxes,
    )


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
