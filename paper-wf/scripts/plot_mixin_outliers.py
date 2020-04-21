"""Plot Mixin Experiment Results"""
import sys

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

STRATEGIES = ["EST", "WGS", "ChIP-Seq"]
COLORS = ["C1", "C2", "C4"]

def main():
    style_use("sra")
    plt.rcParams["figure.figsize"] = 3.335, 2.509

    mixin = (
        pd.read_parquet(snakemake.input[0])
        .set_index(["Prop Contamination", "Contamination Source"])
        .unstack()
    )
    mixin.columns = mixin.columns.droplevel()

    ax = mixin[STRATEGIES].plot(
        marker="o", color=["C1", "C2", "C4"], legend=False
    )  # type: plt.Axes
    ax.set(xlabel="Contamination (prop)", ylabel="Outliers (prop)", yticks=[0.2, .4, .6, .8])
    add_legend(ax)
    # sns.despine(ax=ax, left=True, bottom=True)

    plt.savefig(snakemake.output[0])

def add_legend(ax):
    legend_elements = [
        Line2D(
            [0], [0], marker="o", markersize=3, color=color, label=strategy, markerfacecolor=color, linestyle=""
        )
        for strategy, color in zip(STRATEGIES, COLORS)
    ]

    ax.legend(
        handles=legend_elements,
        loc="upper left",
        # bbox_to_anchor=(0, 0.8, 0.1, 0.1),
        # bbox_transform=ax.transAxes,
    )




if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(input="../../output/library_strategy-wf/contamination_mixin.parquet")

    main()
