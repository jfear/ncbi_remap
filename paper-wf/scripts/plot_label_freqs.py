"""Plot library strategy representation in the SRA.

Bar plot showing the percent of samples for each library strategy.

"""
import os
from dateutil.parser import parse as date_parse

import numpy as np
import pandas as pd
from pymongo import MongoClient

import matplotlib.pyplot as plt
import seaborn as sns


def main():
    pct_labels = get_label_percentages()

    fig = plt.figure(figsize=plt.figaspect(2))
    ax = sns.barplot(
        "pct_samples", "Label", data=pct_labels, color="lightgray", edgecolor="k"
    )  # type: plt.Axes
    ax.set(
        xlim=(0, 40),
        xlabel="Percent of Samples",
        ylabel="SRA Label",
        title="Library Strategy Representation",
    )
    sns.despine(ax=ax, left=True, bottom=True)
    for i, dd in pct_labels.iterrows():
        ax.text(20, i, f"{np.round(dd.pct_samples, 4)}%", va="center", ha="center")

    plt.savefig(snakemake.output[0])


def get_label_percentages():
    client = MongoClient()
    db = client["sramongo"]
    ncbi = db["ncbi"]

    df = (
        pd.DataFrame(
            ncbi.aggregate(
                [
                    {"$match": {"sra_create_date": {"$lte": date_parse("2019-03-01")}}},
                    {"$project": {"_id": False, "srx": True, "library_strategy": True}},
                ]
            )
        )
        .set_index("srx")
        .squeeze()
    )

    client.close()

    return (
        (df.value_counts() / df.shape[0] * 100)
        .rename_axis("Label")
        .rename("pct_samples")
        .to_frame()
        .reset_index()
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(workdir="paper-wf")

    main()
