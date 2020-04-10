"""Plot all features pairwise."""
from pathlib import Path
import tarfile
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

SCATTER_STYLE = dict(s=10, edgecolors="w", linewidths=0.2, rasterized=True)


def main():
    style_use(snakemake.params.get("style", "sra"))

    outliers = pd.read_parquet(snakemake.input.outliers).squeeze()
    features = (
        pd.read_parquet(snakemake.input.features)
        .join(outliers, how="inner")
        .sample(n=1_000)
    )

    # Save each plot seperately
    folder = prep_folder()
    for f1, f2 in combinations(features.columns[:-1], 2):
        joint_plot(f1, f2, features)
        plt.savefig(folder / f"{f1}_vs_{f2}.svg")
        plt.close()

    tar_folder(folder)


def prep_folder():
    folder = Path(snakemake.output[0]).with_suffix("")
    folder.mkdir(exist_ok=True)
    return folder


def joint_plot(x, y, features):
    sns_params = dict(hue="rnaseq_outliers", hue_order=[False, True], palette=["C0", "C3"],)

    inliers = features[~features.rnaseq_outliers]
    outliers = features[features.rnaseq_outliers]

    g = sns.JointGrid(x, y, data=features)
    sns.scatterplot(x, y, data=features, ax=g.ax_joint, **sns_params)
    kdeplot(inliers[x], c="C0", shade=True, ax=g.ax_marg_x)
    kdeplot(outliers[x], c="C3", shade=True, ax=g.ax_marg_x)
    kdeplot(inliers[y], c="C0", shade=True, vertical=True, ax=g.ax_marg_y)
    kdeplot(outliers[y], c="C3", shade=True, vertical=True, ax=g.ax_marg_y)

    return g


def kdeplot(*args, **kwargs):
    try:
        sns.kdeplot(*args, **kwargs)
        kwargs["ax"].legend_.remove()
    except RuntimeError:
        kwargs["kde"] = False
        kwargs["color"] = kwargs["c"]
        del kwargs["c"]
        del kwargs["shade"]
        sns.distplot(*args, **kwargs)
        kwargs["ax"].set(xlabel="", ylabel="")


def tar_folder(folder):
    with tarfile.open(snakemake.output[0], mode="w:gz") as tar:
        for file_name in folder.iterdir():
            tar.add(file_name)


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                features="../../output/library_strategy-wf/prealn_features.parquet",
                outliers="../../output/library_strategy-wf/rnaseq_outliers.parquet",
            ),
            output="../../output/library_strategy-wf/pairlot_features.tgz",
        )

    main()
