"""Plot all features pairwise."""
import sys
from pathlib import Path
import tarfile
from itertools import permutations

import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

NAME_MAPPER = {
    "percent_alignment": r"Reads Aligned (%)",
    "percent_duplication": f"Duplicated Reads (%)",
    "percent_genes_on": f"Gene Expressed (%)",
    "number_junction_reads": r"Reads at Exon Junctions (#)",
    "percent_utr_bases": r"UTR Bases (%)",
    "gene_body_three_prime": "3' Gene Body Coverage (avg)",
    "percent_intergenic_bases": r"Intergenic Bases (%)",
    "percent_reverse": r"Reverse Oriented Reads (%)",
    "number_genic_reads": "Genic Reads (#)",
    "gene_body_middle": "Middle Gene Body Coverage (avg)",
    "number_reads": "Total Reads (#)",
    "median_cv_coverage": "Coef Var Coverage (median)",
    "gene_body_five_prime": "5' Gene Body Coverage (avg)",
    "percent_mrna_bases": r"mRNA Bases (%)",
    "reads_MQ0": "Reads 0 Mapping Quality (#)",
    "percent_rrna_reads": r"rRNA Reads (%)",
    "number_junctions_on": "Exon Junctions Touched (#)",
    "percent_intronic_bases": r"Intronic Bases (%)",
    "number_multimapping_reads": "Multimapping Reads (#)",
    "average_quality": "Quality Score (avg)",
    "number_reads_too_short": "Reads Too Short (#)",
}


def main():
    style_use(snakemake.params.get("style", "sra"))
    plt.rcParams["figure.constrained_layout.use"] = False

    iso = joblib.load(snakemake.input[0])  # type : ncbi_remap.iforest.SraIsolationForest
    features = iso.features.copy()
    features["rnaseq_outliers"] = iso.isoutlier_all

    # Save each plot seperately
    folder = prep_folder()
    for f1, f2 in permutations(features.columns[:-1], 2):
        joint_plot(f1, f2, features)
        plt.savefig(folder / f"{f1}_vs_{f2}.svg")
        plt.close()

    tar_folder(folder)


def prep_folder():
    folder = Path(snakemake.output[0]).with_suffix("")
    folder.mkdir(exist_ok=True)
    return folder


def joint_plot(x, y, features):
    sns_params = dict(
        hue="rnaseq_outliers", hue_order=[False, True], palette=["C0", "C3"], rasterized=True
    )

    inliers = features[~features.rnaseq_outliers]
    outliers = features[features.rnaseq_outliers]

    g = sns.JointGrid(x, y, data=features)
    sns.scatterplot(x, y, data=features, ax=g.ax_joint, **sns_params)
    kdeplot(inliers[x], color="C0", shade=True, ax=g.ax_marg_x)
    kdeplot(outliers[x], color="C3", shade=True, ax=g.ax_marg_x)
    kdeplot(inliers[y], color="C0", shade=True, vertical=True, ax=g.ax_marg_y)
    kdeplot(outliers[y], color="C3", shade=True, vertical=True, ax=g.ax_marg_y)
    g.ax_joint.legend_.remove()
    g.ax_joint.margins(0.01)
    g.ax_joint.set(xlabel=NAME_MAPPER[x], ylabel=NAME_MAPPER[y])

    return g


def kdeplot(*args, **kwargs):
    try:
        sns.kdeplot(*args, **kwargs)
        kwargs["ax"].legend_.remove()
    except RuntimeError:
        kwargs["kde"] = False
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
            input="../../output/library_strategy-wf/isolation_forest.pkl",
            output="../../output/library_strategy-wf/pairlot_features.tgz",
        )

    main()
