"""Plot Shap Interactions"""
import sys
from itertools import permutations
from pathlib import Path
import tarfile

import joblib
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

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

SCATTER_DEFAULTS = dict(cmap="plasma", rasterized=True)


def main():
    # style_use()
    plt.rcParams["figure.dpi"] = 100

    iso = joblib.load(snakemake.input[0])
    folder = prep_folder()

    for f1, f2 in permutations(iso.columns, 2):
        plot(f1, f2, iso)
        plt.savefig(folder / f"{f1}_vs_{f2}.svg")
        plt.close()

    for feature in iso.columns:
        plot_no_color(feature, iso)
        plt.savefig(folder / f"{feature}.svg")
        plt.close()

    tar_folder(folder)


def plot(feature1, feature2, iso):
    # Get data
    idx1 = iso.columns.get_loc(feature1)
    idx2 = iso.columns.get_loc(feature2)
    x = iso.X_test.iloc[:, idx1]
    y = iso.shap_values[:, idx1]
    color = iso.X_test.iloc[:, idx2]
    # Set color scale
    vmin, vmax = color_scale(color)
    norm = mpl.colors.Normalize(vmin, vmax)
    sm = mpl.cm.ScalarMappable(norm, cmap=SCATTER_DEFAULTS["cmap"])

    # Plot
    _, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=plt.figaspect(1 / 2), gridspec_kw=dict(width_ratios=[1, 1, 0.05]),
    )

    ax1.scatter(
        x[iso.isinlier_test],
        y[iso.isinlier_test],
        c=color[iso.isinlier_test],
        vmin=vmin,
        vmax=vmax,
        **SCATTER_DEFAULTS,
    )

    ax2.scatter(
        x[iso.isoutlier_test],
        y[iso.isoutlier_test],
        c=color[iso.isoutlier_test],
        vmin=vmin,
        vmax=vmax,
        **SCATTER_DEFAULTS,
    )

    # Add colorbar and tweak plots
    xlim = get_limits(x)
    ylim = get_limits(y)
    params = dict(xlim=xlim, ylim=ylim, xlabel=NAME_MAPPER[feature1])

    plt.colorbar(sm, label=NAME_MAPPER[feature2], cax=ax3)
    ax1.set(ylabel="SHAP Value", title="Inliers", **params)
    ax2.set(ylabel="", yticklabels=[], title="Outliers", **params)


def plot_no_color(feature1, iso):
    # Get data
    idx1 = iso.columns.get_loc(feature1)
    x = iso.X_test.iloc[:, idx1]
    y = iso.shap_values[:, idx1]

    # Plot
    _, ax = plt.subplots() # type: plt.Figure, plt.Axes
    ax.scatter(
        x[iso.isinlier_test], y[iso.isinlier_test], c="C0", edgecolors="w", **SCATTER_DEFAULTS,
    )

    ax.scatter(
        x[iso.isoutlier_test], y[iso.isoutlier_test], c="C3", edgecolors="w", **SCATTER_DEFAULTS,
    )

    ax.set(xlabel=NAME_MAPPER[feature1], ylabel="SHAP Value")
    sns.despine(ax=ax)


def color_scale(color):
    if 0 <= min(color) and max(color) <= 100:
        return 0, 100
    return np.percentile(color, [25, 75])


def prep_folder():
    folder = Path(snakemake.output[0]).with_suffix("")
    folder.mkdir(exist_ok=True)
    return folder


def tar_folder(folder):
    with tarfile.open(snakemake.output[0], mode="w:gz") as tar:
        for file_name in folder.iterdir():
            tar.add(file_name)


def get_limits(x):
    pad = 0.05 * (max(x) - min(x))
    return min(x) - pad, max(x) + pad


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(input="../../output/library_strategy-wf/isolation_forest.pkl",)

    main()
