"""Plot Feature Interaction Panel"""
import sys

import joblib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, SubplotSpec
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use
from ncbi_remap.iforest import SraIsolationForest

COLORS = ["C0", "C3"]
SELECTED_FEATURES = [
    ("percent_rrna_reads", "percent_mrna_bases"),
    ("percent_alignment", "gene_body_middle"),
    ("number_junction_reads", "number_multimapping_reads"),
    ("percent_intergenic_bases", "percent_duplication"),
    ("number_genic_reads", "gene_body_five_prime"),
    ("percent_intergenic_bases", "percent_intronic_bases"),
]
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

style_use(snakemake.params.get("style", "sra"))
plt.rcParams["figure.figsize"] = (6, 8)

SCATTER_PARAMS = dict(
    # s=plt.rcParams["lines.markersize"] + 2,
    edgecolors=plt.rcParams["scatter.edgecolors"]
)

DIST_PARAMS = dict(
    hist=False,
    kde_kws=dict(shade=True)
)


def main():
    iso = joblib.load(snakemake.input[0])

    _joint = 1
    _marg = 0.2
    _space = 0.001

    fig = plt.figure()
    gs = GridSpec(
        nrows=8,
        ncols=5,
        figure=fig,
        width_ratios=[_joint, _marg, _space, _joint, _marg],
        height_ratios=[_marg, _joint, _space, _marg, _joint, _space, _marg, _joint],
        hspace=0,
        wspace=0
    )
    plot_panel(
        gs[1, 0], gs[0, 0], gs[1, 1], fig, SELECTED_FEATURES[0][0], SELECTED_FEATURES[0][1], iso
    )
    plot_panel(
        gs[1, 3], gs[0, 3], gs[1, 4], fig, SELECTED_FEATURES[1][0], SELECTED_FEATURES[1][1], iso
    )

    plot_panel(
        gs[4, 0], gs[3, 0], gs[4, 1], fig, SELECTED_FEATURES[2][0], SELECTED_FEATURES[2][1], iso
    )
    plot_panel(
        gs[4, 3], gs[3, 3], gs[4, 4], fig, SELECTED_FEATURES[3][0], SELECTED_FEATURES[3][1], iso
    )

    plot_panel(
        gs[7, 0], gs[6, 0], gs[7, 1], fig, SELECTED_FEATURES[4][0], SELECTED_FEATURES[4][1], iso
    )
    plot_panel(
        gs[7, 3], gs[6, 3], gs[7, 4], fig, SELECTED_FEATURES[5][0], SELECTED_FEATURES[5][1], iso
    )

    plt.savefig(snakemake.output[0])


def plot_panel(
    joint: SubplotSpec,
    x_marg: SubplotSpec,
    y_marg: SubplotSpec,
    fig: plt.Figure,
    feature1: str,
    feature2: str,
    iso: SraIsolationForest,
):
    ax_joint = fig.add_subplot(joint)
    ax_x_marg = fig.add_subplot(x_marg)
    ax_y_marg = fig.add_subplot(y_marg)

    X = iso.X_test[feature1]
    Y = iso.X_test[feature2]

    inlier = iso.isinlier_test
    outlier = iso.isoutlier_test

    ax_joint.scatter(X[inlier], Y[inlier], zorder=0, rasterized=True, color=COLORS[0], **SCATTER_PARAMS)
    ax_joint.scatter(X[outlier], Y[outlier], zorder=1, rasterized=True, color=COLORS[1], **SCATTER_PARAMS)
    ax_joint.set(xlabel=NAME_MAPPER[feature1], ylabel=NAME_MAPPER[feature2])
    ax_joint.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax_joint.yaxis.set_major_locator(plt.MaxNLocator(3))

    sns.distplot(X[inlier], ax=ax_x_marg, color=COLORS[0], **DIST_PARAMS)
    sns.distplot(X[outlier], ax=ax_x_marg, color=COLORS[1], **DIST_PARAMS)
    ax_x_marg.set(xlabel="", ylabel="", xticks=[], yticks=[])
    sns.despine(ax=ax_x_marg, left=True, bottom=True)

    sns.distplot(Y[inlier], vertical=True, ax=ax_y_marg, color=COLORS[0], **DIST_PARAMS)
    sns.distplot(Y[outlier], vertical=True, ax=ax_y_marg, color=COLORS[1], **DIST_PARAMS)
    ax_y_marg.set(xlabel="", ylabel="", xticks=[], yticks=[])
    sns.despine(ax=ax_y_marg, left=True, bottom=True)


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(input="../../output/library_strategy-wf/isolation_forest.pkl",)

    main()
