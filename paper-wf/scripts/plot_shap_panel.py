"""Plot Shap dependencies as a panel"""
import sys

import joblib
import matplotlib.pyplot as plt

sys.path.insert(0, "../src")
from ncbi_remap.plotting import style_use

COLORS = "C0", "C3"  # Inliers, Outliers
SELCTED_FEATURES = [
    "percent_rrna_reads",
    "percent_mrna_bases",
    "percent_alignment",
    "gene_body_middle",
    "number_junction_reads",
    "percent_intergenic_bases",
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


def main():
    iso = joblib.load(snakemake.input[0])  # type: SraIsolationForest

    style_use(snakemake.params.get("style", "sra"))
    plt.rcParams["figure.figsize"] = (6, 8)

    fig, axes = plt.subplots(3, 2)  # type: plt.Figure, plt.Axes
    for feature, ax in zip(SELCTED_FEATURES, axes.flat):
        plot_panel(iso, feature, ax)
    fig.text(0, 0.5, "SHAP Value", rotation=90, ha="right", va="center")
    fig.text(0.5, 0, "Feature Value", ha="center", va="top")

    plt.savefig(snakemake.output[0])


def plot_panel(iso, feature: str, ax: plt.Axes):
    idx = iso.columns.get_loc(feature)
    ax.scatter(
        iso.inliers_test[feature],
        iso.shap_values[iso.isinlier_test, idx],
        c=COLORS[0],
        zorder=0,
        rasterized=True,
    )
    ax.scatter(
        iso.outliers_test[feature],
        iso.shap_values[iso.isoutlier_test, idx],
        c=COLORS[1],
        zorder=1,
        rasterized=True,
    )
    ax.set_title(NAME_MAPPER[feature])
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    return ax


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(input="../../output/library_strategy-wf/isolation_forest.pkl",)

    main()
