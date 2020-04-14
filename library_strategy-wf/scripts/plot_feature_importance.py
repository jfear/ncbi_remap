"""Feature Importance"""
import sys
import joblib
import numpy as np
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

    iso = joblib.load(snakemake.input[0])  # type: SraIsolationForest
    mean_shap, columns = iso.mean_shap_values()
    cols = [NAME_MAPPER[x] for x in columns]

    fig, ax = plt.subplots(figsize=plt.figaspect(1.2))
    sns.barplot(mean_shap, cols, color="C0", ax=ax)
    ax.set(xlabel="mean(|SHAP Value|)", ylabel="")

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input="../../output/library_strategy-wf/isolation_forest.pkl", output="",
        )

    main()
