"""Plot all SHAP scores."""
import sys

import joblib
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, "../src")
from ncbi_remap.iforest import SraIsolationForest


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
    "gene_body_middle": "Mid Gene Body Coverage (avg)",
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
    df = get_features_w_shap(iso)
    df_w_flags = add_outliers(df, iso)

    g = sns.FacetGrid(
        df_w_flags,
        col="feature",
        col_wrap=4,
        hue="RNA-Seq Outlier",
        hue_order=[False, True],
        palette=["C0", "C3"],
        sharex=False,
        sharey=False,
    )
    g.map(sns.scatterplot, "feature_score", "shap_value", s=20, rasterized=True)
    g.set_titles("{col_name}")

    g.savefig(snakemake.output[0])


def get_features_w_shap(iso):
    features = iso.X_test
    shap_values = pd.DataFrame(iso.shap_values, index=features.index, columns=features.columns)

    # melt for seaborn compatability
    features_melt = features.reset_index().melt(
        id_vars="srx", var_name="feature", value_name="feature_score"
    )
    shap_values_melt = shap_values.reset_index().melt(
        id_vars="srx", var_name="feature", value_name="shap_value"
    )

    return features_melt.merge(shap_values_melt, on=["srx", "feature"]).replace(NAME_MAPPER)


def add_outliers(df, iso):
    outlier_srx = iso.outliers_test.index
    df["RNA-Seq Outlier"] = False
    df.loc[df.srx.isin(outlier_srx), "RNA-Seq Outlier"] = True
    return df


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(input="../../output/library_strategy-wf/isolation_forest.pkl",)

    main()
