"""Are outliers enriched for a specific Library Selection method?"""
import sys

import joblib
import pandas as pd

sys.path.insert(0, "../src")
from ncbi_remap.stats import run_chisq
from ncbi_remap.iforest import SraIsolationForest


def main():
    # Data Wrangle
    iso = joblib.load(snakemake.input.iso)  # type: SraIsolationForest
    df = (
        pd.read_parquet(snakemake.input.labels)
        .reindex(iso.X_test.index)
        .assign(is_outlier=lambda x: iso.isoutlier_test)
        .drop(columns="library_strategy")
    )
    ct = pd.crosstab(df.library_selection, df.is_outlier)

    # Run chi^2 with post hoc tests
    chi2_res = run_chisq(ct.T)

    # Save results
    joblib.dump(chi2_res, snakemake.output[0])


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                iso="../../output/library_strategy-wf/isolation_forest.pkl",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            ),
        )

    main()
