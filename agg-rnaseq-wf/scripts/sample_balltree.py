import os

import pandas as pd
from sklearn.neighbors import BallTree
from joblib import load, dump


def main():
    df = pd.read_parquet(snakemake.input.selected_values)
    vi = load(snakemake.input.vi)
    searcher = BallTree(df, metric="mahalanobis", VI=vi)

    dump(searcher, snakemake.output[0])

if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                selected_values="../../output/agg-rnaseq-wf/feature_pca_select_components.parquet",
                vi="../../output/agg-rnaseq-wf/feature_inverse_cov.pkl",
            )
        )

    main()
