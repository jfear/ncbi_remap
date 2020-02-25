import os

import pandas as pd
from joblib import dump
from sklearn.covariance import ShrunkCovariance


def main():
    reduced = pd.read_parquet(snakemake.input[0])
    cov = ShrunkCovariance(assume_centered=True)
    cov.fit(reduced)

    dump(cov, snakemake.output.model)
    dump(cov.covariance_, snakemake.output.cov)
    dump(cov.precision_, snakemake.output.inv_cov)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/agg-rnaseq-wf/sample_pca_select_components.parquet",
            output=dict(
                model="../../output/agg-rnaseq-wf/sample_cov_model.pkl",
                cov="../../output/agg-rnaseq-wf/sample_cov.pkl",
                inv_cov="../../output/agg-rnaseq-wf/sample_inverse_cov.pkl",
            ),
        )

    main()
