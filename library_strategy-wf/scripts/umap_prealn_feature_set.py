import os

import joblib
import pandas as pd
import umap
from sklearn.preprocessing import StandardScaler


def main():
    features = pd.read_parquet(snakemake.input[0]).dropna()

    scaler = StandardScaler()
    X = scaler.fit_transform(features)

    reducer = umap.UMAP()
    embeddings = reducer.fit_transform(X)
    df = pd.DataFrame(embeddings, columns=["UMAP1", "UMAP2"], index=features.index)

    # save
    joblib.dump(reducer, snakemake.output.model)
    df.to_csv(snakemake.output.embeddings, sep="\t")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/prealn_feature_set.parquet"
        )

    main()
