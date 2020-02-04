import os

import joblib
import numpy as np
import pandas as pd
import umap

np.random.seed(42)

def main():
    scaled_features = pd.read_parquet(snakemake.input[0])
    X = scaled_features.values

    reducer = umap.UMAP()
    embeddings = reducer.fit_transform(X)
    df = pd.DataFrame(embeddings, columns=["UMAP1", "UMAP2"], index=scaled_features.index)

    # save
    joblib.dump(reducer, snakemake.output.model)
    df.to_csv(snakemake.output.embeddings, sep="\t")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/scaled_prealn_feature_set.parquet"
        )

    main()
