import os

import numpy as np
import pandas as pd
import umap

RANDOM_STATE = np.random.RandomState(42)


def main():
    scaled_features = pd.read_parquet(snakemake.input[0])
    X = scaled_features.values

    reducer = umap.UMAP(random_state=RANDOM_STATE)
    embeddings = reducer.fit_transform(X)

    (
        pd.DataFrame(
            embeddings, columns=["UMAP1", "UMAP2"], index=scaled_features.index
        ).to_parquet(snakemake.output[0])
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/scaled_prealn_feature_set.parquet"
        )

    main()
