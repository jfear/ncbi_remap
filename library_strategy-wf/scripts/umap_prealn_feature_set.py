import os
import pickle

import pandas as pd
from sklearn.preprocessing import StandardScaler
import umap


def main():
    features = pd.read_parquet(snakemake.input[0])

    scaler = StandardScaler()
    X = scaler.fit_transform(features)

    reducer = umap.UMAP()
    embedings = reducer.fit_transform(X)

    # save
    pickle.dump(reducer, open(snakemake.output.model, 'wb'))


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/prealn_feature_set.parquet"
        )

    main()
