import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

np.random.seed(42)

def main():
    scaled_features = pd.read_parquet(snakemake.input[0])
    X = scaled_features.values

    reducer = PCA(n_components=12).fit(X)
    component_names = pd.Index([f"PC{x+1}" for x in range(reducer.n_components_)], name="PCs")

    embeddings = pd.DataFrame(reducer.transform(X), index=scaled_features.index, columns=component_names)

    loadings = pd.DataFrame(reducer.components_, index=component_names, columns=scaled_features.columns)

    explained_variance = pd.DataFrame(
        reducer.explained_variance_ratio_, index=component_names, columns=["explained_variance"]
    )

    embeddings.to_csv(snakemake.output.embeddings, sep="\t")
    loadings.to_csv(snakemake.output.loadings, sep="\t")
    explained_variance.to_csv(snakemake.output.variance, sep="\t")


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/scaled_prealn_feature_set.parquet"
        )

    main()
