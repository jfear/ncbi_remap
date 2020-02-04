import os

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def main():
    features = pd.read_parquet(snakemake.input[0]).dropna()

    scaler = StandardScaler()
    X = scaler.fit_transform(features)

    reducer = PCA(n_components=12).fit(X)
    component_names = pd.Index([f"PC{x+1}" for x in range(reducer.n_components_)], name="PCs")

    embeddings = pd.DataFrame(reducer.transform(X), index=features.index, columns=component_names)

    loadings = pd.DataFrame(reducer.components_, index=component_names, columns=features.columns)

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
            input="../../output/library_strategy-wf/prealn_feature_set.parquet"
        )

    main()
