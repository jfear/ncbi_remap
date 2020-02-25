import os

import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from joblib import dump


def main():
    scaler = StandardScaler()
    reducer = PCA(n_components=8_000, random_state=42)

    df = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    samples = df.index.values
    X = df.values
    del df

    scaled = scaler.fit_transform(X)
    del X

    reduced = reducer.fit_transform(scaled)
    del scaled

    pd.DataFrame(
        reduced, index=samples, columns=[f"PC{x+1}" for x in range(reduced.shape[1])]
    ).to_parquet(snakemake.output.reduced)
    dump(scaler, snakemake.output.scaler)
    dump(reducer, snakemake.output.reducer)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(input="../../output/agg-rnaseq-wf/tpm_gene_counts.tsv")

    main()
