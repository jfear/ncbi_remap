import os

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

import umap

np.random.seed(42)

def main():
    tpm = pd.read_csv(snakemake.input[0], sep="\t", index_col=0)
    scaler = StandardScaler() 
    X = scaler.fit_transform(tpm)

    reducer = umap.UMAP()
    embeddings = reducer.fit_transform(X)

    df = pd.DataFrame(embeddings, columns=["UMAP1", "UMAP2"], index=tpm.index)
    df.to_parquet(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/agg-rnaseq-wf/tpm_gene_counts.tsv"
        )

    main()
