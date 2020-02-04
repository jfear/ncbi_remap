import os

import pandas as pd
from sklearn.preprocessing import StandardScaler


def main():
    df = pd.read_parquet(snakemake.input[0]).dropna()
    scaler = StandardScaler()

    pd.DataFrame(scaler.fit_transform(df), index=df.index, columns=df.columns).to_parquet(
        snakemake.output[0]
    )


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/prealn_feature_set.parquet"
        )

    main()
