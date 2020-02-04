import os

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import sra_umap


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    umap = pd.read_parquet(snakemake.input.embeddings)
    outliers = pd.read_parquet(snakemake.input.outliers)

    plot = sra_umap.Plot(
        embeddings=umap,
        labels=outliers.library_strategy,
        outlier_flag=outliers.library_strategy_flag_outlier,
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        ax_kwargs=dict(snakemake.params.get("ax_kwargs", {})),
    )

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                embeddings="../../output/library_strategy-wf/umap_prealn_features_embeddings.parquet",
                outliers="../../output/library_strategy-wf/prealn_outliers.parquet"
            ),
        )

    main()
