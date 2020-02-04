import os

import pandas as pd
import matplotlib.pyplot as plt

from ncbi_remap.plotting.umap_panel import Plot


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    umap = pd.read_parquet(snakemake.input.umap)
    outliers = pd.read_parquet(snakemake.input.outliers)

    Plot(
        embeddings=umap,
        labels=outliers.library_strategy,
        outlier_flag=outliers.library_strategy_flag_outlier,
        panel_kwargs=snakemake.params.get("panel_kwargs", {}),
        plot_kwargs=snakemake.params.get("plot_kwargs", {}),
    ).panel.savefig(snakemake.output.strategy)

    Plot(
        embeddings=umap,
        labels=outliers.library_selection,
        outlier_flag=outliers.library_selection_flag_outlier,
        panel_kwargs=snakemake.params.get("panel_kwargs", {}),
        plot_kwargs=snakemake.params.get("plot_kwargs", {}),
    ).panel.savefig(snakemake.output.selection)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                outliers="../../output/library_strategy-wf/prealn_outliers.parquet",
                umap="../../output/library_strategy-wf/umap_prealn_features_embeddings.parquet",
            ),
            output=dict(strategy="", selection=""),
        )

    main()
