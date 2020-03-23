import os
import sys

import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, "../src")
from ncbi_remap.plotting.umap_panel import Plot


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    umap = pd.read_parquet(snakemake.input.umap)
    labels = pd.read_parquet(snakemake.input.labels)

    Plot(
        embeddings=umap,
        labels=labels.library_strategy,
        panel_kwargs=snakemake.params.get("panel_kwargs", {}),
        plot_kwargs=snakemake.params.get("plot_kwargs", {}),
    ).panel.savefig(snakemake.output.strategy)

    Plot(
        embeddings=umap,
        labels=labels.library_selection,
        panel_kwargs=snakemake.params.get("panel_kwargs", {}),
        plot_kwargs=snakemake.params.get("plot_kwargs", {}),
    ).panel.savefig(snakemake.output.selection)


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
                umap="../../output/library_strategy-wf/umap_prealn_features_embeddings.parquet",
            ),
            output=dict(strategy="", selection=""),
        )

    main()
