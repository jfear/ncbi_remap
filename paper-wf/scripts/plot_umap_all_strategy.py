import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import sra_umap


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    plot = sra_umap.Plot(
        embeddings_file=snakemake.input.embeddings,
        labels_file=snakemake.input.labels,
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        ax_kwargs=dict(snakemake.params.get("ax_kwargs", {})),
    )

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        # %load_ext autoreload
        # %autoreload 2
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                embeddings="../output/library_strategy-wf/umap_prealn_features_embeddings.tsv",
                labels="../output/library_strategy-wf/sra_strategy_selection.parquet"
            ),
            params=dict(
                plot_kwargs=dict(),
                ax_kwargs=dict(),
            ),
        )

    main()
