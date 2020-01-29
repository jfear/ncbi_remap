import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import umap


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    plot = umap.Plot(
        features_file=snakemake.input.features,
        labels_file=snakemake.input.sra_labels,
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        ax_kwargs=dict(snakemake.params.get("ax_kwargs", {})),
    )

    if "title" in snakemake.params:
        plt.suptitle(snakemake.params.title

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        # %load_ext autoreload
        # %autoreload 2
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                features="../../output/library_strategy-wf/prealn_feature_set.parquet",
                sra_labels="../../output/library_strategy-wf/sra_strategy_selection.parquet,
            )
            params=dict(
                plot_kwargs=dict(),
                ax_kwargs=dict(),
            ),
        )

    main()
