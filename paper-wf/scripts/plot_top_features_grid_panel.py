import os

import matplotlib.pyplot as plt

from ncbi_remap.plotting.top_features_panel import Plot


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    Plot(
        features_file=snakemake.input.features,
        importance_file=snakemake.input.importance,
        labels_file=snakemake.input.labels,
        panel_kwargs=snakemake.params.get("panel_kwargs", None),
        plot_kwargs=snakemake.params.get("plot_kwargs", None),
    ).savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                features="../../output/library_strategy-wf/scaled_prealn_feature_set.parquet",
                importance="../../output/library_strategy-wf/random_forest_feature_importance.tsv",
                labels="../../output/library_strategy-wf/summarized_metadata.parquet"
            )
        )

    main()
