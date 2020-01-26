import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import feature_importance


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    plot = feature_importance.Plot(
        snakemake.input[0],
        strategy_plot_kwargs=dict(snakemake.params.get("strategy_plot_kwargs", {})),
        selection_plot_kwargs=dict(snakemake.params.get("selection_plot_kwargs", {})),
        strategy_axes_kwargs=dict(snakemake.params.get("strategy_axes_kwargs", {})),
        selection_axes_kwargs=dict(snakemake.params.get("selection_axes_kwargs", {})),
    )
    plt.suptitle(snakemake.params.get("title", ""))

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        # %load_ext autoreload
        # %autoreload 2
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/paper-wf/srx_prealn.tsv",
            params=dict(
                style=("sra", "sra_talk"),
            ),
        )

    main()
