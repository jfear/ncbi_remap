import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import library_layout


def main():
    plt.style.use(snakemake.params.style)

    plot = library_layout.Plot(
        snakemake.input[0],
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        annot_kwargs=dict(snakemake.params.get("annot_kwargs", {})),
        ax_kwargs=dict(snakemake.params.get("ax_kwargs", {})),
    )
    plt.suptitle(snakemake.params.title)

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
                title="Library Layout",
                plot_kwargs=dict(color="white", edgecolor="k"),
                ax_kwargs=dict(ylabel="Samples"),
            ),
        )

    main()
