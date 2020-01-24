import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import read_length


def main():
    plt.style.use(snakemake.params.style)

    plot = read_length.Plot(
        snakemake.input[0],
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        ax_kwargs=dict(snakemake.params.get("ax_kwargs", {})),
    )
    plt.setp(plot.ax.get_xticklabels(), rotation=45, ha="right")
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
                title="Read Length",
                plot_kwargs=dict(
                    showfliers=False,
                    color="white",
                    boxprops=dict(edgecolor="k"),
                    medianprops=dict(color="k"),
                    whiskerprops=dict(color="k"),
                    capprops=dict(color="k"),
                ),
                ax_kwargs=dict(ylabel="Avg. Read Length (bp)", xlabel="Year"),
            ),
        )

    main()
