import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import genes_on


def main():
    plt.style.use(snakemake.params.style)

    plot = genes_on.Plot(
        snakemake.input[0],
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        ax_kwargs=dict(snakemake.params.get("ax_kwargs", {})),
    )
    plt.suptitle(snakemake.params.get("title", ""))

    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        # %load_ext autoreload
        # %autoreload 2
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/paper-wf/aln_genic_counts.tsv",
            params=dict(
                title="Gene Diversity",
                style=("sra", "sra_talk"),
                plot_kwargs=dict(),
                ax_kwargs=dict(),
            ),
        )

    main()
