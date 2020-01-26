import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import prealn_features


def main():
    plt.style.use(snakemake.params.get("style", "sra"))
    plt.rcParams["axes.labelsize"] = 15

    plot = prealn_features.Plot(
        snakemake.input[0],
        grid_kwargs=dict(snakemake.params.get("grid_kwargs", {})),
        diag_plot_kwargs=dict(snakemake.params.get("diag_plot_kwargs", {})),
        offdiag_plot_kwargs=dict(snakemake.params.get("offdiag_plot_kwargs", {})),
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
