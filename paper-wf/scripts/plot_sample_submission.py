import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import sample_submission


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    plot = sample_submission.Plot(
        snakemake.input[0],
        bar_plot_kwargs=dict(snakemake.params.get("bar_plot_kwargs", {})),
        reg_plot_kwargs=dict(snakemake.params.get("reg_plot_kwargs", {})),
        bar_ax_kwargs=dict(snakemake.params.get("bar_ax_kwargs", {})),
        reg_ax_kwargs=dict(snakemake.params.get("reg_ax_kwargs", {})),
    )

    if "title" in snakemake.params:
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
                bar_plot_kwargs=dict(lowess=True),
                reg_plot_kwargs=dict(lowess=True),
                bar_ax_kwargs=dict(),
                reg_ax_kwargs=dict(),
            ),
        )

    main()
