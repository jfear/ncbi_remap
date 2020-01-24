import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import sample_submission


def main():
    plt.style.use(snakemake.params.style)

    plot = sample_submission.Plot(
        snakemake.input[0],
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
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
                title="Sample Submission",
                style=("sra", "sra_talk"),
                plot_kwargs=dict(lowess=True),
                ax_kwargs=dict(),
            ),
        )

    main()
