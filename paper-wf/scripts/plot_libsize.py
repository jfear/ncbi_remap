import os

import matplotlib.pyplot as plt
import seaborn as sns

from ncbi_remap.plotting import library_size


def main():
    plt.style.use(snakemake.params.get("style", "sra"))

    plot = library_size.Plot(
        snakemake.input[0],
        plot_kwargs=dict(snakemake.params.get("plot_kwargs", {})),
        joint_ax_kwargs=dict(snakemake.params.get("joint_ax_kwargs", {})),
        libsize_ax_kwargs=dict(snakemake.params.get("libsize_ax_kwargs", {})),
        duplication_ax_kwargs=dict(snakemake.params.get("duplication_ax_kwargs", {})),
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
                title="Library Size",
                style=("sra", "sra_talk"),
                plot_kwargs=dict(),
                joint_ax_kwargs=dict(),
                libsize_ax_kwargs=dict(),
                duplication_ax_kwargs=dict(),
            ),
        )

    main()
