import os

import pandas as pd
import matplotlib.pyplot as plt

from ncbi_remap.plotting import isolation_score

def main():
    plt.style.use(snakemake.params.get("style", "sra"))
    isolation_score.Plot(snakemake.input[0])
    plt.savefig(snakemake.output[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/prealn_outliers.parquet"
        )

    main()
