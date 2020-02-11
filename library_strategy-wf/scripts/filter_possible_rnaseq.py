import os
from pickle import dump

import pandas as pd

OK_TERMS = [
    "3Prime-Seq",
    "CAP-Seq",
    "CAP-Seq|PRO-Seq|RNA-Seq|OTHER",
    "ChIP-Seq|RNA-Seq",
    "EST|RNA-Seq",
    "miRNA-Seq",
    "Nascent-Seq",
    "ncRNA-Seq",
    "ncRNA-Seq|EST",
    "ncRNA-Seq|OTHER",
    "PRO-Seq",
    "RNA-Seq",
    "RNA-Seq|ChIP-Seq",
    "RNA-Seq|EST",
    "RNA-Seq|FL-cDNA",
    "RNA-Seq|miRNA-Seq",
    "RNA-Seq|ncRNA-Seq|miRNA-Seq",
    "RNA-Seq|OTHER",
    "Start-Seq",
]


def main():
    library_strategy = pd.read_parquet(snakemake.input[0]).iloc[:, 0].squeeze()

    srxs = set([k for k, v in library_strategy.items() if v in OK_TERMS])

    dump(srxs, open(snakemake.output[0], "wb"))


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/summarized_metadata.parquet"
        )

    main()
