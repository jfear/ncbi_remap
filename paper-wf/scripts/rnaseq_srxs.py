import os
import pandas as pd


def main():
    srxs = (
        pd.read_parquet(snakemake.input[0])
        .query("Fear_et_al_library_strategy == 'RNA-Seq'")
        .index
    )

    with open(snakemake.output[0], "w") as file_out:
        file_out.write("\n".join(srxs))


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG", False):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/summarized_metadata.parquet"
        )

    main()
