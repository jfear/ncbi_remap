"""Pull out RNA-Seq samples from feature set"""
import pandas as pd


def main():
    rnaseq_srxs = (
        pd.read_parquet(snakemake.input.labels).query("library_strategy == 'RNA-Seq'").index
    )

    (
        pd.read_parquet(snakemake.input.features)
        .reindex(rnaseq_srxs)
        .dropna(how="all")
        .sort_index()
        .to_parquet(snakemake.output[0])
    )


if __name__ == "__main__":
    if "snakemake" not in locals() or not hasattr(snakemake, "scriptdir"):
        from ncbi_remap.mock import MockSnake

        snakemake = MockSnake(
            input=dict(
                features="../../output/library_strategy-wf/prealn_features.parquet",
                labels="../../output/library_strategy-wf/sra_strategy_selection.parquet",
            ),
            output="",
        )

    main()
