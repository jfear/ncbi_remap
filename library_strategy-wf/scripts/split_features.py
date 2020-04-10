"""Pull out RNA-Seq samples from feature set"""
import pandas as pd

MAPPER = {"rnaseq": "RNA-Seq", "est": "EST", "wgs": "WGS", "chip": "ChIP-Seq"}


def main():
    labels = pd.read_parquet(snakemake.input.labels).library_strategy.squeeze()
    features = pd.read_parquet(snakemake.input.features)

    for short_name, long_name in MAPPER.items():
        srxs = labels[labels == long_name].index
        (
            features.reindex(srxs)
            .dropna(how="all")
            .sort_index()
            .to_parquet(snakemake.output[short_name])
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
