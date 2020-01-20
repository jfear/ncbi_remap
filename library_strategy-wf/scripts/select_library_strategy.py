"""Take all the different library strategies and select the best one."""
import os
import pandas as pd


def main():
    libstrat = load_data()

    # remove samples with no data, these have not completed the workflow
    libstrat = libstrat[~libstrat.random_forest.isnull()]

    # Re-annotate samples
    df = libstrat.apply(annot_library_strategy, axis=1)
    df.name = "Fear_et_al_library_strategy"
    df.index.name = "srx"
    df.to_frame().to_parquet(snakemake.output[0])


def load_data():
    # Grab the different versions of library strategy annotation
    sra = pd.read_parquet(snakemake.input["sra"]).library_strategy.rename("sra")
    free_text = pd.read_parquet(snakemake.input["free_text"]).squeeze().rename("author")
    random_forest = (
        pd.read_parquet(snakemake.input["forest"])
        .reset_index()
        .melt(id_vars="srx")
        .groupby("srx")
        .value.value_counts()
        .unstack()
        .idxmax(axis=1)
        .rename("random_forest")
    )

    random_forest_other = (
        pd.read_parquet(snakemake.input["other"])
        .reset_index()
        .melt(id_vars="srx")
        .groupby("srx")
        .value.value_counts()
        .unstack()
        .idxmax(axis=1)
        .rename("random_forest")
    )

    random_forest_combined = pd.concat([random_forest, random_forest_other])

    return pd.concat([sra, free_text, random_forest_combined], sort=True, axis=1)


def annot_library_strategy(x):
    """Selection function"""
    # If all same
    if (x.sra == x.author) & (x.sra == x.random_forest):
        return x.sra

    # If there is author metadata
    if x.author:
        # Re-annotate OTHER samples if author and random forest match
        if (x.sra == "OTHER") & (x.author == x.random_forest):
            return x.author

        # Re-annotate OTHER samples if author has a single value
        if (x.sra == "OTHER") & ~("|" in x.author):
            return x.author

        # Ambiguous so output everything, but keep author upfront.
        if x.sra != x.random_forest:
            string = x.author

            if x.random_forest not in string:
                string = string + "|" + x.random_forest

            if x.sra not in string:
                string = string + "|" + x.sra

            return string

    # If there is no author metadata
    # If SRA/Machine same
    if x.sra == x.random_forest:
        return x.sra

    # EST and RNA-Seq are really similar, so if the SRA uses EST then stick
    # with that.
    if (x.sra == "EST") & (x.random_forest == "RNA-Seq"):
        return "EST"

    # WGA and WGS are really similar, so if the SRA uses WGA then stick
    # with that.
    if (x.sra == "WGA") & (x.random_forest == "WGS"):
        return "WGA"

    # Ambiguous so output everything, but keep machine up front
    return "|".join([x.random_forest, x.sra])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                sra="../../output/library_strategy-wf/munge_library_strategy_from_mongo.parquet",
                free_text="../../output/library_strategy-wf/free_text_library_strategy.parquet",
                forest="../../output/library_strategy-wf/random_forest_library_strategy.parquet",
                other="../../output/library_strategy-wf/random_forest_library_strategy_other.parquet",
            )
        )

    main()
