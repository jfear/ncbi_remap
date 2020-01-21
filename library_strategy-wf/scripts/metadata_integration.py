"""Aggregate library strategy and selection.

I have library strategies from the SRA, free text metadata, and a random
forest classifier. I have library selection only from the SRA and the
classifier.

"""
import os
import pandas as pd


def main():
    metadata = load_data()

    library_strategy = (
        metadata.apply(annot_library_strategy, axis=1)
        .rename("Fear_et_al_library_strategy")
        .rename_axis("srx")
    )

    library_selection = (
        metadata.apply(annot_library_selection, axis=1)
        .rename("Fear_et_al_library_selection")
        .rename_axis("srx")
    )

    pd.concat([library_strategy, library_selection], axis=1, sort=False).apply(
        cross_annotation, axis=1
    ).to_parquet(snakemake.output[0])


def load_data() -> pd.DataFrame:
    # Grab the different versions of library strategy annotation
    sra = (
        pd.read_parquet(snakemake.input["sra"])
        .rename(lambda x: x.replace("library", "sra"), axis=1)
        .replace({"other": "OTHER"})
    )

    free_text = pd.read_parquet(snakemake.input["free_text"]).squeeze().rename("author_strategy")

    random_forest = pd.read_parquet(snakemake.input["forest"]).rename(
        lambda x: x.replace("library_", ""), axis=1
    )

    return pd.concat([sra, free_text, random_forest], sort=True, axis=1, join="inner")


def annot_library_strategy(row: pd.Series) -> str:
    """Aggregates library strategies.
    
    Combines library strategies from the SRA, free text (author), and
    predicted values from random forest. 

    I only replace OTHER if I have support from both free text and random
    forest.

    If ambiguous I concatenate values in the following order:

    `author|predicted|sra`

    The SRA value is always last.

    """
    # If all same
    if row.sra_strategy == row.author_strategy == row.pred_strategy:
        return row.sra_strategy

    # If there is author metadata
    if row.author_strategy:
        # Re-annotate OTHER samples if author and random forest match
        if (row.sra_strategy == "OTHER") & (row.author_strategy == row.pred_strategy):
            return row.author_strategy

        # Re-annotate OTHER samples if author has a single value
        if (row.sra_strategy == "OTHER") & ~("|" in row.author_strategy):
            return row.author_strategy

        # Ambiguous: I am building the string like this because I want a unique
        # list and author may already have a few strategies.
        if row.sra_strategy != row.pred_strategy:
            string = row.author_strategy

            if row.pred_strategy not in string:
                string = string + "|" + row.pred_strategy

            if row.sra_strategy not in string:
                string = string + "|" + row.sra_strategy

            return string

    # If there is no author metadata and they are the same
    if row.sra_strategy == row.pred_strategy:
        return row.sra_strategy

    # EST and RNA-Seq are really similar, so if the SRA uses EST then stick
    # with that.
    if (row.sra_strategy == "EST") & (row.pred_strategy == "RNA-Seq"):
        return "EST"

    # WGA and WGS are really similar, so if the SRA uses WGA then stick
    # with that.
    if (row.sra_strategy == "WGA") & (row.pred_strategy == "WGS"):
        return "WGA"

    # Ambiguous: keep prediction first
    return "|".join([row.pred_strategy, row.sra_strategy])


def annot_library_selection(row: pd.Series) -> str:
    """Aggregates library selection.
    
    Combines library selection from the SRA and predicted values from random
    forest.

    If ambiguous I concatenate values in the following order:

    `predicted|sra`

    The SRA value is always last.

    """
    # If all same
    if row.sra_selection == row.pred_selection:
        return row.sra_selection

    # Ambiguous: keep prediction first
    return "|".join([row.pred_selection, row.sra_selection])


def cross_annotation(row: pd.Series) -> pd.Series:
    """Tweaks based on looking at both strategy and selection."""
    strategy = row.Fear_et_al_library_strategy
    selection = row.Fear_et_al_library_selection

    if ("ChIP-Seq" in strategy) & ~("ChIP" in selection):
        selection = "ChIP|" + selection

    if ~("ChIP-Seq" in strategy) & ("ChIP" in selection):
        strategy = "ChIP-Seq|" + strategy

    _row = row.copy()
    _row.Fear_et_al_library_strategy = strategy
    _row.Fear_et_al_library_selection = selection

    return _row


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input=dict(
                sra="../../output/library_strategy-wf/sra_strategy_selection.parquet",
                free_text="../../output/library_strategy-wf/free_text_library_strategy.parquet",
                forest="../../output/library_strategy-wf/random_forest_predicted_labels.parquet",
            )
        )

    main()
