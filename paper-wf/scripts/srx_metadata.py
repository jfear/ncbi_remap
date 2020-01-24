import os
import pandas as pd

from pymongo import MongoClient


def main():
    df = pd.concat([get_sra_metadata(), get_library_strategy()], axis=1, sort=True).rename_axis("srx")
    df.to_csv(snakemake.output[0], sep="\t")


def get_sra_metadata():
    cols = [
        "library_name",
        "BioSample",
        "BioProject",
        "date_created",
        "platform",
        "instrument_model",
    ]
    return (
        pd.DataFrame(
            ncbi.aggregate(
                [
                    {
                        "$project": {
                            "_id": 0,
                            "srx": 1,
                            "library_name": 1,
                            "BioSample": "$BioSample.accn",
                            "BioProject": "$BioProject.accn",
                            "date_created": "$sra_create_date",
                            "platform": 1,
                            "instrument_model": 1,
                        }
                    }
                ]
            )
        )
        .set_index("srx")
        .reindex(columns=cols)
    )


def get_library_strategy():
    return pd.read_parquet(snakemake.input[0])


if __name__ == "__main__":
    if os.getenv("SNAKE_DEBUG"):
        from ncbi_remap.debug import snakemake_debug

        snakemake = snakemake_debug(
            input="../../output/library_strategy-wf/summarized_metadata.parquet"
        )

    try:
        client = MongoClient()
        db = client["sramongo"]
        ncbi = db["ncbi"]
        main()
    finally:
        client.close()
