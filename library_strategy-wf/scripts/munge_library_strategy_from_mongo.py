import pandas as pd
from pymongo import MongoClient


def main():
    """Create table mapping SRX to library strategy"""
    df = pd.DataFrame(
        list(
            NCBI.aggregate(
                [{"$project": {"_id": 0, "srx": 1, "library_strategy": 1, "library_selection": 1}}]
            )
        )
    ).set_index("srx")

    df.to_parquet(snakemake.output[0])

    # Generate table of summary counts
    library_strategy_counts = df.library_strategy.value_counts().to_frame()
    library_strategy_counts.to_csv(snakemake.output[1], sep="\t")


if __name__ == "__main__":
    try:
        MONGO_CLIENT = MongoClient()
        DB = MONGO_CLIENT["sramongo"]
        NCBI = DB["ncbi"]
        main()
    finally:
        MONGO_CLIENT.close()
