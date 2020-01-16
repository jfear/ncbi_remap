import pandas as pd
from pymongo import MongoClient


def main():
    """Create table mapping SRX to library strategy"""
    libstrat = pd.DataFrame(
        list(
            ncbi.aggregate(
                [{"$project": {"_id": 0, "srx": "$srx", "library_strategy": "$library_strategy"}}]
            )
        )
    ).set_index("srx")

    libstrat.to_parquet(snakemake.output[0])

    # Generate table of summary counts
    vcnts = libstrat.library_strategy.value_counts().to_frame()
    vcnts.to_csv(snakemake.output[1], sep="\t")


if __name__ == "__main__":
    mongoClient = MongoClient(host="localhost", port=27017)
    db = mongoClient["sramongo"]
    ncbi = db["ncbi"]

    try:
        main()
    finally:
        mongoClient.close()
